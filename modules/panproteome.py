import subprocess
import glob
import os
import shutil
import sys
import tempfile
import pandas as pd
from collections import defaultdict
from modules.utils import *
from modules.pfam import Pfam, Hits
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


_PYHMMER_PROFILES = {}


def load_pyhmmer_profiles(pfam_db):
    """
    载入Pfam模型(已hmmpress则直接载入优化后的profile)，同一进程内只载入一次
    """
    if pfam_db not in _PYHMMER_PROFILES:
        from pyhmmer.plan7 import HMMFile
        with HMMFile(pfam_db) as hf:
            if hf.is_pressed():
                _PYHMMER_PROFILES[pfam_db] = list(hf.optimized_profiles())
            else:
                _PYHMMER_PROFILES[pfam_db] = list(hf)
    return _PYHMMER_PROFILES[pfam_db]


class Protein:
    """
    存放蛋白质序列
    """

    def __init__(self, domain: list = None, name: str = "", sequence: str = ""):
        self.name = name
        self.sequence = sequence
        self.domain = domain if domain is not None else []

    def __str__(self):
        return self.name

    def hmm_profile(self, scan_out, scan_type):
        aHits = Hits(scan_out, scan_type)
        self.domain = aHits  # domain是很多pfam的组合



class Proteome(dict):
    """
    存放蛋白质组
    """

    def __init__(self, fasta_name, size: int = None):
        super().__init__()
        self.name = os.path.basename(fasta_name).replace('.faa', '')
        self.size = size
        self.file = fasta_name  # 带路径的file

    def __str__(self):
        return self.name

    def add_seq(self):
        fasta_seqs = gen_seqs_with_headers(self.file)
        for seqid, seq in fasta_seqs.items():
            seqid = seqid.split(" ")[0]  # 去除faa文件里的功能描述部分
            # for protein in self:  # protein是Protein对象
            #     if protein.name == seqid:
            #         protein.sequence = seq
            # # 如果没有找到匹配的Protein，创建新的Protein并添加到proteins列表中
            # new_protein = Protein(name=seqid, sequence=seq)
            # self.append(new_protein)
            if seqid in self:
                # 更新已存在的Protein对象
                self[seqid].sequence = seq
            else:
                # 创建并添加新的Protein对象
                new_protein = Protein(name=seqid, sequence=seq)
                self[seqid] = new_protein


    def read_fasta(self, fasta_name, extract_ids):
        fasta_file = gen_seqs_with_headers(fasta_name, extract_ids)
        self.name = os.path.basename(fasta_name).replace('.faa', '')
        if not extract_ids: # 提取序列
            for seqid, seq in fasta_file.items():
                seqid = seqid.split(" ")[0]  # 去除faa文件里的功能描述部分
                validate_fasta_ids(seqid, self.name)
                aProtein = Protein(name=seqid, sequence=seq)
                self[seqid] = aProtein
        elif extract_ids: # 仅提取ID
            for seqid in fasta_file:
                seqid = seqid.split(" ")[0]
                validate_fasta_ids(seqid, self.name)
                aProtein = Protein(name=seqid)
                self[seqid] = aProtein
        self.size = len(self)

    def search_mmseq_pfam_domain(self, mmseq_db, out_dir, threads, evalue=1e-5):
        """
        使用mmseqs注释pfam
        """
        res = os.path.join(out_dir, f'result_{self.name}.pfam')
        tmp_dir = os.path.join(out_dir, f'{str(self.name)}_tmp')
        mmseqs_cmd = ' '.join(
            ['mmseqs', 'easy-search', self.file, mmseq_db, res, tmp_dir,
             '-e', str(evalue), '--threads', str(threads),
             '--format-output', 'query,target,qlen,qstart,qend,tlen,tstart,tend,alnlen,bits,evalue,gapopen,fident'])
        proc = subprocess.run(mmseqs_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f'mmseqs easy-search failed for {self.name}: {proc.stderr[-2000:]}')

        if os.path.exists(res):
            # 处理mmseqs的结果
            necessary_columns = ['query', 'target', 'qlen', 'qstart', 'qend']
            grouped_data = pd.read_csv(res, sep='\t', header=None,
                                       names=['query', 'target', 'qlen', 'qstart', 'qend',
                                              'tlen', 'tstart', 'tend', 'alnlen', 'bits',
                                              'evalue', 'gapopen', 'fident'],
                                       usecols=necessary_columns).groupby('query')

            for query, group in grouped_data:
                aProtein = Protein(name=str(query))
                aProtein.hmm_profile(group, 'mmseqs')
                self[query] = aProtein
            os.remove(res)
            shutil.rmtree(tmp_dir)
        return

    def search_pfam_domain_batch(self, threads, pfam_db, evalue='1e-5', chunks=None, cpu=1):
        """
        每个蛋白质组只做一次批量hmmscan，不再逐条蛋白启动进程。
        hmmscan对每条查询序列都要遍历整个数据库且多线程扩展性差，
        所以把序列分成threads块，每块一个hmmscan进程(--cpu 1)并行运行，再按蛋白合并domtblout结果。
        """
        names = list(self.keys())
        cpu = max(1, int(cpu))
        if chunks is None:
            chunks = max(1, int(threads) // cpu)  # 默认块数 = 总线程数 / 每块cpu
        n_chunks = max(1, min(int(chunks), len(names)))
        chunks = [names[i::n_chunks] for i in range(n_chunks)]  # 轮转分块，使各块负载均衡

        def run_chunk(idx, chunk_names, tmp_dir):
            in_fa = os.path.join(tmp_dir, f'in_{idx}.faa')
            domtbl = os.path.join(tmp_dir, f'out_{idx}.domtbl')
            with open(in_fa, 'w') as f:
                for name in chunk_names:
                    f.write(f'>{name}\n{self[name].sequence}\n')
            cmd = ['hmmscan', '--cpu', str(cpu), '-E', str(evalue), '--domE', str(evalue),
                   '--noali', '-o', os.devnull, '--domtblout', domtbl, pfam_db, in_fa]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0:
                raise RuntimeError(f'hmmscan failed for {self.name} (chunk {idx}): {res.stderr}')
            grouped_chunk = defaultdict(list)
            with open(domtbl) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    grouped_chunk[line.split()[3]].append(line)
            return grouped_chunk

        grouped = {}
        with tempfile.TemporaryDirectory(prefix='hmmscan_') as tmp_dir:
            with ThreadPoolExecutor(max_workers=n_chunks) as executor:
                futures = [executor.submit(run_chunk, i, c, tmp_dir) for i, c in enumerate(chunks) if c]
                for future in futures:
                    grouped.update(future.result())  # 任一块失败都会抛出异常
        for name, protein in self.items():
            protein.hmm_profile(grouped.get(name, []), scan_type='hmmscan')
        return

    def search_pfam_domain_pyhmmer(self, threads, pfam_db, evalue='1e-5'):
        """
        使用pyhmmer在内存中做hmmscan，Pfam模型只载入一次并缓存，不再重复读库和启动子进程
        """
        import pyhmmer
        profiles = load_pyhmmer_profiles(pfam_db)
        alphabet = pyhmmer.easel.Alphabet.amino()
        names = list(self.keys())
        queries = []
        for name in names:
            seq = pyhmmer.easel.TextSequence(name=name.encode(), sequence=self[name].sequence)
            queries.append(seq.digitize(alphabet))
        e = float(evalue)
        results = {}
        for hits in pyhmmer.hmmer.hmmscan(queries, profiles, cpus=int(threads), E=e, domE=e):
            qname = hits.query.name
            qname = qname.decode() if isinstance(qname, bytes) else qname
            qlen = len(hits.query)
            rows = []
            for hit in hits:
                if not hit.reported:
                    continue
                acc = hit.accession
                acc = acc.decode() if isinstance(acc, bytes) else acc
                for domain in hit.domains:
                    if not domain.reported:
                        continue
                    rows.append((acc, domain.alignment.target_from, domain.alignment.target_to, qlen))
            results[qname] = rows
        for name, protein in self.items():
            protein.hmm_profile(results.get(name, []), scan_type='pyhmmer')
        return

    def write_out_pfam(self, out_dir):  # self是一个faa文件
        proteome_domain = defaultdict(lambda: defaultdict(list))
        for protein in self.values():  # 每个orf
            domain_list = protein.domain
            if domain_list:
                for domain_o in domain_list:
                    proteome_domain[protein.name]["Domain"].append(domain_o.id)
                    proteome_domain[protein.name]["LenCov"].append(domain_o.percent)
            else:
                proteome_domain[protein.name]["Domain"] = list()
                proteome_domain[protein.name]["LenCov"] = list()
            protein.sequence = ""  # 每个蛋白质做完后就把序列的属性删掉
        FileOperator(f'{self.name}.pfam', out_dir, "json", proteome_domain).write()


class Panproteome(list):
    """
    存放泛基因组
    """

    def __init__(self, f, outdir, threads, db, method, evalue, hmm_chunks=None, hmm_cpu=1,
                 hmm_engine='pyhmmer'):
        super().__init__()

        # 获取所有需要处理的faa文件列表
        faa_files = sorted(glob.glob(os.path.join(f, '*.faa')))
        empty_files = [f for f in faa_files if os.path.getsize(f) == 0]
        if empty_files:
            message(text="Empty FASTA files detected", label="Error")
            for f in empty_files:
                message(text=f"- {os.path.basename(f)}", label="FILE")
            sys.exit(1)

        # 判断已经完成的faa
        completed_proteomes = [os.path.splitext(os.path.basename(file))[0] for file in
                               glob.glob(os.path.join(outdir, '*.pfam'))]
        message(text=f"{len(completed_proteomes)} already processed. Skipping...", label='Information')

        # 创建进度条
        with tqdm(total=len(faa_files), desc="Processing proteomes", unit="proteome") as pbar:
            for faa_file in faa_files:
                aProteome = Proteome(fasta_name=faa_file)
                self.append(aProteome)

                if aProteome.name in completed_proteomes:
                    # 已处理过的逻辑
                    aProteome.read_fasta(fasta_name=faa_file, extract_ids=True)
                    json_data = FileOperator(f'{aProteome}.pfam', outdir, "json")
                    json_data.read()
                    for protein in aProteome.values():
                        protein: Protein
                        domain_list = []
                        for i in range(len(json_data.data[protein.name]['Domain'])):
                            pfam_id = json_data.data[protein.name]['Domain'][i]
                            percent = json_data.data[protein.name]['LenCov'][i]
                            domain_list.append(Pfam(pfam_id=pfam_id, percent=percent))
                        protein.domain = domain_list
                    pbar.set_postfix_str(f"Skipped {aProteome.name}", refresh=False)
                else:
                    if method == 'hmmscan':
                        aProteome.read_fasta(fasta_name=faa_file, extract_ids=False)
                        if hmm_engine == 'pyhmmer':
                            aProteome.search_pfam_domain_pyhmmer(threads=threads, pfam_db=db, evalue=evalue)
                        else:
                            aProteome.search_pfam_domain_batch(threads=threads, pfam_db=db, evalue=evalue,
                                                               chunks=hmm_chunks, cpu=hmm_cpu)
                    elif method == 'mmseqs-search':
                        aProteome.read_fasta(fasta_name=faa_file, extract_ids=True)
                        aProteome.search_mmseq_pfam_domain(mmseq_db=db, out_dir=outdir,
                                                           threads=threads, evalue=evalue)
                    aProteome.write_out_pfam(out_dir=outdir)
                    pbar.set_postfix_str(f"Processed {aProteome.name}", refresh=False)

                # 更新进度条
                pbar.update(1)


    def add_proteome_sequence(self):
        for proteome in self:
            proteome.add_seq()

    def remove_redundant_sequences(self, outdir):
        pfam_dict = {}
        for proteome in self:
            pfam_set = set()
            for protein in proteome:
                if sum([i.percent for i in protein.domain]) <= 0.6:  # 在长度上判断
                    pfam_ = "None"
                else:
                    pfam_ = ','.join(sorted([i.id for i in protein.domain]))
                pfam_set.add(pfam_)
            # 如果这个pfam不存在set中
            pfam_key = frozenset(pfam_set)
            if pfam_key not in pfam_dict:
                pfam_dict[pfam_key] = []
            pfam_dict[pfam_key].append(proteome)
        # 随机选择每个pfam组中的一个对象
        unique_proteomes = []
        redundancy_string = ''
        for group in pfam_dict.values():
            # 按照 size 降序排列，如果 size 相同则按名称升序排列
            sorted_group = sorted(group, key=lambda x: (-x.size, x.name))
            represent = sorted_group[0]
            unique_proteomes.append(represent)
            redundancy_string += f"* {represent.name}\n" + ' '.join([i.name for i in group]) + '\n'
        FileOperator('redundancy_infomation.txt', dir_=outdir, data=redundancy_string).write()
        self.clear()
        self.extend(unique_proteomes)
        return
