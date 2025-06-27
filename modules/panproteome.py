import subprocess
import glob
import os
import shutil
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, Manager
from modules.utils import *
from modules.pfam import Pfam, Hits
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


class Protein:
    """
    存放蛋白质序列
    """

    def __init__(self, domain: list = [], name: str = "", sequence: str = ""):
        self.name = name
        self.sequence = sequence
        self.domain = domain

    def __str__(self):
        return self.name

    def hmm_profile(self, scan_out, scan_type):
        aHits = Hits(scan_out, scan_type)
        self.domain = aHits  # domain是很多pfam的组合

    def _hmm_scan(self, pfamDB, evalue='1e-5'):
        """
        对每条蛋白序列使用hmmscan进行Pfam注释
        :param evalue: evalue阈值
        """
        in_temp = make_temp_file(prefix='in_', close=False)
        in_temp.write(f'>{self.name}\n{self.sequence}')
        in_temp.close()
        out_temp = make_temp_file(prefix='out_', close=True)
        cmd = ["hmmscan", "-E", str(evalue), "--domE", str(evalue), "--domtblout", out_temp.name, pfamDB, in_temp.name]
        # cmd2 = ['hmmscan', '--cut_ga', '--domtblout', out_temp.name, pfamDB, in_temp.name]
        cap = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cap.communicate()
        if cap.returncode != 0:
            # some errors happened
            for e in cap.stderr:
                print(e)
        else:
            # remove input temporary sequence file
            os.unlink(in_temp.name)
        self.hmm_profile(out_temp.name, scan_type='hmmscan')
        return out_temp.name

    def identify_pfam(self, pfamDB, evalue='1e-5'):
        out_hmmscan = self._hmm_scan(pfamDB, evalue)
        os.unlink(out_hmmscan)
        # 通过进程池 Protein.domain 似乎无法更新，所以把结果put出来在进程池外再更新
        # queue.put((self.name, self.domain))


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
        mmseqs_cmd = subprocess.Popen(mmseqs_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        mmseqs_cmd.communicate()

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

    def search_pfam_domain(self, threads, pfam_db, evalue='1e-5'):
        """
        搜索Pfam-A.hmm
        """
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # 提交所有任务到线程池，并创建一个进度条
            futures = {executor.submit(protein.identify_pfam, pfam_db, evalue):
                           protein for protein in self.values()}
            # 使用 tqdm 创建进度条
            for future in as_completed(futures):
                try:
                    result = future.result()  # This will raise an exception if the function raised one
                except Exception as e:
                    print(f"Error processing {futures[future]}: {e}")
            # 等待所有任务完成（虽然 as_completed 已经做了这个，但我们可以显式地处理）
            for future in futures.keys():
                future.result()

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

    def __init__(self, f, outdir, threads, db, method, evalue):
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
                        aProteome.search_pfam_domain(threads=threads, pfam_db=db, evalue=evalue)
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
