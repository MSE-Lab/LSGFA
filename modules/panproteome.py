import subprocess
import glob
import os
import random
from collections import defaultdict
from multiprocessing import Pool, Manager
from pyfasta import Fasta
from modules.utils import *
from modules.pfam import *
from tqdm import tqdm

# pfamDB = os.path.join(os.getcwd(), 'modules', 'database', 'Pfam-A.hmm')
pfamDB = '/media/disk2/biodatabases/Pfam/Pfam-A.hmm'  # 浪潮
# pfamDB = '/home/biodbs/Pfam35.0/Pfam-A.hmm'  # 集群


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

    def _hmm_profile(self, scan_out):
        aHits = Hits(scan_out)
        self.domain = aHits  # domain是很多pfam的组合

    def _hmm_scan(self, evalue='1e-5'):
        """
        对每条蛋白序列使用hmmscan进行Pfam注释
        :param evalue: evalue阈值
        """
        in_temp = make_temp_file(prefix='in_', close=False)
        in_temp.write(f'>{self.name}\n{self.sequence}')
        in_temp.close()
        out_temp = make_temp_file(prefix='out_', close=True)
        cmd = ["hmmscan", "-E", evalue, "--domE", evalue, "--domtblout", out_temp.name, pfamDB, in_temp.name]
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
        self._hmm_profile(out_temp.name)
        return out_temp.name

    def identify_pfam(self, queue, evalue='1e-5'):
        out_hmmscan = self._hmm_scan(evalue)
        os.unlink(out_hmmscan)
        # 通过进程池 Protein.domain 似乎无法更新，所以把结果put出来在进程池外再更新
        queue.put((self.name, self.domain))


class Proteome(list):
    """
    存放蛋白质组
    """

    def __init__(self, fasta_name, name: str = None, size: int = None):
        super().__init__()
        self.name = name
        self.size = size
        self._read_fasta(fasta_name)

    def __str__(self):
        return self.name

    def _read_fasta(self, fasta_name):
        fasta_file = Fasta(fasta_name)
        for seqid, seq in fasta_file.items():
            seqid = seqid.split(" ")[0]  # 去除faa文件里的功能描述部分
            aProtein = Protein(name=seqid, sequence=seq)
            self.append(aProtein)
        self.name = os.path.basename(fasta_name).replace('.faa', '')
        self.size = len(self)
        os.remove(f"{fasta_name}.flat")
        os.remove(f"{fasta_name}.gdx")

    def search_pfam_domain(self, threads=60, evalue='1e-5'):
        """
        搜索Pfam-A.hmm
        """
        processes = Pool(processes=threads)
        aQueue = Manager().Queue()
        for protein in self:
            processes.apply_async(protein.identify_pfam, args=(aQueue, evalue))
        processes.close()
        identified_pfams = dict()
        ntd = 0
        while True:
            a = aQueue.get(timeout=None)
            identified_pfams[a[0]] = a[1]
            ntd += 1
            if ntd == self.size:
                # 当从进程池队列中get的结果数量等于Proteome的size，即蛋白质序列数量时停止get结果
                break
        for protein in self:
            protein.domain = identified_pfams[protein.name]
        return

    def write_out_pfam(self, out_dir):  # self是一个faa文件
        proteome_domain = defaultdict(lambda: defaultdict(list))
        for protein in self:  # 每个orf
            domain_list = protein.domain
            if domain_list:
                for domain_o in domain_list:
                    proteome_domain[protein.name]["Domain"].append(domain_o.id)
                    proteome_domain[protein.name]["LenCov"].append(domain_o.percent)
            else:
                proteome_domain[protein.name]["Domain"] = list()
                proteome_domain[protein.name]["LenCov"] = list()
        FileOperator(f'{self.name}.pfam', out_dir, "json", proteome_domain).write()


class Panproteome(list):
    """
    存放泛基因组
    """

    def __init__(self, f):
        super().__init__()
        for faa_file in sorted(glob.glob(os.path.join(f, '*.faa'))):
            aProteome = Proteome(fasta_name=faa_file)
            self.append(aProteome)

    def _identify_pfam(self, threads, outdir):
        """
        对每个基因组进行Pfam的鉴定
        存在鉴定结果的，读取鉴定结果，进行Protein实例化
        :param threads: 线程
        :param outdir: 输出目录
        """

        proteome: Proteome
        completed_proteomes = [os.path.splitext(os.path.basename(file))[0] for file in
                      glob.glob(os.path.join(outdir, '*.pfam'))]
        message(text=f"{len(completed_proteomes)} already processed. Skipping...", label='Information')
        for proteome in self:
            if proteome.name in completed_proteomes:
                json_data = FileOperator(f'{proteome}.pfam', outdir, "json")
                json_data.read()
                for protein in proteome:
                    protein: Protein
                    domain_list = []
                    for i in range(len(json_data.data[protein.name]['Domain'])):
                        pfam_id = json_data.data[protein.name]['Domain'][i]
                        percent = json_data.data[protein.name]['LenCov'][i]
                        domain_list.append(Pfam(pfam_id=pfam_id, percent=percent))
                    protein.domain = domain_list
            else:
                message(text=f'identify Pfam for {proteome.name}')
                proteome.search_pfam_domain(threads=threads)
                proteome.write_out_pfam(out_dir=outdir)

    def put_pfam_file(self, threads, outdir):
        # 并行的运行hmmscan为每个proteome.faa鉴定pfam
        self._identify_pfam(threads=threads, outdir=outdir)

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
        FileOperator('redundancy_infomation.txt', dir_ = outdir, data=redundancy_string).write()
        return unique_proteomes
