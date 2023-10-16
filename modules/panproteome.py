import glob
import os.path
import shutil
import subprocess
from multiprocessing import Pool, Manager
from collections import defaultdict

import igraph

from pyfasta import Fasta

from modules.pfam import *
from modules.build_graph import *

# pfamDB = os.path.join(os.getcwd(), 'modules', 'Pfam-A.hmm')
pfamDB = '/media/disk2/biodatabases/Pfam/Pfam-A.hmm'


class Protein:

    def __init__(self, name: str = "", sequence: str = "", domain: list = None):
        self.name = name
        self.sequence = sequence
        self.domain = domain

    def __str__(self):
        return self.name

    def _hmm_profile(self, scan_out):
        aHits = Hits(scan_out)
        hits_clean = aHits.ana_relations()

        self.domain = hits_clean

    def _hmm_scan(self, evalue='1e-5'):
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
        # 通过进程池 Protein.domain 似乎无法更新，所以吧结果put出来在进程池外再更新
        queue.put((self.name, self.domain))


class Proteome(list):

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
            aProtein = Protein(name=seqid, sequence=seq)
            self.append(aProtein)
        self.name = os.path.basename(fasta_name).replace('.faa', '')
        self.size = len(self)
        # os.remove(f"{fasta_name}.flat")
        # os.remove(f"{fasta_name}.gdx")

    def search_pfam_domain(self, threads=60, evalue='1e-5'):
        """
        搜索Pfam-A.hmm
        :return:
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

    def _write_out_pfam(self, out_dir):  # self是一个faa文件
        # 用于输出pfam的文件
        proteome_domain = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for protein in self:  # 每个orf
            domain_list = protein.domain  # 提取蛋白的domain
            if domain_list is not None:
                combined = ";".join(sorted([d.id for d in domain_list]))  # 组合
                for domain_o in domain_list:
                    proteome_domain[combined][protein.name]["Domain"].append(domain_o.id)
                    proteome_domain[combined][protein.name]["LenCov"].append(domain_o.percent)
            else:  # 注释不到pfam时就返回空列表
                proteome_domain["None"][protein.name]["Domain"] = list()
                proteome_domain["None"][protein.name]["LenCov"] = list()
        FileOperator(f'{self.name}.pfam', out_dir, "json", proteome_domain).write()


class Panproteome(list):
    completed_file = "completed_proteomes.txt"

    def __init__(self, f):
        super().__init__()
        for faa_file in glob.glob(os.path.join(f, '*.faa')):
            aProteome = Proteome(fasta_name=faa_file)
            self.append(aProteome)

    def _identify_pfam(self, threads, outdir):

        proteome: Proteome
        try:
            with open(self.completed_file, "r") as file:
                completed_proteomes = file.read().splitlines()
        except FileNotFoundError:
            completed_proteomes = []

        for proteome in self:
            if proteome.name in completed_proteomes:
                print(f"Proteome {proteome.name} already processed. Skipping...")
            else:
                message(text=f'identify Pfam for {proteome.name}')
                proteome.search_pfam_domain(threads=threads)
                proteome._write_out_pfam(out_dir=outdir)
                # 将已完成的proteome_id记录到文件中
                with open(self.completed_file, "a") as file:
                    file.write(f"{proteome.name}\n")

    @staticmethod
    def backup_completed_file():
        if os.path.exists(Panproteome.completed_file):
            shutil.copy(Panproteome.completed_file, "completed_backup.txt")
            os.remove(Panproteome.completed_file)

    def put_pfam_file(self, threads, outdir):
        # 并行的运行hmmscan为每个proteome.faa鉴定pfam
        self._identify_pfam(threads=threads, outdir=outdir)

    @staticmethod
    def make_pfam_graph(outdir, d_co_index=0.5, d_co_len=0.8):
        # 作图
        pfam = PGraph(outdir)
        pfam.get_full_connected_edges()
        pfam.get_append_edges(d_co_index, d_co_len)
        vs = pfam.genes
        es = pfam.related_edges

        g = igraph.Graph()
        g.add_vertices(vs)
        g.add_edges(es)
        print(g.ecount())
        g.write_gml(os.path.join(outdir,'pfam_graph.gml'))
        return g
