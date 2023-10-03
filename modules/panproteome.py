import os
import glob
import shutil
import json
from multiprocessing import Pool, Manager
# import pandas as pd

from pyfasta import Fasta

from modules.pfam import *
from modules.utils import *
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

    def _write_out_pfam(self, outdir):  # self是一个faa文件
        with open(os.path.join(os.getcwd(), outdir, f'{self.name}.pfam'), 'w') as out:
            proteome_dic = {}
            for protein in self:  # 每个orf
                try:
                    domains = {}
                    for d in protein.domain:
                        if d.id in domains:
                            p_ = domains[d.id]
                            domains[d.id] = p_ + d.percent
                        else:
                            domains[d.id] = d.percent
                except TypeError:
                    domains = 'None'
                protein_dict = {protein.name: domains}  # {'orf1': {'pf1': 0.5, 'pf2': 0.3}}
                proteome_dic.update(protein_dict)
            proteome_json = json.dumps(proteome_dic, sort_keys=True, indent=4, separators=(',', ': '))
            out.write(proteome_json)


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
                proteome._write_out_pfam(outdir=outdir)

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

    def make_pfam_graph(self, outdir, intersect_t):

        members = glob.glob(os.path.join(os.getcwd(), outdir, '*.pfam'))
        vs_es = get_edges(members, intersect_t)
        vs = vs_es[0]
        es = vs_es[1]

        g = build_graph(vs, es)
        partitions = graph_split(g)
        # 得到的partitions是一个嵌套列表，其中每个列表表示一个clade，每个clade里包含基因的编号
        return partitions
