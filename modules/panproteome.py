import os
import glob
from multiprocessing import Pool, Manager

from pyfasta import Fasta

from modules.pfam import *
from modules.utils import *
from build_graph import *

pfamDB = os.path.join(os.getcwd(), 'modules', 'Pfam-A.hmm')


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

    def _hmm_scan(self, evalue):
        in_temp = make_temp_file(prefix='in_', close=False)
        in_temp.write(f'>{self.name}\n{self.sequence}')
        in_temp.close()
        out_temp = make_temp_file(prefix='out_', close=True)
        cmd = ["hmmscan", "-E", evalue, "--domE", evalue, "--domtblout", out_temp.name, pfamDB, in_temp.name]
        # cmd2 = ['hmmscan', '--cut_ga', '--domtblout', out_temp.name, pfamDB, in_temp.name]
        cap = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cap.wait()
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
        self.name = os.path.basename(fasta_name).split('/')[-1]
        self.size = len(self)

    def search_pfam_domain(self, treads=60, evalue='1e-5'):
        """
        搜索Pfam-A.hmm
        :return:
        """
        processes = Pool(processes=treads)
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


class Panproteome(list):

    def __init__(self, f):
        super().__init__()
        self.member = glob.glob(os.path.join(f, '*.faa'))
        self.proteomes = []
        self._load_proteomes()
        self.domains = []
        self._get_all_domains()

    def _load_proteomes(self):  # 实例化为Proteome对象
        for file in self.member:
            aProteome = Proteome(file)
            self.proteomes.append(aProteome)

    def _get_all_domains(self):
        for proteome in self.proteomes:
            proteome.search_pfam_domain()

        all_domains = []
        for proteome in self.proteomes:
            for protein in proteome:
                all_domains.append(protein.domain)
        self.domains = all_domains

    def graph(self):
        vs_es = get_edges(self.domains)
        vs = vs_es[0]
        es = vs_es[1]
        g = build_graph(vs, es)
        return g
