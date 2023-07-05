import os
import glob
import subprocess
from multiprocessing import Pool, Manager
from pyfasta import Fasta
from tempfile import NamedTemporaryFile

pfamDB = os.path.join('/Users/zhixiaoyang/PycharmProjects/LSGFA/modules/Pfam-A.hmm')


class Protein:

    def __init__(self, name: str = "", sequence: str = ""):
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return self.name


class Proteome(list):

    def __init__(self, name: str = None, size: int = None):
        super().__init__()
        self.name = name
        self.size = size

    def __str__(self):
        return self.name

    def read_fasta(self, fasta_name):
        self.name = os.path.basename(fasta_name).split('/')[-1]
        fasta_file = Fasta(fasta_name)
        for seqid, seq in fasta_file.items():
            aProtein = Protein(name=seqid, sequence=seq)
            self.append(aProtein)
        self.size = len(self)

    @staticmethod
    def hmmscan(protein, evalue):
        in_temp = NamedTemporaryFile(mode='w', delete=False)
        in_temp.write(f'>{protein.name}\n{protein.sequence}')
        in_temp.close()
        cmd = ["hmmscan", "-E", "1e-5", "--domE", "1e-5", "--domtblout", "out1", pfamDB, in_temp.name]
        cap = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cap.wait()
        print(cap.returncode)
        for e in cap.stderr:
            print(e)

    @staticmethod
    def get_results(queue):
        results = list()
        return results

    def search_domain(self, treads=60, evalue='1e-5'):
        """
        搜索Pfam-A.hmm
        :return:
        """
        processes = Pool(processes=treads)
        aQueue = Manager().Queue()
        for protein in self:
            self.hmmscan(protein, evalue)
            break
            # processes.apply_async(self.hmmscan, args=(aQueue, protein, evalue))
        processes.close()
        result = list()
        # ntd = 0
        # while True:
        #     a = aQueue.get(timeout=None)
        #     result.append(a)
        #     print(a)
        #     ntd += 1
        #     if ntd == self.size:
        #         break
        return result


class Panproteome(list):

    def __init__(self, f):
        super().__init__()
        self.member = glob.glob(os.path.join(f, '*.faa'))

    def func(self):
        pass
