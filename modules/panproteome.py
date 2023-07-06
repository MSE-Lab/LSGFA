import os
import glob
import subprocess
from multiprocessing import Pool, Manager
from pyfasta import Fasta

from modules.pfam import *
from modules.utils import *

pfamDB = os.path.join(os.getcwd(), 'modules', 'Pfam-A.hmm')


class Protein:

    def __init__(self, name: str = "", sequence: str = "", domain: Pfam = None):
        self.name = name
        self.sequence = sequence
        self.domain = domain

    def __str__(self):
        return self.name

    def _hmmscan(self, evalue):
        in_temp = make_temp_file(prefix='in_', close=False)
        in_temp.write(f'>{self.name}\n{self.sequence}')
        in_temp.close()
        out_temp = make_temp_file(prefix='out_', close=True)
        cmd = ["hmmscan", "-E", evalue, "--domE", evalue, "--domtblout", out_temp.name, pfamDB, in_temp.name]
        cap = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cap.wait()
        if cap.returncode != 0:
            # some errors happened
            for e in cap.stderr:
                print(e)
        else:
            # remove input temporary sequence file
            os.unlink(in_temp.name)
        return out_temp.name

    @staticmethod
    def _domain(hmmscan_result):
        return

    def pfam(self, queue, evalue='1e-5'):
        out_hmmscan = self._hmmscan(evalue)
        domain = self._domain(out_hmmscan)
        self.domain = domain
        queue.put(out_hmmscan)


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
            processes.apply_async(protein.pfam, args=(aQueue, evalue))
        processes.close()
        result = list()
        ntd = 0
        while True:
            a = aQueue.get(timeout=None)
            result.append(a)
            print(a)
            ntd += 1
            if ntd == self.size:
                break
        return result


class Panproteome(list):

    def __init__(self, f):
        super().__init__()
        self.member = glob.glob(os.path.join(f, '*.faa'))

    def func(self):
        pass
