
import os
from multiprocessing import Pool
from pyfasta import Fasta


class Proteome(Fasta):

    def __init__(self, fasta_name):
        super().__init__(fasta_name)

    def do_hmmscan(self):
        return
