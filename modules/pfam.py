
import os
from tempfile import TemporaryDirectory

from multiprocessing import Pool
from pyfasta import Fasta


class Hit:

    def __init__(self):
        pass

    def compare(self, other):
        pass


class Protein:

    def __init__(self):
        pass


class Proteome(Fasta):

    def __init__(self, fasta_name):
        super().__init__(fasta_name)

    def do_hmmscan(self):
        return
