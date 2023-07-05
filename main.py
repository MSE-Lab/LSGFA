from modules.panproteome import *


if __name__ == '__main__':
    p = Proteome()
    p.read_fasta(fasta_name="./testdata/A.faa")

    f = p.search_domain(treads=1)
