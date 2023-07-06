from modules.panproteome import *


if __name__ == '__main__':
    p = Proteome(fasta_name="./testdata/B.faa")

    f = p.search_pfam_domain(treads=5)
