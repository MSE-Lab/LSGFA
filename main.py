from modules.panproteome import *

if __name__ == '__main__':
    faas = './testdata'
    pp = Panproteome(faas)
    pp.make_pfam_graph(threads=60, outdir=faas)
