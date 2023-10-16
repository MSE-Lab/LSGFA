import json

from modules.panproteome import *

if __name__ == '__main__':
    start = time.time()
    faas = './testdata'
    pp = Panproteome(faas)
    pp.make_pfam_graph(outdir=faas)
    end = time.time()
    print('total time ',end - start)


a = './testdata_3/GCA_001751275.1.pfam'
with open(a, 'r') as f:
    json_data = json.load(f)