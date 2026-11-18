"""pyhmmer 与 hmmscan 批量版注释结果一致性测试: python tests/test_hmm_consistency.py [N条蛋白]"""
import os, sys, itertools, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.panproteome import Proteome

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(ROOT, 'modules', 'Pfam-A.hmm')
FAA = os.path.join(ROOT, 'testdata', 'Staphylococcus', 'GCA_000236925.1.faa')


def load(n, tmp):
    sub = os.path.join(tmp, 'GCA_000236925.1.faa')
    seqs, name = {}, None
    with open(FAA) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if len(seqs) >= n:
                    break
                name = line[1:].split()[0]; seqs[name] = ''
            elif name:
                seqs[name] += line
    with open(sub, 'w') as o:
        for k, v in seqs.items():
            o.write(f'>{k}\n{v}\n')
    p = Proteome(sub); p.read_fasta(sub, extract_ids=False)
    return p


def result(p):
    return {k: [(d.id, round(d.percent, 6)) for d in v.domain] for k, v in p.items()}


def main(n=30):
    if not os.path.exists(DB):
        print('SKIP: Pfam-A.hmm not found'); return 0
    with tempfile.TemporaryDirectory() as tmp:
        a = load(n, tmp); a.search_pfam_domain_pyhmmer(2, DB)
        b = load(n, tmp); b.search_pfam_domain_batch(2, DB, chunks=2, cpu=1)
        ra, rb = result(a), result(b)
    ok = ra == rb
    print('IDENTICAL' if ok else 'DIFFERENT', f'({n} proteins, {sum(1 for v in ra.values() if v)} with domain)')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 30))
