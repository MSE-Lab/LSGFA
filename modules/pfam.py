import os


class Pfam:

    def __init__(self, name, id, length, start, end, hitlength):
        self.name = name
        self.id = id
        self.length = length
        self.start = start
        self.end = end
        self.hitlength = hitlength

    def relation(self, other):
        # 用于计算两个pfam之间的关系
        return


class Hits(list):

    def __init__(self, hmmscan_out):
        super().__init__()
        with open(hmmscan_out) as f:
            content = f.readlines()
            hits = [i for i in content if not i.startswith("#")]
            self.raw = hits
        if self.raw:
            self._parse()

    def _parse(self):
        for i in self.raw:
            ilist = [_ for _ in i.split(' ') if _ != '']
            pfam_name = ilist[0]
            pfam_id = ilist[1]
            pfam_length = ilist[2]
            pfam_start = ilist[17]
            pfam_end = ilist[18]
            pfam_hit_length = int(pfam_end)-int(pfam_start)
            aPfam = Pfam(name=pfam_name, id=pfam_id, length=pfam_length,
                         start=pfam_start, end=pfam_end, hitlength=pfam_hit_length)
            self.append(aPfam)

    def ana_relations(self):
        # 用于判断这个Hit里要保留哪个pfam

        pass