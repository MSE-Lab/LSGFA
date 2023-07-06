import os


class Pfam:

    def __init__(self, name, id, length):
        self.name = name
        self.id = id
        self.length = length

    def relation(self, other):
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
            ilist = [_ for _ in i.split(' ') if _ is not '']
            pfam_name = ilist[0]
            pfam_id = ilist[1]
            pfam_length = ilist[3]
            aPfam = Pfam(name=pfam_name, id=pfam_id, length=pfam_length)
            self.append(aPfam)

    def ana_relations(self):
        pass