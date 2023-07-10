import os
from itertools import combinations


class Pfam:

    def __init__(self, name, pfid, length, orf, start, end, hitlength):
        self.name = name
        self.id = pfid
        self.length = length
        self.orf = orf
        self.start = start
        self.end = end
        self.hitlength = hitlength

    def relation(self, other):
        # 处理两个Pfam之间的关系
        if self.hitlength >= other.hitlength:
            long = self
            short = other
        else:
            long = other
            short = self
        if long.start >= short.end or long.end <= short.start:  # 独立的两个
            return True
        else:  # 当两个为包含、相等或overlap时
            return False


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
            pfam_length = int(ilist[2])
            pfam_orf = ilist[3]
            pfam_start = int(ilist[17])
            pfam_end = int(ilist[18])
            pfam_hit_length = pfam_end-pfam_start
            aPfam = Pfam(name=pfam_name, pfid=pfam_id, length=pfam_length, orf=pfam_orf,
                         start=pfam_start, end=pfam_end, hitlength=pfam_hit_length)
            self.append(aPfam)

    def ana_relations(self):
        """
        处理一个Hit里的所有Pfam，保留相互独立的pfam，且每个pfam都是最长的
        :return: 一个list，里面的元素为Pfam对象
        """
        # 对列表按长度进行降序排列
        sorted_list = sorted(self, key=lambda x: x.hitlength, reverse=True)
        cleaned_list = [sorted_list[0]]
        for i in range(1, len(sorted_list)):  # 从第二个开始比较
            current_pfam = sorted_list[i]
            independent = True
            for j in range(len(cleaned_list)):  # 当前的pfam与clean_list中的pfam对比
                existing_pfam = cleaned_list[j]
                if not Pfam.relation(current_pfam, existing_pfam):  # 与已有的Pfam对象不独立
                    if current_pfam.hitlength > existing_pfam.hitlength:
                        cleaned_list[j] = current_pfam  # 替换为较长的Pfam对象
                    else:
                        independent = False
                    break
            if independent:
                cleaned_list.append(current_pfam)
        return cleaned_list
