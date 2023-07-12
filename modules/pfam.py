class Pfam:

    def __init__(self, pfam_name, pfam_id, length, orf, start, end, hit_len):
        self.name = pfam_name
        self.id = pfam_id
        self.length = length
        self.orf = orf
        self.start = start
        self.end = end
        self.hit_len = hit_len

    def __str__(self):
        return self.name

    def related_to(self, other):
        # 处理两个Pfam之间的关系
        if self.hit_len >= other.hit_len:
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
            split_line = [_ for _ in i.split(' ') if _ != '']
            pfam_name = split_line[0]
            pfam_id = split_line[1]
            pfam_length = int(split_line[2])
            pfam_orf = split_line[3]
            pfam_start = int(split_line[17])
            pfam_end = int(split_line[18])
            pfam_hit_length = pfam_end - pfam_start
            aPfam = Pfam(pfam_name=pfam_name,
                         pfam_id=pfam_id,
                         length=pfam_length,
                         orf=pfam_orf,
                         start=pfam_start,
                         end=pfam_end,
                         hit_len=pfam_hit_length)
            self.append(aPfam)

    def ana_relations(self):
        """
        处理一个Hit里的所有Pfam，保留相互独立的pfam，且每个pfam都是最长的
        :return: 一个list，里面的元素为Pfam对象
        """
        # 对列表按长度进行降序排列
        sorted_list = sorted(self, key=lambda x: x.hit_len, reverse=True)
        try:
            cleaned_list = [sorted_list[0]]
        except IndexError:
            # hmmscan 扫描结果为空
            return None
        for i in range(1, len(sorted_list)):  # 从第二个开始比较
            current_pfam = sorted_list[i]
            independent = True
            for j in range(len(cleaned_list)):  # 当前的pfam与clean_list中的pfam对比
                existing_pfam = cleaned_list[j]
                if not current_pfam.related_to(existing_pfam):  # 与已有的Pfam对象不独立
                    if current_pfam.hit_len > existing_pfam.hit_len:
                        cleaned_list[j] = current_pfam  # 替换为较长的Pfam对象
                    else:
                        independent = False
                    break
            if independent:
                cleaned_list.append(current_pfam)
        return cleaned_list
