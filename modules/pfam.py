class Pfam:

	def __init__(self, pfam_id, percent, pfam_name=None, start=None, end=None, hit_len=None):
		self.name = pfam_name
		self.id = pfam_id
		self.start = start
		self.end = end
		self.hit_len = hit_len
		self.percent = percent

	def __str__(self):
		return self.name


class Hits(list):
	# 存放每条序列hit到的内容，每个元素是一个pfam的对象

	def __init__(self, hmmscan_out: str = ""):
		super().__init__()
		self.raw = []
		self._handle_raw(hmmscan_out)
		self.ana_relations()
		self.sum_percent()

	def _handle_raw(self, hmmscan_out):
		with open(hmmscan_out) as f:
			content = f.readlines()  # hits是一个list，存放读取未处理的内容
			hits = [i for i in content if not i.startswith("#")]

		for i in hits:
			split_line = [_ for _ in i.split(' ') if _ != '']
			pfam_id = split_line[1]
			pfam_start = int(split_line[17])
			pfam_end = int(split_line[18])
			pfam_hit_length = pfam_end - pfam_start
			seq_len = int(split_line[5])
			pfam_percent = pfam_hit_length / seq_len
			aPfam = Pfam(pfam_id=pfam_id,
						 start=pfam_start,
						 end=pfam_end,
						 hit_len=pfam_hit_length,
						 percent=pfam_percent)
			self.raw.append(aPfam)

	@staticmethod
	def related_to(pfam_a, pfam_b):
		# 判断两个Pfam之间的重叠关系
		if pfam_a.hit_len >= pfam_b.hit_len:
			long = pfam_a
			short = pfam_b
		else:
			long = pfam_b
			short = pfam_a
		if long.start >= short.end or long.end <= short.start:  # 独立的两个
			return True
		else:  # 当两个为包含、相等或overlap时
			return False

	def ana_relations(self):
		"""
        处理一个Hit里的所有Pfam，保留相互独立的pfam，且每个pfam都是最长的
        :return: 一个list，里面的元素为Pfam对象
        """
		# 对列表按长度进行降序排列
		sorted_list = sorted(self.raw, key=lambda x: x.hit_len, reverse=True)
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
				if not self.related_to(current_pfam, existing_pfam):  # 与已有的Pfam对象不独立
					if current_pfam.hit_len > existing_pfam.hit_len:
						cleaned_list[j] = current_pfam  # 替换为较长的Pfam对象
					else:
						independent = False
					break
			if independent:
				cleaned_list.append(current_pfam)
		self.raw = cleaned_list

	def sum_percent(self):
		# 用于处理相同pfam的percent
		domains = {}
		for i in self.raw:
			if i.id in domains:
				p_ = i.percent
				domains[i.id] = p_ + i.percent
			else:
				domains[i.id] = i.percent
		for k, v in domains.items():
			aPfam = Pfam(pfam_id=k,
						 percent=v)
			self.append(aPfam)
