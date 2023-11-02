import math
import os
import igraph
from igraph import Graph
from modules.utils import *
from collections import Counter
import leidenalg as la


class CommonD:
	# 用于存放同样的domain
	def __init__(self, name: str = "", genes: list = None, method: str = ""):
		self.name = name
		self.genes = genes
		self.method = method

	def __str__(self):
		return f'{self.name}.fa'

	def __len__(self):
		return len(self.genes)

	def __repr__(self):
		return self.name

	def write_seqs(self, seqInfo, out_path):
		fasta = '\n'.join([f'>{gene}\n{seqInfo[gene]}' for gene in self.genes])
		FileOperator(name=f'{self.name}.fa', dir_=out_path, data=fasta).write()

	def homology_search(self, query_dir, db_dir, res_dir):
		# 同源序列的搜索
		cmd_o = CmdManger(process=self.method)
		if self.method == "mmseqs":
			db = os.path.join(query_dir, f'{self.name}.fa')
		else:
			db = os.path.join(db_dir, self.name)
		query = os.path.join(query_dir, f'{self.name}.fa')  # 要进行比较的fasta文件
		res = os.path.join(res_dir, f'{self.name}-{self.name}.txt')  # 比较结果
		if '.fa' == db[:-3]:  # 因为mmseq不需要db，它的db是自己
			db_cmd = None
		else:
			cmd_o.make_db(input_name=query, db=db)
			db_cmd = cmd_o.cmd
		cmd_o.homology_searching(query=query, db=db, out_name=res)
		search_cmd = cmd_o.cmd
		return db_cmd, search_cmd

	def parse_homology_search(self, res_path):  # 处理同源搜索的结果
		paried_dict = {}  # 存放两向的结果
		file = FileOperator(name=f'{self.name}-{self.name}.txt', dir_=res_path)
		file.read()
		ids = set()  # 存放在identity上满足条件的id
		for line in file.data:  # 处理结果
			try:
				row = line.strip("\n").split("\t")
				id1 = row[0]
				id2 = row[1]
				if id1 != id2:
					ID = f'{id1}\t{id2}'  # 因为做query和做db的结果不同，所以要计算两项的值
					ID_reverse = f'{id2}\t{id1}'
					bitscore = float(row[11])
					ident = float(row[2])
					# e_value = -math.log10(float(row[10]))
					if ident >= 0:
						ids.add(id1)
						ids.add(id2)
						if ID not in paried_dict and ID_reverse not in paried_dict:
							paried_dict[ID] = bitscore
						elif ID not in paried_dict and ID_reverse in paried_dict:
							paried_dict[ID_reverse] = (bitscore + paried_dict[ID_reverse]) / 2
						elif ID in paried_dict and ID_reverse not in paried_dict:
							paried_dict[ID] = (bitscore + paried_dict[ID]) / 2  # 计算两项的平均值
			except (IndexError, ValueError):
				sys.stderr.write(
					"\nERROR: Query or hit sequence ID in BLAST results file was missing or incorrectly formatted.\n")
				raise
		paired_genes = ids	# 两两配对的gene
		edges = []
		weights = []
		for id,weight in paried_dict.items():
			edges.append(id.split("\t"))  # 获得一对边的列表
			weights.append(weight)  # 获取权重
		if paried_dict:	# 如果dict内容存在
			return paired_genes, weights, edges
		else:
			return paired_genes, None, None


class PfamG(list):
	# 存放图的list，里面每个元素是一个CommonD
	def __init__(self, graph: Graph = None, method: str = ""):
		super(PfamG, self).__init__()
		self.graph = graph
		self.method = method
		self._get_pfam_cluster()
		self._get_pfam_cluster_append()

	def _get_pfam_combination(self):
		# 得到所有pfam的组合
		return sorted(set(self.graph.vs['pfam']))

	def _get_pfam_cluster_append(self):
		# 把连通图间连接的点当作是一个cc
		for n, cc in enumerate(self.graph.components()):
			if len(cc) > 1:
				loci = self.graph.vs[cc]
				c_d = CommonD(name=f'CC{n:0>7}', genes=loci['name'], method=self.method)
				# 初始化
				self.append(c_d)

	def _get_pfam_cluster(self):
		# 连通图单独计算
		for n, pfam_c in enumerate(self._get_pfam_combination()):
			if pfam_c == 'None':
				print(f'PD{n:0>7}', len(self.graph.vs.select(pfam=pfam_c)))
			genes = [node['name'] for node in self.graph.vs.select(pfam=pfam_c)]  # 获得当前pfam的名字
			common_domain = CommonD(name=f'PD{n:0>7}', genes=genes, method=self.method)
			# 初始化
			self.append(common_domain)

	def homology_search_commands(self, query_path, db_path, res_path, seq_info):
		db_cmds = list()
		search_cmds = list()
		for cluster in self:
			cluster: CommonD
			cluster.write_seqs(seqInfo=seq_info, out_path=query_path)
			db_cmd, search_cmd = cluster.homology_search(query_dir=query_path, db_dir=db_path, res_dir=res_path)
			# 解包元组
			db_cmds.append(db_cmd)
			search_cmds.append(search_cmd)
		return db_cmds, search_cmds

	def la_graph(self, res_dir, graph_file_name, max_genome):
		genes_all = set()  # 存放所有的基因
		paired_genes_all = set()
		g = igraph.Graph()
		edges = []
		weights = []
		for cluster in self:
			cluster: CommonD
			genes_all = genes_all.union(set(cluster.genes))
			paired_genes, weight, edge_c = cluster.parse_homology_search(res_dir)
			paired_genes_all = paired_genes_all.union(paired_genes)
			if weight is not None:	 # 如果边的权重为0，则不添加边
				edges.extend(edge_c)
				weights.extend(weight)
		g.add_vertices(list(genes_all))	 # 添加节点
		g.add_edges(edges)  # 添加边
		g.es['weight'] = weights
		g.write_gml(graph_file_name)	 # 生成图
		partition = la.find_partition(g, partition_type=la.RBERVertexPartition, weights='weight',
									resolution_parameter=1)
		print(f'partition numbers: {len(partition)}')
		SOG = 0
		HOG = 0
		for i in partition:
			genome = [n['name'].split('|')[0] for n in g.vs[i]]
			genome_set = set(genome)
			if len(genome) >= max_genome:
				if len(genome_set) == len(genome):
						SOG += 1
				else:
					HOG += 1
		print(f'SOG: {SOG}')
		print(f'HOG: {HOG}')
		return partition


class OGs:
	# 用于存放OG的属性和方法
	def __init__(self, name: str = "", seqs: list = None):
		self.name = name
		self.genes = seqs
		self.type = None

	def __len__(self):
		return len(self.genes)

	def __repr__(self):
		return self.name

	def _genomes(self):
		return set([g.split("|")[0] for g in self.genes])

	def len_genome(self):
		return len(self._genomes())

	def to_fasta(self, seq_path, seq_info):
		fasta = '\n'.join([f'>{gene}\n{seq_info[gene]}' for gene in self.genes])
		fasta_o = FileOperator(name=f'{self.name}.fa', dir_=seq_path, data=fasta)
		fasta_o.write()

	def mark_og(self, max_genome):
		# 给og打标记
		if self.len_genome() >= max_genome:
			if len(self.genes) == self.len_genome():
				self.type = "SOG"
			else:
				self.type = "HOG"
		else:
			self.type = "Others"  # og内的基因组数量比输入的基因组数量少


class Clusters(list):
	def __init__(self, cluster_file_name: str = "", max_genome: int = None):
		super(Clusters, self).__init__()
		self._read_mcl_cluster(cluster_file_name)
		self.max = max_genome

	def _read_mcl_cluster(self, file):
		cluster_file = FileOperator(name=file)
		cluster_file.read()
		for n, line in enumerate(cluster_file.data):  # n取位置，line取内容
			og_name = f'OG{n + 1:0>7}'
			og = OGs(name=og_name, seqs=line.strip("\n").split("\t"))
			self.append(og)

	def write_og(self, seq_path, seq_info):
		# 输出og的文件
		new_og = ""
		for og in self:
			og: OGs
			og.mark_og(self.max)  # 打标记
			genes = "\t".join(og.genes)
			new_og += f'{og.name}\t{og.type}\t{genes}\n'
			if og.type == "SOG":
				out_path = os.path.join(seq_path, 'SOG')
			else:
				out_path = os.path.join(seq_path, 'Others')  # type为Other和HOG的
			os.makedirs(out_path, exist_ok=True)
			og.to_fasta(out_path, seq_info)
		o = FileOperator(name='Orthogroups.txt', dir_=seq_path, data=new_og)
		o.write()
