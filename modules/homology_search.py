import math
import os
import igraph
from igraph import Graph
from modules.utils import *
from collections import Counter
from pyfasta import Fasta
from modules.panproteome import Protein


class Ngroup:
	# 用于存放None这部分内容
	# Ngroup是一个dict，里面存放的是gene_id:seq
	def __init__(self, fasta_name, method: str = "", graph: Graph = None):
		super().__init__()
		self.content = []
		self.name = os.path.basename(fasta_name).split('.')[0]
		self.method = method
		self.graph = graph
		for gene_id, seq in Fasta(fasta_name).items():
			a_protein = Protein(name=gene_id, sequence=seq)
			self.content.append(a_protein)

	def __len__(self):
		return len(self.content)

	# def write_seqs(self, seqInfo, out_path):  # 写序列文件
	# 	fasta = '\n'.join([f'>{gene}\n{seqInfo[gene]}' for gene in self.genes])
	# 	FileOperator(name=f'{self.name}.fa', dir_=out_path, data=fasta).write()

	def homology_search(self, query_dir, db_dir, res_dir, threads):
		# 返回搜索的命令
		cmd_o = CmdManger(process=self.method, thread=threads)
		if self.method == "mmseqs":
			db = os.path.join(query_dir, f'{self.name}.fa')
		else:
			db = os.path.join(db_dir, self.name)
		query = os.path.join(query_dir, f'{self.name}.fa')  # 要进行比较的fasta文件
		res = os.path.join(res_dir, 'none_pfam_blast.txt')  # 比较结果
		if db[:-3] == '.fa':  # 因为mmseq不需要db，它的db是自己
			db_cmd = None
		else:
			cmd_o.make_db(input_name=query, db=db)
			db_cmd = [cmd_o.cmd]
		cmd_o.homology_searching(query=query, db=db, out_name=res)
		search_cmd = [cmd_o.cmd]
		return db_cmd, search_cmd

	def build_homology_graph(self, res_file):  # 用于处理blast的文件，构建网络
		group_edges = []  # 存放有hit的结果
		file = FileOperator(name=res_file)
		file.read()
		for line in file.data:  # 处理结果
			row = line.strip("\n").split("\t")
			id1 = row[0]
			id2 = row[1]  # 获取每行blast的两个id
			if id1 != id2:  # 当和其它序列有hit时
				group_edges.append(tuple(sorted([id1,id2])))

		vs = [i.name for i in self.content]  # 添加点（string
		ng_group = igraph.Graph(directed=False)
		ng_group.add_vertices(vs)  # 添加点，点是protein的name
		ng_group.vs['object'] = self.content  # 给点添加属性，object属性是protein的对象
		ng_group.add_edges(group_edges)
		self.graph = ng_group

	def get_partition_genes(self):
		# 获取全连通图
		partition_genes = []
		for cc in self.graph.components():  # 获取每个社区内的蛋白
			community_subgraph = self.graph.subgraph(cc)
			genes_in_cc = [node['object'] for node in community_subgraph.vs]  # 获取cc内的蛋白对象
			partition_genes.append(genes_in_cc)
		return list(partition_genes)  # 返回的是一个list of list，里面每个元素是一个社区内的所有节点


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
			genes = [node['name'] for node in self.graph.vs.select(pfam=pfam_c)]  # 获得当前pfam的名字
			common_domain = CommonD(name=f'PD{n:0>7}', genes=genes, method=self.method)
			# 初始化
			self.append(common_domain)

	def homology_search_commands(self, query_path, db_path, res_path, seq_info):
		db_cmds = list()
		search_cmds = list()
		for cluster in self:
			cluster: Ngroup
			cluster.write_seqs(seqInfo=seq_info, out_path=query_path)
			db_cmd, search_cmd = cluster.homology_search(query_dir=query_path, db_dir=db_path, res_dir=res_path)
			# 解包元组
			db_cmds.append(db_cmd)
			search_cmds.append(search_cmd)
		return db_cmds, search_cmds

	def mcl_abc(self, res_dir, abc_file_name, threads, inflation, out_name):
		genes_all = set()
		paired_genes_all = set()
		abc_all = ""
		edges = []

		for cluster in self:
			cluster: CommonD
			genes_all = genes_all.union(set(cluster.genes))
			paired_genes, abc, edge_c = cluster.parse_homology_search(res_dir)
			paired_genes_all = paired_genes_all.union(paired_genes)
			if abc is not None:
				abc_all += abc + '\n'
				edges.extend(edge_c)
		abc_o = FileOperator(name=os.path.basename(abc_file_name), dir_=os.path.dirname(abc_file_name), data=abc_all)
		abc_o.write()
		mcl_o = CmdManger(thread=str(threads))
		mcl_o.mcl(abc_file_name, inflation, out_name)
		mcl_cmd = mcl_o.cmd
		call_mcl_cmd = CallCmd([mcl_cmd], process_info="MCL cluster", parallel=False)
		call_mcl_cmd.processing()  # mcl聚类

		g = igraph.Graph()
		g.add_vertices(list(genes_all))
		g.add_edges(edges)
		n = 0
		# 以下部分是输出ssn网络CC的节点信息
		node_string = ''
		for cc in g.components():
			cc_num = f'CC{n:0>7}:'
			node_names = ' '.join([g.vs[node]['name'] for node in cc])
			n += 1
			node_string += f'{cc_num} {node_names}\n'
		cc_file = FileOperator(name='cc_node.txt', dir_=res_dir, data=node_string)
		cc_file.write()


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
