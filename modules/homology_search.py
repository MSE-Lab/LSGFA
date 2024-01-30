import os
import igraph
from igraph import Graph
from modules.utils import *
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

	def homology_search(self, query_dir, db_dir, res_dir, threads):
		# 返回搜索的命令
		cmd_o = CmdManger(process=self.method, thread=threads)
		if self.method == "mmseqs":
			db = os.path.join(query_dir, f'{self.name}.fa')
		else:
			db = os.path.join(db_dir, self.name)
		query = os.path.join(query_dir, f'{self.name}.fa')  # 要进行比较的fasta文件
		res = os.path.join(res_dir, f'{self.name}.txt')  # 比较结果
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
