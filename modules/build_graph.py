from collections import defaultdict, deque
from itertools import combinations, product
from igraph import Graph
import leidenalg as la
import json


def co_index(lista, listb):
	# 用于判断两个ORF之间Pfam的关联度
	set_len = len(set(lista) & set(listb))
	if set_len == 0:
		co_number = 0
	else:
		co_number = min(set_len/len(lista), set_len/len(listb))
	return co_number


def co_len(orf_a, orf_b):
	# 用于判断两个ORF之间共享的pfam之间hit到的最短长度
	len_a = 0
	len_b = 0
	co_pfam_list = list(set(orf_a.keys() & orf_b.keys()))
	for pf in co_pfam_list:
		len_a = len_a + orf_a[pf]
		len_b = len_b + orf_b[pf]
	min_len = min(len_a, len_b)
	return min_len


def get_orf_dic(members):
	# 从.pfam文件中获取内容，处理为字典
	# members是.Pfam的文件

	orf_dic = {}
	for member in members:
		with open(member, 'r') as f:
			data = json.load(f)
			orf_dic.update(data)
	return orf_dic


def get_vertices(orf_dic):
	# 用于获取点，即orf的名字
	vertices = list(key for key in orf_dic.keys())
	return vertices


def get_edges(orf_dic, d_co_index = 0.5, d_co_len = 0.5):
	# 用于获取边
	edges_list = set()

	# 获取不同类型的pfam对应的orf
	pf_dic = {}
	for orf,p in orf_dic.items():
		# o为  'orf1'
		# p为  {'PF1': 0.1, 'PF2': 0.18, 'PF3': 0.5}}
		try:
			pf_type = tuple(sorted(list(p.keys())))
			pf_dic.setdefault(pf_type, []).append(orf)
		except AttributeError:  # 当注释结果为None时
			pass

	# 对不同类型的pfam组合的处理
	pf_list = list(sorted(pf_dic.keys()))

	# 先处理pfam相同的类
	for p_t in pf_list:  # p_t是pfam组合的元组
		orf_list = pf_dic[p_t]  # ['orf1', 'orf2']
		orf_len = len(orf_list)
		if orf_len >1:  # 该类只有一个orf则不予考虑
			for i in range(orf_len):
				orf_i = orf_dic[orf_list[i]]  # orf_i的结构为{'pf1':0.1,'pf2':0.2}
				for j in range(i+1, orf_len):
					orf_j = orf_dic[orf_list[j]]
					if co_len(orf_i, orf_j) > d_co_len:  # 如果长度满足设定值，则添加边
						edges_list.add(tuple(sorted([orf_list[i], orf_list[j]])))
	# 然后对不同pfam类型间的orf进行处理
	for i in range(len(pf_list)):
		pf_i = pf_list[i]
		for j in range(i+1, len(pf_list)):
			pf_j = pf_list[j]

			if co_index(pf_i, pf_j) >= d_co_index:  # 当co_index不满足设定值时，就不予考虑
				orf_list_i = pf_dic[pf_i]  # ['orf1', 'orf2']
				orf_list_j = pf_dic[pf_j]
				for n in range(len(orf_list_i)):
					orf_i = orf_dic[orf_list_i[n]]  # orf_i的结构为{'pf1':0.1,'pf2':0.2}
					for m in range(len(orf_list_j)):
						orf_j = orf_dic[orf_list_j[m]]
						if co_len(orf_i, orf_j) > d_co_len:  # 如果长度满足设定值，则添加边
							edges_list.add(tuple(sorted([orf_list_i[n], orf_list_j[m]])))
	return edges_list


def build_graph(orf_id, edges):
	# 用于建立网络
	g = Graph()
	g.add_vertices(orf_id)
	g.add_edges(list(edges))
	g.write_graphml('OGGraph.graphml')
	return g


def graph_split(self):
	# 对网络进行社区发现
	subgraphs = self.components().subgraphs()
	communities = []
	for subgraph in subgraphs:
		partition = la.find_partition(subgraph, la.ModularityVertexPartition)
		for community in partition:
			community_names = [subgraph.vs[node]["name"] for node in community]
			communities.append(community_names)
	return communities
