from itertools import combinations
from igraph import Graph


def get_edges(members):

	pf_dic = {}
	orf_list = []

	for member in members:
		with open(member) as f:
			contents = f.readlines()
		content_list = [content.split('\t') for content in contents]

		for i in range(len(content_list)):
			orf = content_list[i][0]
			orf_list.append(orf)
			pf = list(set(content_list[i][1].rstrip().split(';')))
			for j in pf:
				if j != 'None':
					if j not in pf_dic:
						pf_dic[j] = [orf]
					else:
						pf_dic[j].append(orf)

	edges_set = set()
	for pfam in pf_dic.keys():
		if len(pf_dic[pfam]) > 1:
			edges_set = edges_set | set(combinations(sorted(pf_dic[pfam]), 2))
		# 对pfam对应的orf列表组合，然后取并集，得到所有边的情况
		# sorted后可以消除(a,b)和(b,a)不一致的情况
	edges = list(edges_set)  # 边
	vertices = list(set(orf_list))  # 点
	return vertices, edges


def build_graph(orf_id, edges):
	"""
	建立网络
	:return: g，网络文件
	"""
	g = Graph()
	g.add_vertices(orf_id)
	g.add_edges(list(edges))
	g.write_graphml('OGGraph.graphml')
	return g
