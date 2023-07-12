from itertools import combinations
from igraph import Graph

def get_edges(self):  # 这个self是一个list，每个元素的一个pfam的对象

	# 建立各个pfam_id对应的orf_id, pfam_id为key
	pf_dic = {}
	for i in range(len(self)):
		if self[i].id not in pf_dic:
			pf_dic[self[i].id] = [self[i].orf]
		else:
			pf_dic[self[i].id].append(self[i].orf)

	edges_set = set()
	for pfam in pf_dic.keys():
		if len(pf_dic[pfam]) > 1:
			edges_set = edges_set | set(combinations(sorted(self[pfam]), 2))
		# 对pfam对应的orf列表组合，然后取并集，得到所有边的情况
		# sorted后可以消除(a,b)和(b,a)不一致的情况
	edges = list(edges_set)  # 边
	vertices = pf_dic.keys()  # 点
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
