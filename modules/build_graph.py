from collections import defaultdict, deque
from itertools import combinations, product
from igraph import Graph
import leidenalg as la


def co_index(lista, listb):
	# 用于判断两个ORF之间Pfam的关联度
	set_len = len(set(lista) & set(listb))
	if set_len == 0:
		co_number = 0
	else:
		co_number = min(set_len/len(lista), set_len/len(listb))
	return co_number

def get_edges(members, intersect_t):
	# 用于获取点和边
	# members是.Pfam的文件

	orf_list = []
	pf_dic = defaultdict(list)
	edges_list = []

	for member in members:
		with open(member) as f:
			for line in f:
				orf, pfs = line.rstrip().split('\t')
				orf_list.append(orf)
				pfs = tuple(sorted(set(pfs.split(';'))))
				if pfs != ('None',):
					pf_dic[pfs].append(orf)

	# 判断两个orf之间是否可以连线
	pf_list = deque(sorted(pf_dic.keys()))

	while pf_list:
		# 相同pfam的彼此相连
		pf_i = pf_list.popleft()
		orf_i = pf_dic[pf_i]
		if len(orf_i) > 1:
			combination = list(combinations(sorted(orf_i), 2))
			edges_list.extend(combination)
		# 不同pfam的考虑阈值
		for pf_j in pf_list:
			orf_j = pf_dic[pf_j]
			if co_index(pf_i, pf_j) > intersect_t:
				combination = list(product(orf_i, orf_j))
				edges_list.extend(combination)

	# 边为：edges_list
	vertices = list(set(orf_list))  # 点
	return vertices, edges_list


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
