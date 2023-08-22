from itertools import combinations
from igraph import Graph
import leidenalg as la


def co_index(lista, listb):
	# 用于判断两个ORF之间Pfam的关联度
	set_len = len(set(lista) & set(listb))
	if set_len == 0:
		co_number = 0
	else:
		co_number = (set_len/len(lista) + set_len/len(listb))*0.5
	return co_number


def get_edges(members):
	# 用于获取点和边
	# members是.Pfam的文件

	orf_pfam_dic = {}
	orf_list = []
	edges_list = []

	for member in members:
		with open(member) as f:
			contents = f.readlines()
		content_list = [content.split('\t') for content in contents]

		# 建立orf对应的pfam的字典
		for i in range(len(content_list)):
			orf = content_list[i][0]
			orf_list.append(orf)
			pf = list(set(content_list[i][1].rstrip().split(';')))
			if pf != ['None']:
				orf_pfam_dic[orf] = pf

	# 判断两个orf之间是否可以连线
	key_list = sorted(list(orf_pfam_dic.keys()))
	n = len(key_list)
	for i in range(n):
		lista = orf_pfam_dic[key_list[i]]
		for j in range(i+1, n):
			listb = orf_pfam_dic[key_list[j]]
			if co_index(lista, listb) >= 0.6:
				edges_list.append((key_list[i], key_list[j]))

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
