import warnings
from itertools import combinations
import igraph
from modules import panproteome
import leidenalg as la
from collections import Counter
from modules.utils import *

warnings.filterwarnings("ignore")


class DomainType:
    # 用于存放Domain的类型
    def __init__(self, proteins: list, name: str = ""):
        self.name = name  # pfam_type
        self.proteins = proteins  # 该pfam_type所包含的protein对象
        self.domain = self._get_domains()

    def __str__(self):
        return self.name

    # def __repr__(self):
    #     return self.name

    def _get_domains(self):
        return set(self.name.split(","))

    def sharing_domain_loci(self, sharing_domain, domain_length_cov):
        # 用于判断是否满足长度阈值
        matching_protein = 0
        for protein in self.proteins:  # overlap的pfam
            domain_cov_sum = sum([pfam.percent for pfam in protein.domain if pfam.id in sharing_domain])
            if domain_cov_sum >= domain_length_cov:
                matching_protein += 1
        return matching_protein


class PGraph:
    def __init__(self, proteomes: panproteome):
        super(PGraph, self).__init__()
        self.domain_type = []  # 保存pfam_type属性
        self.connection = {}  # 用字典来保存connection属性，key是边的两个节点，value是权重
        self.graph = None  # 存储该PGraph的图
        domain_dict = {}
        protein_num = 0
        for proteome in proteomes:
            for protein in proteome:  # 对protein重新分类，实例化DomainType
                protein_num += 1
                if sum([i.percent for i in protein.domain]) <= 0.6:
                    pfam_ = "None"
                else:
                    pfam_ = ','.join(sorted([i.id for i in protein.domain]))
                if pfam_ in domain_dict:
                    domain_dict[pfam_].append(protein)
                else:
                    domain_dict[pfam_] = [protein]
        for k, v in domain_dict.items():
            aDomainType = DomainType(name=k, proteins=v)
            self.domain_type.append(aDomainType)
        message(text=f"Genes number: {protein_num}", label='Information')
        message(text=f"None pfam genes number: {len(domain_dict['None'])}", label='Information')
        message(text=f"DomainType number: {len(self.domain_type)}", label='Information')
        self._get_edges()  # 获取边的信息

    def get_domain_type(self):
        return list(i.name for i in self.domain_type)

    def _compared_domain_component_pairwise(self):
        domain_components = [k for k in self.domain_type if k.name != 'None']
        return list(combinations(domain_components, 2))

    def _get_edges(self, sharing_lencov=0.5):  # overlap的pfam占各自序列的0.5以上
        edges = dict()
        for pfam_type1, pfam_type2 in self._compared_domain_component_pairwise():  # 两个不同domain的CC的组合
            sharing_domains = self.sharing_domain(pfam_type1.domain, pfam_type2.domain)  # CC间有共同的PF
            if sharing_domains is not None:  # 如果有共享的pfam
                matching_protein_1 = pfam_type1.sharing_domain_loci(sharing_domains, sharing_lencov)
                matching_protein_2 = pfam_type2.sharing_domain_loci(sharing_domains, sharing_lencov)
                # 计算权重
                weight = (matching_protein_1*matching_protein_2)/(len(pfam_type1.proteins)*len(pfam_type2.proteins))
                edges_comb = tuple([pfam_type1.name, pfam_type2.name])
                edges[edges_comb] = weight
        self.connection = edges

    @staticmethod
    def sharing_domain(pf_type1, pf_type2):
        # 判断两个pf_type间是否有重合
        sharing_domain = list(set(pf_type1) & set(pf_type2))
        len_sharing_domain = len(sharing_domain)
        return sharing_domain if len_sharing_domain > 0 else None

    def generate_graph(self):  # 构建网路
        domain_type_graph = igraph.Graph()
        vs = self.domain_type
        es = list(self.connection.keys())
        weigth = list(self.connection.values())
        domain_type_graph.add_vertices(vs)  # 添加点
        domain_type_graph.vs['domain_type'] = [i.name for i in self.domain_type]  # 给点添加属性
        es_index = [(domain_type_graph.vs.find(domain_type=edge[0]).index,
                     domain_type_graph.vs.find(domain_type=edge[1]).index)
                    for edge in es]  # 构建边的列表
        domain_type_graph.add_edges(es_index)
        domain_type_graph.es['weight'] = weigth   # 给节点添加pfam的属性
        return domain_type_graph

    def la_find_partition(self):  # 社区发现
        self.graph = self.generate_graph()
        partition = la.find_partition(self.graph, partition_type=la.CPMVertexPartition,
                                      weights='weight',
                                      resolution_parameter=0.9)
        message(text=f"Partition number: {len(partition)}", label='Information')
        partition_genes = []
        for community in partition:  # 获取每个社区内的蛋白
            community_subgraph = self.graph.subgraph(community)
            protein_lists = [node['name'].proteins for node in community_subgraph.vs]
            genes_in_community = [protein for proteins in protein_lists for protein in proteins]
            partition_genes.append(genes_in_community)
        return partition_genes

    def partition_p_og(self, partition, max_genome, out_dir):
        result = ''
        scc_num = 0
        sum_p_og = 0
        for community in partition:
            community_subgraph = self.graph.subgraph(community)
            mode_in_community = community_subgraph.vs['domain_type']
            protein_lists = [node['name'].proteins for node in community_subgraph.vs]
            genes_in_community = [protein.name for proteins in protein_lists for protein in proteins]
            genomes_list = [n.split('|')[0] for n in genes_in_community]
            genomes_set = set(genomes_list)  # 判断该community是否是核心
            if len(genomes_set) >= max_genome:
                p_OG = min(Counter(genomes_list).values())  # 获取每个PG的大小
                sum_p_og += p_OG
            else:
                p_OG = 0
            result += f'SCC{scc_num:0>7}\t{len(mode_in_community)}\t{len(genes_in_community)}\t{p_OG}\n'
            scc_num += 1
        message(text=f'Number of Communities:{scc_num+1}', label='Information')
        message(text=f'Number of Potential_og:{sum_p_og}', label='Information')

        with open(os.path.join(out_dir, 'partition_Potential_og.txt'), 'w') as f:
            title = '#S_ID\tMode_num\tS_size\tPotential_og\n'
            f.write(title)
            f.write(result)
