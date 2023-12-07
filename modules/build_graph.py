import glob
import os
import warnings
from itertools import combinations, product
import igraph
import numpy as np
from modules.utils import *

warnings.filterwarnings("ignore")


class DomainType:
    # 用于存放Domain的类型
    def __init__(self, name: str = "", domain_data: dict = None):
        self.name = name
        self.domain_data = domain_data  # 存放domain的长度信息
        self.domain = self._get_domains()
        # self.add_loci = []

    def __str__(self):
        return self.name

    def __repr__(self):
        print(self.domain)

    def __len__(self):
        return len(self.get_sequences_ids())

    def _get_domains(self):
        return set(self.name.split(";"))

    def get_sequences_ids(self):
        return list(self.domain_data.keys())

    def sharing_domain_loci(self, sharing_domain, domain_length_cov):
        # 用于判断是否满足长度阈值
        add_loci = []
        for loci, domain_cov_info in self.domain_data.items():
            domain_cov_sum = np.sum(
                np.array(domain_cov_info['LenCov'])[np.isin(np.array(domain_cov_info['Domain']), sharing_domain)])
            if domain_cov_sum >= domain_length_cov:
                # self.add_loci.append(loci)
                add_loci.append(loci)
        return add_loci


class PGraph(dict):
    # 存放Pfam的连通图
    def __init__(self, pfam_res_dir: str = "", domains: dict = None):
        super(PGraph, self).__init__()
        members = glob.glob(os.path.join(pfam_res_dir, '*.pfam'))
        self.genes = []
        self.node_attribute = []
        self.related_edges = []
        self.domains = domains
        for m in members:   # 读取pfam的结果
            json_data = FileOperator(os.path.basename(m), pfam_res_dir, "json")
            json_data.read()
            for pfam_component, seq_info in json_data.data.items():
                for gene, value in seq_info.items():
                    if sum(value['LenCov']) <= 0.6:
                        pfam_ = "None"
                    else:
                        pfam_ = pfam_component
                    self.setdefault(pfam_, dict()).update({gene: value})
        self._generate_domain_info()

    def _compared_domain_component_pairwise(self):
        domain_components = [k for k in self.keys() if k != 'None']
        return list(combinations(domain_components, 2))

    def _generate_domain_info(self):
        domain_component = {}
        for pfam_component, seq_info in self.items():
            pfam_ = DomainType(name=pfam_component, domain_data=seq_info)   # 实例化DomainType
            domain_component[pfam_component] = pfam_
        self.domains = domain_component

    def get_full_connected_edges(self):
        # 获得全连通图，给同样domain的节点添加pfam的信息
        for pfam_, pfam_info in self.domains.items():
            if pfam_ != 'None':
                self.genes.extend(pfam_info.get_sequences_ids())
                # self.related_edges.extend(list(combinations(pfam_info.get_sequences_ids(), 2)))
                self.node_attribute.extend([pfam_] * len(pfam_info))
            else:
                self.genes.extend(pfam_info.get_sequences_ids())
                self.node_attribute.extend([pfam_] * len(pfam_info))

    def get_append_edges(self, sharing_lencov=0.5):
        for pf_1, pf_2 in self._compared_domain_component_pairwise():
            # 两个不同domain的CC的组合
            sharing_domains = self.sharing_domain(pf_1, pf_2)  # CC间有共同的PF
            if sharing_domains is not None:  # 如果有共享的pfam
                pf1_o = self.domains[pf_1]  # 索引两个domain对应的DomainType对象
                pf2_o = self.domains[pf_2]
                add_loci_1 = pf1_o.sharing_domain_loci(sharing_domains, sharing_lencov)
                add_loci_2 = pf2_o.sharing_domain_loci(sharing_domains, sharing_lencov)
                # 判断有重合的两个pf_type间是否满足长度阈值
                self.related_edges.extend(product(add_loci_1, add_loci_2))

    @staticmethod
    def sharing_domain(pf_type1, pf_type2):
        # 判断两个pf_type间是否有重合
        domain_1 = pf_type1.split(';')
        domain_2 = pf_type2.split(';')
        sharing_domain = list(set(domain_1) & set(domain_2))
        len_sharing_domain = len(sharing_domain)
        return sharing_domain if len_sharing_domain > 0 else None

    @staticmethod
    def generate_final_graph(vs, es, **kwargs):
        g = igraph.Graph()
        g.add_vertices(vs)
        g.add_edges(es)
        g.vs['pfam'] = kwargs['pfam']   # 给节点添加pfam的属性
        return g
