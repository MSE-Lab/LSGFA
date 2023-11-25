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
        self.add_loci = []

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
        for loci, domain_cov_info in self.domain_data.items():
            domain_cov_sum = np.sum(
                np.array(domain_cov_info['LenCov'])[np.isin(np.array(domain_cov_info['Domain']), sharing_domain)])
            if domain_cov_sum >= domain_length_cov:
                self.add_loci.append(loci)


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
                self.setdefault(pfam_component, dict()).update(seq_info)
        self._generate_domain_info()

    def _compared_domain_component_pairwise(self):
        domain_components = [k for k in self.keys() if len(k.split(";")) > 1]
        return list(combinations(domain_components, 2))

    def _generate_domain_info(self):
        domain_component = {}
        for pfam_component, seq_info in self.items():
            pfam_ = DomainType(name=pfam_component, domain_data=seq_info)   # 实例化DomainType
            domain_component[pfam_component] = pfam_
        self.domains = domain_component

    def get_full_connected_edges(self):
        for pfam_, pfam_info in self.domains.items():
            self.genes.extend(pfam_info.get_sequences_ids())
            self.node_attribute.extend([pfam_] * len(pfam_info))

    def get_append_edges(self, sharing_cov=0.5, len_cov=0.8):
        for pf_1, pf_2 in self._compared_domain_component_pairwise():
            # 两个不同domain的CC的组合
            sharing_domains = self.sharing_domain(pf_1, pf_2, sharing_cov)  # CC间有共同的PF
            if sharing_domains is not None:  # 如果有共享的pfam
                pf1_o = self.domains[pf_1]  # 索引两个domain对应的DomainType对象
                pf2_o = self.domains[pf_2]
                pf1_o.sharing_domain_loci(sharing_domains, len_cov)
                pf2_o.sharing_domain_loci(sharing_domains, len_cov)  # 判断有重合的两个pf_type间是否满足长度阈值
                self.related_edges.extend(product(pf1_o.add_loci, pf2_o.add_loci))

    @staticmethod
    def sharing_domain(pf_type1, pf_type2, domain_sharing_cov):
        # 判断有重合的两个pf_type间是否满足数量阈值
        domain_1 = pf_type1.split(';')
        domain_2 = pf_type2.split(';')
        sharing_domain = list(set(domain_1) & set(domain_2))
        sharing_cov = min([len(sharing_domain) / len(domain_1), len(sharing_domain) / len(domain_2)])
        return sharing_domain if sharing_cov >= domain_sharing_cov else None

    @staticmethod
    def generate_final_graph(vs, es, **kwargs):
        g = igraph.Graph()
        g.add_vertices(vs)
        g.add_edges(es)
        g.vs['pfam'] = kwargs['pfam']   # 给节点添加pfam的属性
        return g
