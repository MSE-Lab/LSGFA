import os
import glob
from itertools import combinations, product
import numpy as np
import warnings
from modules.utils import *

warnings.filterwarnings("ignore")


class DomainType:
    def __init__(self, name: str = "", domain_data: dict = None):
        self.name = name
        self.domain_data = domain_data
        self.domain = self._get_domains()
        self.add_loci = []

    def __str__(self):
        return self.name

    def __repr__(self):
        print(self.domain)

    def __len__(self):
        return len(self.domain)

    def _get_domains(self):
        return set(self.name.split(";"))

    def _get_sequences_ids(self):
        return list(self.domain_data.keys())

    def sharing_domain_loci(self, sharing_domain, domain_length_cov):
        # 两个orf间重叠的pfam
        for loci, domain_cov_info in self.domain_data.items():
            # 用于判断用于计算重叠的pfam的长度之和
            domain_cov_sum = np.sum(
                np.array(domain_cov_info['LenCov'])[np.isin(np.array(domain_cov_info['Domain']),
                                                            sharing_domain)])
            if domain_cov_sum >= domain_length_cov:
                # ？？？？？
                self.add_loci.append(loci)

    def generate_full_graph(self):
        # domain相同的序列，内部自己形成全连通的图
        return list(combinations(self._get_sequences_ids(), 2))


class PGraph(dict):
    def __init__(self, pfam_res_dir: str = "", domains: dict = None):
        super(PGraph, self).__init__()
        members = glob.glob(os.path.join(pfam_res_dir, '*.pfam'))
        self.genes = []
        self.related_edges = []
        self.domains = domains
        for m in members:  # m表示一个注释pfam的文件
            json_data = FileOperator(os.path.basename(m), pfam_res_dir, "json")
            json_data.read()
            for pfam_component, seq_info in json_data.data.items():
                self.setdefault(pfam_component, dict()).update(seq_info)
                self.genes.extend(seq_info.keys())
        self._generate_domain_info()

    def _compared_domain_component_pairwise(self):
        domain_components = [k for k in self.keys() if len(k.split(";")) > 1]
        # ？？？？
        return list(combinations(domain_components, 2))

    def _generate_domain_info(self):
        # 获取domain的信息
        domain_component = {}
        for pfam_component, seq_info in self.items():
            pfam_ = DomainType(name=pfam_component, domain_data=seq_info)
            domain_component[pfam_component] = pfam_
        self.domains = domain_component

    def get_full_connected_edges(self):
        for pfam_, pfam_info in self.domains.items():
            if pfam_ != 'None':
                self.related_edges.extend(pfam_info.generate_full_graph())

    def get_append_edges(self, sharing_cov, len_cov):
        for pf_1, pf_2 in self._compared_domain_component_pairwise():
            sharing_domains = self.sharing_domain(pf_1, pf_2, sharing_cov)
            if sharing_domains is not None:
                pf1_o = self.domains[pf_1]
                pf2_o = self.domains[pf_2]
                pf1_o.sharing_domain_loci(sharing_domains, len_cov)
                pf2_o.sharing_domain_loci(sharing_domains, len_cov)
                self.related_edges.extend(product(pf1_o.add_loci, pf2_o.add_loci))

    @staticmethod
    def sharing_domain(pf_type1, pf_type2, domain_sharing_cov):
        # 比较两个domain间共享pfam的长度
        domain_1 = pf_type1.split(';')
        domain_2 = pf_type2.split(';')
        sharing_domain = list(set(domain_1) & set(domain_2))
        sharing_cov = min([len(sharing_domain) / len(domain_1), len(sharing_domain) / len(domain_2)])
        return sharing_domain if sharing_cov >= domain_sharing_cov else None
