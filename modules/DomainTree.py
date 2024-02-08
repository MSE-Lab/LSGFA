#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2023-12-13 15:43
# @Author  : zhaoyu
# @File    : DomainTree.py
import os
from modules.utils import *
import subprocess
from multiprocessing import Pool
import re
from math import fsum
import glob
from ete3 import Tree

#GLOBAL VARIABLES. MODIFY IF NEEDED
sep = '|'
gsep = re.escape(sep)  # 对|进行转义


class Split:  # 用于存储分裂方式
    def __init__(self):
        self.vecs = None
        self.branch_length = None  # 支长
        self.support = None  # 支持度
        self.name = None
        self.label = []


class DomainTree:
    def __init__(self, name, proteins_list):  # 一个DomainTree存放一个community的tree
        self.name = name  # 储存该树的编号
        self.proteins_list = proteins_list  # 存放蛋白对象
        self.gene_num = len(proteins_list)  # 记录有几个基因
        self.index = None  # 用于存放cd-hit结果
        self.leaves = None
        self.newick = None
        self.splits = []  # 里面存放的是split的class
        self.ortho = []  # 用于记录子树
        self.costs = {}  # 用于储存每个基因的内部分数

    def get_newick(self, newick):
        self.newick = newick.strip('\n')  # 如果有多个树，就分来做
        self.leaves = self.get_leaves(newick)  # 叶节点的list
        for leaf in self.leaves:  # 初始化，假定不存在内部同源，即该物种只有一个基因
            self.costs[leaf] = 1.0

    @staticmethod
    def get_leaves(newick):
        # 正则表达获取newick文件中的叶节点名称，就是所有参与画树的基因
        pattern = "(?<=[,\(])\w.*?(?=[,:;\)])"
        leaves = re.findall(pattern, newick)
        return leaves

    @staticmethod
    def complement(sub, whole):
        """Return elements in Whole that are not present in Sub"""
        complement = set(whole) - set(sub)
        return list(complement)

    def get_position(self):
        # 解析树的文件结构
        newick = self.newick  # 调用Tree对象的属性
        position = {}
        ids = 0
        pos = 0
        closed = []
        for word in newick:  # 对树进行解析
            if word == '(':
                ids += 1  # 表示第几个(
                position[ids] = [pos]  # 括号在这行字符串中的位置
            elif word == ')':
                idc = ids
                while idc in closed and idc != 0:
                    idc = idc - 1
                position[idc].append(pos + 1)  # ensure the string incldues the closing parenthesis
                closed.append(idc)
            pos += 1
        # 得到的结果position为：{1: [0, 3092], 2: [1, 511], 3: [2, 267]}
        # 表示第1个括号，前括号在字符串中的位置是0，后括号在字符串中的位置是3092，以此类推
        # 一个括号表示一棵子树
        return position

    def find_child(self, position):
        inspected = []  # 已经检查过的子树
        miss_bval = 0  # 没有支持度的
        for Key in position.keys():
            r_vec = self.newick[position[Key][0]: position[Key][1]]  # 根据P的位置获得子树
            vec = sorted(self.get_leaves(r_vec))  # 获取这棵子树的叶子（基因名）
            covec = sorted(self.complement(vec, self.leaves))  # 获取除这棵子树之外的其它叶子（基因
            if vec not in inspected and covec not in inspected:  # 如果这种分裂方式没有检查过
                my_splits = Split()  # 初始化对象
                my_splits.vecs = [vec, covec]  # 存储这一分裂方式
                exp = re.escape(r_vec) + r'([0-9\.]*:[0-9\.]+)'
                branch_val = re.findall(exp, self.newick)  # 获取这个子树父节点的值
                # (GCA_000716675.1|ORF_01017:0.59570,GCA_000429085.1|ORF_03862:0.71067)0.989:0.43446)
                # BranchVal = ['0.989:0.43446']
                try:
                    my_splits.support = branch_val[0].split(':')[0]  # 支持度：0.989
                    my_splits.branch_length = branch_val[0].split(':')[1]  # 支长：0.43446
                except:
                    miss_bval += 1  # 什么情况下会出现在这情况？
                self.splits.append(my_splits)  # 在Tree中储存这种分裂方式
                inspected.append(vec)
                inspected.append(covec)  # 记录已经讨论过的分裂方式
        return inspected, miss_bval

    def find_leaves(self, inspected, miss_bval):
        for leaf in self.leaves:  # Splits leading to each terminal are included.
            # 对每个基因进行判断
            vec = [leaf]
            covec = sorted(self.complement(vec, self.leaves))
            if vec not in inspected:
                inspected.append(leaf)
                my_splits = Split()
                my_splits.vecs = [vec, covec]
                exp = re.escape(leaf) + r'\:([0-9\.]+)'
                branch_val = re.findall(exp, self.newick)  # 检查每个叶子（基因）的值
                # (GCA_000716675.1|ORF_01017:0.59570,GCA_000429085.1|ORF_03862:0.71067)0.989:0.43446)
                # BranchVal = ['0.59570']
                # BranchVal = ['0.71067']
                try:
                    my_splits.branch_length = branch_val[0]
                except:
                    miss_bval += 1
                self.splits.append(my_splits)

    def split_decomposition(self):
        """Add a list of splits class objects to myPhylo Class object"""
        # Part I. Where we identify matching parenthesis in the newick.
        position = self.get_position()
        # Part II: Where we use string operations to identify components parts of each split.
        inspected, miss_bval = self.find_child(position)
        self.find_leaves(inspected, miss_bval)

    @staticmethod
    def spp_in_list(alist):
        """Return the species from a list of sequence identifiers"""
        spp = []
        for i in alist:
            spp.append(i.split(sep)[0])  # 以|分割，获取物种名，例如：GCA_000716675.1
        return spp

    @staticmethod
    def largest_box(lo_l):
        # 返回一个嵌套列表内最长的那些列表，他们的其它子列表不返回
        nr = []
        for L in lo_l:
            score = 0
            for J in lo_l:
                if set(L).issubset(J):
                    score += 1
            if score < 2:
                nr.append(L)
        return nr

    def orthologs(self, min_taxa, bsupport=0):
        """This function returns populates the list of orthologs in the PhyloClass object"""
        ortho_branch = []
        # if in-paralogs are to be included, update cost value of each terminal.
        for S in self.splits:  # 遍历每个分裂方式
            if S.support in [None, ''] or float(S.support) >= bsupport:  # 如果这种分裂方式不满足要求
                for i_split in S.vecs:  # 遍历分裂的两部分
                    otus = self.spp_in_list(i_split)  # 返回物种的列表
                    if len(set(otus)) == 1 and len(otus) > 1:
                        # 如果这种分裂内只有一个物种，且该物种有多条序列
                        # 即存在内部同源
                        # find splits representing in-paralogs and update costs
                        for leaf in i_split:  # 对该分裂内的基因遍历
                            i_cost = 1.0 / len(otus)
                            if i_cost < self.costs[leaf]:  # 如果存在多个内同源
                                self.costs[leaf] = i_cost  # 更新该基因的分数
                                # Reduce cost value of inparlogue copies in poportion
                                # to the number of inparalogs inplied by this split.
        for S in self.splits:
            if S.support in [None, ''] or float(S.support) >= bsupport:  # 如果这种分裂方式不满足要求
                for i_split in S.vecs:  # 遍历分裂的两部分
                    otus = self.spp_in_list(i_split)  # 返回物种的列表
                    c_count = fsum(self.costs[i] for i in i_split)  # 计算这一分裂方式的分数总和
                    if len(set(otus)) == c_count and c_count >= min_taxa:
                        if i_split not in ortho_branch:  # 如果这种分裂方式没有被检查过
                            ortho_branch.append(i_split)
        ortho_branch = self.largest_box(ortho_branch)  # 取多种分裂方式的最大子树
        self.ortho = ortho_branch  # 添加属性记录

    def get_fasta(self, group):
        sub_dict = dict()
        sub_proteins = []
        for g in self.proteins_list:
            if g.name in group:
                sub_proteins.append(g)
        for pro_ in sub_proteins:
            sub_dict[pro_.name] = pro_.sequence
        return sub_dict

class Trees(list):  # 用于存放所有的DomainTree
    def __init__(self, partitions, prefix=''):
        super().__init__()
        tree_num = 0
        for community in partitions:
            aDomainTree = DomainTree(f'{prefix}tree{tree_num:0>7}', community)
            self.append(aDomainTree)  # 实例化
            tree_num += 1
        self.tree_num = tree_num

    def put_out_fasta(self, out_dir):  # 输出每个tree的fasta
        for tree in self:
            tree: DomainTree
            fasta = '\n'.join([f'>{protein.name}\n{protein.sequence}' for protein in tree.proteins_list])
            fasta_o = FileOperator(name=f'{tree.name}.fa', dir_=out_dir, data=fasta)
            fasta_o.write()

    @staticmethod
    def run_cd_hit(fasta_file, result_dir, thread):
        # cd_hit的命令行
        result_name = os.path.join(result_dir, f'{os.path.basename(fasta_file)}')
        cd_hit_cmd = ' '.join(['cd-hit', '-i', fasta_file, '-o', result_name, '-c', '0.95', '-T', str(thread), '-d', '0'])
        cap = subprocess.Popen(cd_hit_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        cap.communicate()
        if cap.returncode != 0:
            for e in cap.stderr:
                print(e)

    @staticmethod
    def run_mafft(tree_name, fasta_dir, aln_dir, th):
        fasta_file = os.path.join(fasta_dir, f'{tree_name}.fa')
        aln_file = os.path.join(aln_dir, f'{tree_name}.aln')
        mafft_cmd = ' '.join(['mafft', '--anysymbol', '--auto', '--quiet', '--thread', str(th), fasta_file, '>', aln_file])
        cap = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        cap.communicate()
        if cap.returncode != 0:
            for e in cap.stderr:
                print(e)

    @staticmethod
    def cd_hit(fasta_dir, cd_hit_dir, threads):
        # 并行处理
        processes = Pool(processes=threads)
        fasta_files = glob.glob(os.path.join(fasta_dir, '*.fa'))
        for fa in fasta_files:
            processes.apply_async(Trees.run_cd_hit, args=(fa, cd_hit_dir, 4))
        processes.close()
        processes.join()

    def alignment_tree(self, cd_hit_dir, aln_dir, threads):
        processes = Pool(processes=threads)
        ls_name = []
        m_name = []
        for tree in self:
            if tree.gene_num >= 3000 or 100 > tree.gene_num:
                ls_name.append(tree.name)
            if 3000 > tree.gene_num >= 100:
                m_name.append(tree.name)
        for name in m_name:  # python的多任务并行
            processes.apply_async(Trees.run_mafft, args=(name, cd_hit_dir, aln_dir, 4))
        processes.close()
        processes.join()
        for name in ls_name:  # mafft自己的并行
            Trees.run_mafft(name, cd_hit_dir, aln_dir, threads)

    @staticmethod
    def run_fasttree(aln_file, tree_dir):
        tree_file = os.path.join(tree_dir, f'{os.path.basename(aln_file)}.tree')
        fasttree_cmd = ' '.join(['fasttree', aln_file, '>', tree_file])
        # fasttree_cmd = ' '.join(['/opt/miniconda3/bin/fasttree', aln_file, '>', tree_file])
        cap = subprocess.Popen(fasttree_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        cap.communicate()
        if cap.returncode != 0:
            for e in cap.stderr:
                print(e)

    @staticmethod
    def build_tree(aln_dir, tree_dir, threads):
        processes = Pool(processes=threads)
        aln_files = glob.glob(os.path.join(aln_dir, '*.aln'))
        for aln in aln_files:
            processes.apply_async(Trees.run_fasttree, args=(aln, tree_dir))
        processes.close()
        processes.join()

    def get_index(self, cd_reslut_dir):
        # 用于处理cd_hit的结果文件
        for tree in self:
            tree:DomainTree
            index_file = os.path.join(cd_reslut_dir, f'{tree.name}.fa.clstr')
            with open(index_file, 'r') as f:
                contents = f.read().split('>Cluster')[1:]
            index_dic = dict()
            for line in contents:
                gene_list = re.findall(r'>\s*([^,\n]+)', line)
                if len(gene_list) > 1:
                    value = list()
                    for gene in gene_list:
                        if gene.endswith('*'):
                            key = gene.split('...')[0]
                            value.append(key)
                        else:
                            value.append(gene.split('...')[0])
                    index_dic[key] = value
            tree.index = index_dic  # 将检索信息储存在属性中

    def edit_tree(self, tree_raw_dir, tree_dir):
        for tree in self:
            tree:DomainTree
            raw_file = os.path.join(tree_raw_dir, f'{tree.name}.aln.tree')
            tree_file = os.path.join(tree_dir, f'{tree.name}.tree')
            tree_raw = Tree(raw_file)
            if tree.index != None:
                for key, values in tree.index.items():
                    tree_add = Tree()
                    for value in values:
                        new_leaf = tree_add.add_child(name=value, dist=0)  # 构建该节点对应的子树
                    key_node = tree_raw.search_nodes(name=key)[0]
                    key_node.add_child(tree_add)
                tree_raw.write(outfile=tree_file)

    def ana_newick(self, tree_dir, out_dir, og_dir, max_genome):  # 这部分是调用类里的函数，暂时没有使用
        # 用于读取树和输出树
        ort_list = ''
        Total = 0
        no_tree = 0
        for tree in self:
            tree: DomainTree
            name = tree.name  # 获取文件名
            tree_file = os.path.join(tree_dir, f'{os.path.basename(tree.name)}.aln.tree')
            try:
                with open(tree_file, 'r') as f:  # 读取tree文件
                    newick = f.read()
                tree.get_newick(newick)  # 创建对象的树属性
                tree.split_decomposition()  # 拆分子树
                tree.orthologs(max_genome)  # 查找最大子树
                ort_num = 0
                for group in tree.ortho:  # 拆分得到的每棵子树，list
                    g_seq = ''
                    f_name = f'{name}_{ort_num},'
                    g = ','.join(group).strip(',')  # 输出子树的基因名
                    ort_list += f_name + g + '\n'  # 所有子树的文件
                    group_seq_dict = tree.get_fasta(group)
                    for id, seq in group_seq_dict.items():
                        g_seq += f'>{str(id)}\n{str(seq)}\n'
                    g_seq_file = os.path.join(og_dir, f'{name}_{ort_num}.fasta')
                    with open(g_seq_file, 'w') as f:
                        f.write(g_seq)
                    Total += 1
                    ort_num += 1
            except FileNotFoundError:
                no_tree += 1
        orthogroups_file = os.path.join(out_dir, 'orthogroups.csv')
        with open(orthogroups_file, 'w') as f:
            f.write(ort_list)
        message(text=f'Total orthogroups found: {Total}', label='Infomation')
        message(text=f'{no_tree} tree_files cannot be found.', label='Infomation')
