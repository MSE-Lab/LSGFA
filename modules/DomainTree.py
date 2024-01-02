#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2023-12-13 15:43
# @Author  : zhaoyu
# @File    : DomainTree.py
import os
from modules.utils import *
import subprocess
from multiprocessing import Pool


class DomainTree:
    def __init__(self, name, proteins_list):  # 一个DomainTree存放一个community的tree
        self.name = name
        self.proteins_list = proteins_list
        self.gene_num = len(proteins_list)

    def write_fasta(self, out_dir):  # 输出这棵树的fasta文件
        fasta = '\n'.join([f'>{protein.name}\n{protein.sequence}' for protein in self.proteins_list])
        fasta_o = FileOperator(name=f'{self.name}.fa', dir_=out_dir, data=fasta)
        fasta_o.write()

    def alignment(self, fasta_file, aln_file):
        if self.gene_num > 1:
            mafft_cmd = ' '.join(['mafft', '--anysymbol', '--auto', '--quiet', fasta_file, '>', aln_file])
            cap = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            cap.communicate()
            if cap.returncode != 0:
                # some errors happened
                for e in cap.stderr:
                    print(e)

    @staticmethod
    def fasttree(aln_file, out_file):
        fasttree_cmd = ' '.join(['/opt/miniconda3/bin/fasttree', aln_file, '>', out_file])
        cap = subprocess.Popen(fasttree_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        cap.communicate()
        if cap.returncode != 0:
            # some errors happened
            for e in cap.stderr:
                print(e)


class Trees(list):  # 用于存放所有的DomainTree
    def __init__(self, partitions):
        super().__init__()
        tree_num = 0
        for community in partitions:
            aDomainTree = DomainTree(f'tree{tree_num:0>7}', community)
            self.append(aDomainTree)  # 实例化
            tree_num += 1
        self.tree_num = tree_num

    def put_out_fasta(self, out_dir, threads):  # 输出每个tree的fasta
        for tree in self:
            tree.write_fasta(out_dir)
        # processes = Pool(processes=threads)
        # for tree in self:
        #     processes.apply_async(tree.write_fasta, args=(out_dir,))
        # processes.close()
        # processes.join()

    def alignment_tree(self, fasta_dir, aln_dir, threads):
        processes = Pool(processes=threads)
        for tree in self:
            tree: DomainTree
            fasta_file = os.path.join(fasta_dir, f'{tree.name}.fa')
            aln_file = os.path.join(aln_dir, f'{tree.name}.aln')
            processes.apply_async(tree.alignment, args=(fasta_file, aln_file,))
        processes.close()
        processes.join()

    def build_tree(self, aln_dir, tree_dir, threads):
        processes = Pool(processes=threads)
        for tree in self:
            tree: DomainTree
            if tree.gene_num > 1:
                aln_file = os.path.join(aln_dir, f'{tree.name}.aln')
                tree_file = os.path.join(tree_dir, f'{tree.name}.tree')
                processes.apply_async(tree.fasttree, args=(aln_file, tree_file,))
        processes.close()
        processes.join()
