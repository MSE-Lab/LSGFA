#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2023-07-11 17:27
# @Author  : zhaoyu
# @File    : main_s.py
import os
from optparse import OptionGroup, OptionParser
from modules.homology_search import *
from modules.panproteome import *
from modules.DomainTree import *
import modules.upho as upho
import leidenalg as la
import shutil

global OUT_DIR, THREADS


def get_parameters():
    usage = "LSGFA.py -i <input dir> <Options>"
    opt = OptionParser(usage=usage)
    group0 = OptionGroup(opt, "General options")
    group1 = OptionGroup(opt, "homology search")
    group0.add_option(
        '-i', '--in', type=str, dest='input_dir', default=False,
        help='specify the directory including all genome files')
    group0.add_option(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f"specify a output directory default: {os.getcwd()}")
    group0.add_option(
        '-x', '--extension', type=str, dest='extension', default='fna',
        help=f"Extension of genome files; default: fna")
    group0.add_option(
        '-f', '--force', type=str, dest='re_run', default=False,
        help=f"Re-perform the homology search; default: False")
    group1.add_option(
        '-s', '--search_method', dest='search_method', choices=['diamond', 'mmseqs', 'blastp'], default='diamond',
        help='Homologs searching methods: blastp, mmseqs, diamond. Both mmseqs and diamond are sensitive mode ('
             'default: diamond)')
    group1.add_option(
        '-e', '--evalue', type=float, dest='e_values', default=1e-5,
        help="cut-off of E value")
    group1.add_option(
        '-t', '--threads', type=int, dest='search_threads', default=8,
        help='Homologs searching threads. default: 8')
    opt.add_option_group(group0)
    opt.add_option_group(group1)
    options, args = opt.parse_args()
    input_dir = options.input_dir
    output_dir = options.output_dir
    extension = options.extension
    search_method = options.search_method
    e_values = options.e_values
    search_threads = options.search_threads
    re_run = options.re_run
    parameters_dict = dict(
        input_dir=input_dir, output_dir=output_dir, extension=extension, search_method=search_method, e_values=e_values,
        search_threads=search_threads, re_run=re_run)
    return parameters_dict


def make_working_dir(re_run):
    output_dirs = ['pfam', 'graph', 'alignment', 'tree_raw', 'tree',
                   'cd_hit_result', 'query', 'orthogroups', 'none_pfam',
                   'none_pfam/fa_file']
    if re_run:
        for dir_name in output_dirs:
            try:
                shutil.rmtree(os.path.join(OUT_DIR, dir_name))
            except FileNotFoundError:
                pass
        try:
            os.remove(os.path.join(OUT_DIR, 'orthogroups.csv'))
        except FileNotFoundError:
            pass
    [os.makedirs(dir_, exist_ok=True) for dir_ in output_dirs]


def upho_tree(input_dir, max_genome, tree_dir, og_dir):
    tree_file = glob.glob(os.path.join(tree_dir, '*.aln.tree'))  # 读取tree文件
    orthogroups_file = os.path.join(OUT_DIR, 'orthogroups.csv')  # 用于存储子树的内容
    OrtList = open(orthogroups_file, 'a')  # 写入
    upho.upho_main(tree_file, max_genome, OrtList)
    OrtList.close()

    out_file = os.path.join(OUT_DIR, 'OG_no_subsets.txt')  # 去重复
    upho.No_OG_subsets(orthogroups_file, out_file)

    no_file = os.path.join(OUT_DIR, 'UPhO_nr_orthogroups.csv')  # 去冗余和内部同源
    upho.No_Same_OG_Intesec(out_file, no_file)
    os.remove(out_file)

    message(text="Proceeding to create a fasta file for each ortholog")  # 输出每个子树的fasta文件
    upho.gfr_main(query=no_file, outdir=og_dir, prefix='upho', input_path=input_dir)


@time_used(f'[{timing()}]All to all blast for none_PFAM.')
def aTa_blast(Ngroup, query, none_pfam):
    blast_db = os.path.join(none_pfam, 'db')
    res_dir = os.path.join(none_pfam, 'res')
    [os.makedirs(dir_, exist_ok=True) for dir_ in [blast_db, res_dir]]
    db_cmds, search_cmds = Ngroup.homology_search(query_dir=query, db_dir=blast_db,
                                                  res_dir=res_dir, threads=THREADS)
    if None not in db_cmds:
        db_cmds = CallCmd(cmd_list=db_cmds, process_info="Building database", threads=THREADS)
        db_cmds.processing()
    search_cmd = CallCmd(search_cmds, process_info="Homology searching", threads=THREADS)
    search_cmd.processing()
    Ngroup.build_homology_graph(os.path.join(res_dir, f'{Ngroup.name}.txt'))
    none_partition = Ngroup.get_partition_genes()
    return none_partition


@time_used(f"[{timing()}]Whole processing Done")
def main():
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters['input_dir']
    global OUT_DIR, THREADS
    OUT_DIR = parameters['output_dir']
    THREADS = parameters['search_threads']
    search_method = parameters['search_method']
    re_run = parameters['re_run']

    pfam_dir = os.path.join(OUT_DIR, 'pfam')
    graph_dir = os.path.join(OUT_DIR, 'graph')
    aln_dir = os.path.join(OUT_DIR, 'alignment')
    tree_raw_dir = os.path.join(OUT_DIR, 'tree_raw')
    tree_dir = os.path.join(OUT_DIR, 'tree')
    query_dir = os.path.join(OUT_DIR, 'query')
    og_dir = os.path.join(OUT_DIR, 'orthogroups')
    no_pfam = os.path.join(OUT_DIR, 'none_pfam')
    no_query = os.path.join(no_pfam, 'fa_file')
    cd_hit_dir = os.path.join(no_pfam, 'cd_hit_result')

    make_working_dir(re_run)
    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    pp = Panproteome(input_genomes_dir)  # 初始化
    message(text=f'Start with PFAM annotation ...', label='PROCESS')
    pp.put_pfam_file(threads=100, outdir=pfam_dir)  # pfam注释

    message(text=f'Start building the graph ...', label='PROCESS')
    pfam_graph = PGraph(pp, no_query)  # 初始化，提出需要做blast的文件
    pfam_graph.generate_graph()  # 构建网络
    pfam_graph.put_graph_file(graph_dir)  # 输出网络相关文件
    partitions = pfam_graph.la_find_partition()  # 社区发现
    # blast后再画树的部分
    message(text=f'Start building None_pfam Type Tree ...', label='PROCESS')
    for fa_file in glob.glob(os.path.join(no_query, '*.fa')):  # 循环对需要blast的文件处理
        fa_blast = Ngroup(fa_file, method=search_method)
        fa_partitions = aTa_blast(fa_blast, no_query)
        partitions.extend(fa_partitions)

    # 画树
    message(text=f'Start building Domain Type Tree ...', label='PROCESS')
    trees = Trees(partitions)
    trees.put_out_fasta(query_dir)
    message(text='Put_out_fasta Done.', label='PROCESS')
    # 在比对之前先去冗余
    trees.cd_hit(query_dir, cd_hit_dir, THREADS)
    # 画树
    trees.alignment_tree(query_dir, aln_dir, THREADS)
    message(text='Alignment Done.', label='PROCESS')
    trees.build_tree(aln_dir, tree_raw_dir, THREADS)
    message(text='Build tree Done.', label='PROCESS')
    # 编辑树文件，添加叶子
    trees.edit_tree(tree_raw_dir, tree_dir)

    # upho拆分树
    upho_tree(query_dir, max_genome, tree_dir, og_dir)
    message(text='Analyse tree Done.', label='PROCESS')


if __name__ == '__main__':
    main()
