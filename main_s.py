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

global OUT_DIR, PFAM_DIR, GRAPH_DIR, ALN_DIR, TREE_DIR, \
    QUERY_DIR, OG_DIR, THREADS, SEQ_INFO, INFLATION


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
    group1.add_option(
        '-I', '--inflation', type=str, dest='inflation_co', default="1.5",
        help='MCL inflation, default: 1.5')
    opt.add_option_group(group0)
    opt.add_option_group(group1)
    options, args = opt.parse_args()
    input_dir = options.input_dir
    output_dir = options.output_dir
    extension = options.extension
    search_method = options.search_method
    e_values = options.e_values
    search_threads = options.search_threads
    inflation_co = options.inflation_co
    re_run = options.re_run
    parameters_dict = dict(
        input_dir=input_dir, output_dir=output_dir, extension=extension, search_method=search_method, e_values=e_values,
        search_threads=search_threads, inflation_co=inflation_co, re_run=re_run)
    return parameters_dict


def make_working_dir(re_run):
    global PFAM_DIR, GRAPH_DIR, ALN_DIR, TREE_DIR, QUERY_DIR, OG_DIR
    PFAM_DIR = os.path.join(OUT_DIR, 'pfam')
    GRAPH_DIR = os.path.join(OUT_DIR, 'graph')
    ALN_DIR = os.path.join(OUT_DIR, 'alignment')
    TREE_DIR = os.path.join(OUT_DIR, 'tree')
    QUERY_DIR = os.path.join(OUT_DIR, 'query')
    OG_DIR = os.path.join(OUT_DIR, 'orthogroups')
    if re_run:
        try:
            shutil.rmtree(GRAPH_DIR)
            shutil.rmtree(ALN_DIR)
            shutil.rmtree(TREE_DIR)
            shutil.rmtree(QUERY_DIR)
            shutil.rmtree(OG_DIR)
            os.remove(os.path.join(OUT_DIR, 'orthogroups.csv'))
        except FileNotFoundError:
            message(text='star re_run')
    [os.makedirs(dir_, exist_ok=True) for dir_ in [PFAM_DIR, GRAPH_DIR, QUERY_DIR, ALN_DIR, TREE_DIR, OG_DIR]]


@time_used(f'[{timing()}]Pfam annotation')
def pfam_annotation(input_genomes):
    pp = Panproteome(input_genomes)
    message(text=f'Start with PFAM annotation ...', label='PROCESS')
    pp.put_pfam_file(threads=100, outdir=PFAM_DIR)   # pfam注释
    return pp


@time_used(f'[{timing()}]Building pfam graph')
def build_pfam_graph(pp, max_genome):
    message(text=f'Start building the graph ...', label='PROCESS')
    pfam_graph = PGraph(pp)   # 初始化
    basic_graph = pfam_graph.generate_graph()
    basic_graph.write_gml(os.path.join(GRAPH_DIR, 'pfam_graph.gml'))
    basic_graph.write_ncol(os.path.join(GRAPH_DIR, 'pfam_graph.txt'), names='domain_type')
    partitions = pfam_graph.la_find_partition()   # 社区发现
    # pfam_graph.partition_p_og(partition, max_genome, out_dir)
    return partitions


def upho_tree(input_dir, max_genome):
    tree_file = glob.glob(os.path.join(TREE_DIR, '*.aln.tree'))
    orthogroups_file = os.path.join(OUT_DIR, 'orthogroups.csv')
    OrtList = open(orthogroups_file, 'a')
    upho.upho_main(tree_file, max_genome, OrtList)
    OrtList.close()

    out_file = os.path.join(OUT_DIR, 'OG_no_subsets.txt')
    upho.No_OG_subsets(orthogroups_file, out_file)

    no_file = os.path.join(OUT_DIR, 'UPhO_nr_orthogroups.csv')
    upho.No_Same_OG_Intesec(out_file, no_file)
    os.remove(out_file)

    message(text="Proceeding to create a fasta file for each ortholog")
    upho.gfr_main(query=no_file, outdir=OG_DIR, prefix='upho', input_path=input_dir)


@time_used(f'[{timing()}]Building Domain Type Tree')
def build_domain_tree(input_dir, partitions, threads, max_genome):
    message(text=f'Start building Domain Type Tree ...', label='PROCESS')
    trees = Trees(partitions)
    trees.put_out_fasta(QUERY_DIR)
    message(text='Put_out_fasta Done.', label='PROCESS')
    trees.alignment_tree(QUERY_DIR, ALN_DIR, threads)
    message(text='Alignment Done.', label='PROCESS')
    trees.build_tree(ALN_DIR, TREE_DIR, threads)
    message(text='Build tree Done.', label='PROCESS')

    # trees.ana_newick(TREE_DIR, OUT_DIR, OG_DIR, max_genome)
    # message(text='ana_newick Done.', label='PROCESS')

    upho_tree(input_dir, max_genome)


# @time_used(f'[{timing()}]Homology searching')
# def og_searching(p_graph):
#     db_cmds, search_cmds = p_graph.homology_search_commands(query_path=QUERY_DIR, db_path=DB_DIR,
#                                                             res_path=RES_DIR, seq_info=SEQ_INFO)
#     if None not in db_cmds:
#         db_cmds = CallCmd(cmd_list=db_cmds, process_info="Building homology searching database", parallel=True,
#                           threads=THREADS)
#         db_cmds.parallel_process()
#     search_cmd = CallCmd(search_cmds, process_info="Homology searching", parallel=True, threads=THREADS)
#     search_cmd.parallel_process()


# @time_used(f"[{timing()}]MCL")
# def running_mcl(p_graph):
#     abc_name = os.path.join(OUT_DIR, "mcl.abc")
#     out_name = os.path.join(OUT_DIR, "mcl.cluster.txt")
#     p_graph.mcl_abc(res_dir=RES_DIR, abc_file_name=abc_name, threads=THREADS, out_name=out_name,
#                     inflation=INFLATION)


@time_used(f"[{timing()}]Writing orthogroups")
def write_og_files(mcl_cluster_file, og_path, max_g, seq_info):
    clusters = Clusters(cluster_file_name=mcl_cluster_file, max_genome=max_g)  # 初始化
    clusters.write_og(seq_path=og_path, seq_info=seq_info)


def get_resolution_profile(simplify_pfam_graph, max_genome, outfile_name):
    optimiser = la.Optimiser()
    profile = optimiser.resolution_profile(simplify_pfam_graph, la.CPMVertexPartition,
                                           resolution_range=(0, 1),
                                           weights='weight')
    result = ''
    for p in profile:  # 每个参数的结果循环
        p_sub = p.subgraphs()  # 生成网络
        p_OG = 0
        for p_s in p_sub:  # 每个社区计算
            # mode = '  '.join(p_s.vs['name'])
            genes_in_community = [item for sublist in p_s.vs['genes'] for item in sublist]
            genomes_list = [n.split('|')[0] for n in genes_in_community]
            genomes_set = set(genomes_list)  # 判断该community是否是核心
            if len(genomes_set) >= max_genome:
                p_OG += min(Counter(genomes_list).values())  # 获取每个PG的大小
            else:
                p_OG += 0
        result += f'{p.resolution_parameter}\t{len(p.subgraphs())}\t{p_OG}\n'
    with open(outfile_name, 'w') as f:
        title = '#resolution_parameter\tcommunity_num\tPotential_og\n'
        f.write(title)
        f.write(result)


@time_used(f"[{timing()}]Whole processing Done")
def main():
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters['input_dir']
    global OUT_DIR, THREADS, SEQ_INFO, INFLATION
    OUT_DIR = parameters['output_dir']
    THREADS = parameters['search_threads']
    search_method = parameters['search_method']
    INFLATION = parameters['inflation_co']
    re_run = parameters['re_run']
    make_working_dir(re_run)
    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    # pfam graph construction
    pp = pfam_annotation(input_genomes_dir)
    partitions = build_pfam_graph(pp, max_genome)
    # 每个partition画树
    build_domain_tree(input_genomes_dir, partitions, THREADS, max_genome)


if __name__ == '__main__':
    main()
