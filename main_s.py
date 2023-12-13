#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2023-07-11 17:27
# @Author  : zhaoyu
# @File    : main_s.py
from optparse import OptionGroup, OptionParser
from modules.homology_search import *
from modules.panproteome import *
import leidenalg as la

global OUT_DIR, DB_DIR, RES_DIR, QUERY_DIR, OG_DIR, THREADS, SEQ_INFO, INFLATION


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
    parameters_dict = dict(
        input_dir=input_dir, output_dir=output_dir, extension=extension, search_method=search_method, e_values=e_values,
        search_threads=search_threads, inflation_co=inflation_co)
    return parameters_dict


def make_working_dir():
    global DB_DIR, RES_DIR, QUERY_DIR, OG_DIR
    DB_DIR = os.path.join(OUT_DIR, 'db')
    RES_DIR = os.path.join(OUT_DIR, 'res')
    QUERY_DIR = os.path.join(OUT_DIR, 'query')
    OG_DIR = os.path.join(OUT_DIR, 'Orthogroups')
    [os.makedirs(dir_, exist_ok=True) for dir_ in [DB_DIR, RES_DIR, QUERY_DIR, OG_DIR]]


@time_used(f'[{timing()}]Building pfam graph')
def build_pfam_graph(input_genomes):
    pp = Panproteome(input_genomes)
    pp.put_pfam_file(threads=100, outdir=OUT_DIR)   # pfam注释
    simplify_pfam_graph = pp.make_pfam_graph(OUT_DIR)    # 生成graph
    genome_info = pp.make_sequences_info()  # 获取蛋白质id及其对应的序列
    return simplify_pfam_graph, genome_info


@time_used(f'[{timing()}]Homology searching')
def og_searching(p_graph):
    db_cmds, search_cmds = p_graph.homology_search_commands(query_path=QUERY_DIR, db_path=DB_DIR,
                                                            res_path=RES_DIR, seq_info=SEQ_INFO)
    if None not in db_cmds:
        db_cmds = CallCmd(cmd_list=db_cmds, process_info="Building homology searching database", parallel=True,
                          threads=THREADS)
        db_cmds.parallel_process()
    search_cmd = CallCmd(search_cmds, process_info="Homology searching", parallel=True, threads=THREADS)
    search_cmd.parallel_process()


@time_used(f"[{timing()}]MCL")
def running_mcl(p_graph):
    abc_name = os.path.join(OUT_DIR, "mcl.abc")
    out_name = os.path.join(OUT_DIR, "mcl.cluster.txt")
    p_graph.mcl_abc(res_dir=RES_DIR, abc_file_name=abc_name, threads=THREADS, out_name=out_name,
                    inflation=INFLATION)


@time_used(f"[{timing()}]Writing orthogroups")
def write_og_files(mcl_cluster_file, og_path, max_g, seq_info):
    clusters = Clusters(cluster_file_name=mcl_cluster_file, max_genome=max_g)  # 初始化
    clusters.write_og(seq_path=og_path, seq_info=seq_info)


@time_used(f"[{timing()}]Find communities")
def la_find_partition(simplify_pfam_graph, max_genome, outfile_name):
    partition = la.find_partition(simplify_pfam_graph, partition_type=la.CPMVertexPartition, weights='weight',
                                  resolution_parameter=0.9)
    result = ''
    scc_num = 0
    for community in partition:
        community_subgraph = simplify_pfam_graph.subgraph(community)
        mode_in_community = community_subgraph.vs['name']
        # mode = '  '.join(mode_in_community)
        genes_in_community = [item for sublist in community_subgraph.vs['genes'] for item in sublist]
        genomes_list = [n.split('|')[0] for n in genes_in_community]
        genomes_set = set(genomes_list)  # 判断该community是否是核心
        if len(genomes_set) >= max_genome:
            p_OG = min(Counter(genomes_list).values())  # 获取每个PG的大小
        else:
            p_OG = 0
        result += f'SCC{scc_num:0>7}\t{len(mode_in_community)}\t{len(genes_in_community)}\t{p_OG}\n'
        scc_num += 1
    with open(outfile_name, 'w') as f:
        title = '#S_ID\tMode_num\tS_size\tPotential_og\n'
        f.write(title)
        f.write(result)

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
    # make_working_dir()
    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    print(f'genomes Numbers: {max_genome}')

    # pfam graph construction
    s_pg, s_info = build_pfam_graph(input_genomes_dir)
    # la_find_partition(s_pg, max_genome, os.path.join(OUT_DIR, 's_pg_count_og.txt'))
    get_resolution_profile(s_pg, max_genome, os.path.join(OUT_DIR, 'resolution_profile.txt'))

    # test_graph, SEQ_INFO = build_pfam_graph(input_genomes_dir)
    # homology searching
    # p_graph = PfamG(test_graph, method=search_method)
    # og_searching(p_graph)
    # mcl
    # running_mcl(p_graph)
    # writing OGs
    # out_name = os.path.join(OUT_DIR, "mcl.cluster.txt")
    # write_og_files(mcl_cluster_file=out_name, og_path=OG_DIR,
    #                max_g=max_genome,
    #                seq_info=SEQ_INFO)


if __name__ == '__main__':
    main()
