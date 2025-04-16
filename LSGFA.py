#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Time    ：2024/10/30
# @Author  ：zhaoyu
# @File    ：LSGFA.py
import os.path
from modules.panproteome import *
from modules.build_graph import *
import shutil
import argparse
from modules.homology_search import *
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import modules.none_mapping as mapping


def get_parameters():
    parser = argparse.ArgumentParser(
        description='This program can get the Pfam Network.',
        formatter_class=argparse.RawTextHelpFormatter)

    # 基本参数
    parser.add_argument(
        '-i', '--in', type=str, dest='input_dir',
        help='The directory including all genome files')
    parser.add_argument(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f'Specify a output directory default: {os.getcwd()}')

    # 重做参数
    parser.add_argument('-f', action='store_true',
                        help='Re-perform whole process(including pfam annotation)')
    parser.add_argument('-fg', action='store_true',
                        help='Re-perform the graph search(not including pfam annotation)')
    parser.add_argument("-fb", action="store_true",
                        help="Re-perform the homology search(blast)")

    # 停断参数
    parser.add_argument('--pfam', action='store_true',
                        help='Stop at pfam annotation)')
    parser.add_argument('--pg', action='store_true',
                        help='Stop at pfam graph search')

    # 个性化参数
    parser.add_argument('-r', action='store_true',
                        help='Run deduplication')
    parser.add_argument(
        '-t', '--threads', type=int, dest='threads', default=8,
        help='Hmmscan threads. default: 8')
    parser.add_argument(
        '-db', dest='pfam_db', default=os.path.join(os.getcwd(), 'modules', 'database', 'Pfam-A.hmm'),
        help=f"Pfam database path. default: {os.path.join(os.getcwd(), 'modules', 'database', 'Pfam-A.hmm')}")
    parser.add_argument(
        '-search',
        choices=['hmmscan', 'mmseqs-search'], default='hmmscan',  # 只允许这几个选项
        help='Select a method for blasting. default=hmmscan')
    parser.add_argument(
        '-ssn',
        choices=['1', '2', '3'], default='3',  # 只允许这几个选项
        help='Select a method for build sequence similarity network (SSN).\n'
             'default = 3\n'
             '1 = Reciprocal best hit (RBH)\n'
             '2 = Specify a reciprocal hit above the threshold (identity >= 40).\n'
             '3 = reciprocal hits above the minimum reciprocal hit threshold.\n')
    parser.add_argument(
        '-blast',
        choices=['diamond', 'mmseqs-search'], default='mmseqs-search',  # 只允许这几个选项
        help='Select a method for blasting. default=mmseqs-search')
    parser.add_argument(
        '-id', dest='identity', type=int, default=40,
        help='The identity of homology search [0-100], default = 40.')
    parser.add_argument(
        '-c', dest='coverage', type=int, default=50,
        help='The coverage of homology search [0-100], default = 50.')
    parser.add_argument(
        '-e', dest='evalue', default=1e-5,
        help='The coverage of homology search [0-1], default = 1e-5.')
    parser.add_argument(
        '-inflation', dest='inflation', default=1.5,
        help='Inflation (varying this parameter affects granularity) [1.2-5.0], default = 1.5.')
    parser.add_argument(
        '-partition_type',
        choices=['1', '2', '3', '4', '5'], default='4',  # 只允许这几个选项
        help='Select a method for la.partition_type.\n'
             'default = 4\n'
             '1 = ModularityVertexPartition\n'
             '2 = RBConfigurationVertexPartition\n'
             '3 = RBERVertexPartition\n'
             '4 = CPMVertexPartition\n'
             '5 = SurpriseVertexPartition\n')
    parser.add_argument(
        '-rp', dest='resolution_parameter', default=0.9,
        help='[0-1.0], default = 0.9.\n'
             'Some methods accept resolution parameters,\n'
             'such as RBConfigurationVertexPartition, RBERVertexPartition and CPMVertexPartition. \n'
             'The larger the resolution_parameter, the more subgraphs will be obtained.')

    args = parser.parse_args()  # general options
    return args


def make_working_dir(output_dirs, delate=True):
    if delate:
        for dir_name in output_dirs:
            try:
                shutil.rmtree(dir_name)
            except FileNotFoundError:
                pass
    # 重新创建目录
    [os.makedirs(dir_, exist_ok=True) for dir_ in output_dirs]


def work_flow(input_file, blast_out_dir, threads, sub_cc_dir,  method, num, inflation, ssn):
    result_files = []
    agroup = DomainGroup(input_file)  # cc文件
    if len(agroup) <= 2:  # 如果只有两条序列，就不做blast和mcl了
        abc_file = None
        shutil.copyfile(input_file, os.path.join(sub_cc_dir, f'{agroup.name}_1.faa'))
    else:
        abc_file = os.path.join(blast_out_dir, f'{agroup.name}.abc')

        if not os.path.exists(abc_file):
            blast_res = os.path.join(blast_out_dir, f'{agroup.name}.txt')
            if method == 'diamond':
                agroup.make_db(blast_out_dir, threads)  # make db
                split_files = agroup.split_file(blast_out_dir)
                for split_file in split_files:
                    result_file = agroup.homology_abc(split_file, blast_out_dir, threads, method, num)
                    result_files.append(result_file)
                agroup.merge_files(blast_res, result_files)
                if len(split_files) != 1:
                    for file in split_files + result_files:
                        os.remove(file)
            elif method == 'mmseqs-search':
                result_file = agroup.homology_abc(input_file, blast_out_dir, threads, method, num)
        try:
            agroup.filter_abc(blast_res, abc_file, ssn)
            agroup.mcl_identity(abc_file, blast_out_dir, threads, inflation=inflation)  # 处理blast结果
            agroup.mcl_cc_file(sub_cc_dir)  # mcl聚类
        except FileNotFoundError as e:
            print(f"Error: {e}. \n"
                  f"Some error happend, there is no {result_file}.")


# def work_flow(input_file, blast_out_dir, threads, sub_cc_dir, identity, coverage, method, num, inflation):
#     result_files = []
#     agroup = DomainGroup(input_file)  # cc文件
#     if len(agroup) <= 2:  # 如果只有两条序列，就不做blast和mcl了
#         result_file_path = None
#         shutil.copyfile(input_file, os.path.join(sub_cc_dir, f'{agroup.name}_1.faa'))
#     else:
#         result_file_path = os.path.join(blast_out_dir, f'{agroup.name}.abc')
#
#         if not os.path.exists(result_file_path):
#             if method == 'diamond':
#                 agroup.make_db(blast_out_dir, threads)  # make db
#                 split_files = agroup.split_file(blast_out_dir)
#                 for split_file in split_files:
#                     result_file = agroup.homology_abc(split_file, blast_out_dir, threads, method, num, identity, coverage)
#                     result_files.append(result_file)
#                 agroup.merge_files(result_file_path, result_files)
#                 if len(split_files) != 1:
#                     for file in split_files + result_files:
#                         os.remove(file)
#             elif method == 'mmseqs-search':
#                 result_file_path = agroup.homology_abc(input_file, blast_out_dir, threads, method, num, identity, coverage)
#         try:
#             agroup.mcl_identity(result_file_path, blast_out_dir, threads, inflation=1.5)  # 处理blast结果
#         except FileNotFoundError as e:
#             print(f"Error: {e}. \n"
#                   f"Some error happend, there is no {result_file}.")
#         agroup.mcl_cc_file(sub_cc_dir)  # 建立rbh


def parallel_workflow(faa_files, blast_out_dir, threads, sub_cc_dir, method, inflation, ssn, num=None):
    total_tasks = len(faa_files)  # 添加进度条
    make_working_dir([blast_out_dir, sub_cc_dir])
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(work_flow, faa, blast_out_dir, str(threads), sub_cc_dir,
                                   method, num, inflation, ssn): faa for faa in faa_files}
        with tqdm(total=total_tasks) as pbar:
            for future in as_completed(futures):
                try:
                    future.result()  # This will raise an exception if the function raised one
                except Exception as e:
                    print(f"Error processing {futures[future]}: {e}")
                finally:
                    pbar.update(1)
        # 等待所有任务完成（虽然 as_completed 已经做了这个，但我们可以显式地处理）
        for future in futures.keys():
            future.result()  # 确保每个任务都已完成，这里也会捕获异常

    # 定义要删除的文件和目录的条件
    delete_conditions = [lambda x: x.startswith('tmp'),
                         lambda x: x.endswith('txt'),
                         lambda x: x.endswith('dmnd')]

    for item in os.listdir(blast_out_dir):
        item_path = os.path.join(blast_out_dir, item)
        if any(condition(item) for condition in delete_conditions):
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
            elif os.path.isfile(item_path):
                os.remove(item_path)


class GeneGroup:
    """
    GeneGroup作为存放RBH后获得的子家族
    每个GeneGroup都是一个gene的cluster
    """

    def __init__(self, name: str = "", marker: str = '', species=None):
        self.name = name  # sub_cc_id
        self.marker = marker
        self.species = species


def define_marker(percent):
    marker = ''
    if 0.99 <= percent <= 1:
        marker = 'Core'
    elif 0.95 <= percent < 0.99:
        marker = 'Soft_core'
    elif 0.15 <= percent < 0.95:
        marker = 'Shell'
    elif 0 <= percent < 0.15:
        marker = 'Cloud'
    return marker


def cc_list(cc_dir, species_num, hom_dir):
    # 获取路径中的faa文件，处理为cc:gene_num的字典
    members = sorted(glob.glob(os.path.join(cc_dir, '*.faa')))
    object_list = []
    for member in members:
        cc_gene = gen_seqs_with_headers(member, extract_ids=True)  # 该cc内的基因
        cc_name = os.path.basename(member).split('.')[0]
        cc_species = set([i.split('|')[0] for i in cc_gene])
        marker_ = define_marker(len(cc_species) / int(species_num))
        aGeneGroup = GeneGroup(name=cc_name, marker=marker_, species=cc_species)
        object_list.append(aGeneGroup)

        # 将结果添加到文件中
        cc_result = f'{cc_name}\n{",".join(cc_gene)}\n'
        with open(os.path.join(hom_dir, 'sub_cc_list.txt'), 'a') as f:
            f.write(cc_result)
    return object_list


# 提出不同marker的cc
def put_out_file(object_list, out_dir):
    marker_dict = {'Core': [], 'Soft_core': [], 'Shell': [], 'Cloud': []}
    for obj in object_list:
        obj: GeneGroup
        marker_dict[obj.marker].append(obj.name)
    for k, v in marker_dict.items():
        ccs = '\n'.join(v)
        with open(os.path.join(out_dir, f'{str(k)}_genes_list.txt'), 'w') as f:
            f.write(ccs)
    summary = f'Core genes (99% <= strains <= 100%)\t{len(marker_dict["Core"])}\n' \
              f'Soft core genes (95% <= strains < 99%)\t{len(marker_dict["Soft_core"])}\n' \
              f'Shell genes (15% <= strains < 95%)\t{len(marker_dict["Shell"])}\n' \
              f'Cloud genes (0% <= strains < 15%)\t{len(marker_dict["Cloud"])}\n'
    with open(os.path.join(out_dir, 'summary.txt'), 'w') as f:
        f.write(summary)


@time_used(f"[{timing()}]Whole processing Done")
def main():
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters.input_dir
    OUT_DIR = parameters.output_dir
    THREADS = parameters.threads

    identity = parameters.identity
    coverage = parameters.coverage
    evalue = parameters.evalue
    inflation = float(parameters.inflation)
    partition_type = parameters.partition_type
    resolution_parameter = float(parameters.resolution_parameter)
    ssn = parameters.ssn

    pfam_dir = os.path.join(OUT_DIR, 'pfam')  # 存放pfam注释结果
    graph_dir = os.path.join(OUT_DIR, 'graph')  # 按pfam建graph的信息
    query_dir = os.path.join(OUT_DIR, 'graph_cc')  # pfam分出来
    no_pfam = os.path.join(OUT_DIR, 'none_pfam')  # 未注释到pfam信息的处理结果
    homo_dir = os.path.join(OUT_DIR, 'homology_search')  # 使用blast对cc进一步拆分
    pangenome_dir = os.path.join(OUT_DIR, 'pangenome')  # pangenome信息

    # 判断是否进行重做
    if parameters.fg:  # 从PGraph重做
        message(text=f"Re-perform the graph search(not including annotation).", label='PROCESS')
        output_dirs = [graph_dir, query_dir, no_pfam,
                       homo_dir, pangenome_dir]
        make_working_dir(output_dirs)
    if parameters.fb:
        message(text=f"Re-perform the homology search(blast).", label='PROCESS')
        output_dirs = [homo_dir, pangenome_dir]
        make_working_dir(output_dirs)
    if parameters.f:
        message(text=f"Re-perform whole process(including pfam annotation).", label='PROCESS')
        output_dirs = [pfam_dir,
                       graph_dir, query_dir, no_pfam,
                       homo_dir, pangenome_dir]
        make_working_dir(output_dirs)
    else:
        message(text=f"The program will continue what was left unfinished.", label='PROCESS')
        output_dirs = [pfam_dir, graph_dir, query_dir, no_pfam,
                       homo_dir, pangenome_dir]
        make_working_dir(output_dirs, delate=False)

    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    if not os.path.exists(os.path.join(graph_dir, 'cc_infomation.txt')):  # 如果按pfam划分这一步已经做完了
        pp = Panproteome(input_genomes_dir, pfam_dir, THREADS, parameters.pfam_db, parameters.search, evalue)  # 初始化及pfam注释
        message(text=f'Start with PFAM annotation ...', label='PROCESS')
        if parameters.r:
            message(text=f"Run deduplication.", label='PROCESS')
            pp.remove_redundant_sequences(outdir=graph_dir)
        message(text=f"After removing the redundancy, there are {len(pp)} genomes left.", label='Information')
        if parameters.pfam:
            message(text='Pfam annotate Done.', label='PROCESS')
        else:
            # 输出domain聚类的cc
            message(text='Start adding sequences ...', label='PROCESS')
            pp.add_proteome_sequence()
            message(text=f'Start building the graph ...', label='PROCESS')
            pfam_graph = PGraph(pp, no_pfam)  # 初始化，提出需要做blast的文件
            del pp  # 删除pp释放内存

            message(text='Start generating graph ...', label='PROCESS')
            pfam_graph.generate_graph()  # 构建网络
            message(text='Start puting graph file ...', label='PROCESS')
            pfam_graph.put_graph_file(graph_dir)  # 输出网络相关文件
            message(text='Start finding partition ...', label='PROCESS')
            pfam_graph.la_find_partition(graph_dir, query_dir, partition_type=partition_type,
                                         resolution_parameter=resolution_parameter)  # 社区发现
            del pfam_graph  # 删除pfam_graph
    # 输出Pfam的相关统计信息
    count_all_pfam(pfam_dir)
    if parameters.pfam:
        message(text='Pfam annotate Done.', label='PROCESS')
    elif parameters.pg:
        message(text='Pfam graph search Done.', label='PROCESS')
    else:
        message(text='Start clustering ...', label='PROCESS')
        with open(os.path.join(graph_dir, 'cc_infomation.txt'), 'r') as file:
            redundant_num = file.readline().rstrip().split(' ')[-1]
        unused_sequences_file = os.path.join(no_pfam, 'unused_sequences.faa')
        blast_dir = os.path.join(homo_dir, 'blast')
        sub_cc = os.path.join(homo_dir, 'sub_cc')
        make_working_dir([blast_dir, sub_cc], delate=False)
        if not os.path.exists(unused_sequences_file):  # 如果还没有做mapping
            faa_files = sorted(glob.glob(os.path.join(query_dir, '*.fa')))
            parallel_workflow(faa_files, blast_dir, THREADS, sub_cc,
                              parameters.blast, inflation, ssn, redundant_num)  # cc内部
            none_dict = gen_seqs_with_headers(os.path.join(no_pfam, 'none_pfam.fa'))  # none去mapping已有的cc
            message(text='Start mapping ...', label='PROCESS')
            mapping.mapping_cc_flow(os.path.join(homo_dir, 'sub_cc'), no_pfam, none_dict,
                                    THREADS, redundant_num, parameters.blast)
            message(text='Start clustering None Pfam sequences...', label='PROCESS')
        # 没有mapping的结果自己做组内cluster
            mapping.none_mmseqs_cluster(unused_sequences_file, no_pfam, THREADS,
                                        sub_cc, identity, coverage)

        # for infile in glob.glob(os.path.join(blast_dir, '*.dmnd')):  # 删除diamond的db和mmseq的tmp
        #     os.remove(infile)
        # for dir_name in os.listdir(blast_dir):
        #     if os.path.isdir(dir_name) and dir_name.startswith('tmp_'):
        #         shutil.rmtree(dir_name)

        # count_pangenome.py 统计泛基因组的信息
        object_list = cc_list(sub_cc, redundant_num, homo_dir)
        put_out_file(object_list, pangenome_dir)
    message(text='Analyse Done.', label='PROCESS')


if __name__ == '__main__':
    main()
