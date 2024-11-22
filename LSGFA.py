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

global OUT_DIR, THREADS


def get_parameters():
    parser = argparse.ArgumentParser(
        description='This program can get the Pfam Network. ')

    # PGraph.py部分的参数
    parser.add_argument(
        '-i', '--in', type=str, dest='input_dir', default=False,
        help='The directory including all genome files')
    parser.add_argument(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f'Specify a output directory default: {os.getcwd()}')

    parser.add_argument('-f', action='store_true',
                        help='Re-perform whole process(including pfam annotation)')
    parser.add_argument('-fg', action='store_true',
                        help='Re-perform the graph search(not including pfam annotation)')
    parser.add_argument("-fb", action="store_true",
                        help="Re-perform the homology search(blast)")

    parser.add_argument('--pfam', action='store_true',
                        help='Stop at pfam annotation)')
    parser.add_argument('--pg', action='store_true',
                        help='Stop at pfam graph search')

    parser.add_argument('-r', action='store_true',
                        help='Run deduplication')
    parser.add_argument(
        '-t', '--threads', type=int, dest='threads', default=8,
        help='Hmmscan threads. default: 8')

    # Pfam_cc.py部分的参数
    parser.add_argument(
        '-id', dest='identity', type=int, default=40,
        help='The identity of homology search, default = 40.')
    parser.add_argument(
        '-c', dest='coverage', type=int, default=50,
        help='The coverage of homology search, default = 50.')
    # parser.add_argument("--no_cc_output", action="store_true",
    #                     help="Do not output cc file")

    args = parser.parse_args()  # general options
    return args


def make_working_dir(output_dirs):
    for dir_name in output_dirs:
        try:
            shutil.rmtree(dir_name)
        except FileNotFoundError:
            pass
    # 重新创建目录
    [os.makedirs(dir_, exist_ok=True) for dir_ in output_dirs]


def work_flow(input_file, blast_out_dir, threads, sub_cc_dir, identity, coverage):
    result_files = []
    agroup = DomainGroup(input_file)  # cc文件
    result_file_path = os.path.join(blast_out_dir, f'{agroup.name}.txt')  # blast的结果
    if not os.path.exists(result_file_path):
        agroup.make_db(blast_out_dir, threads)  # make db
        split_files = agroup.split_file(blast_out_dir)
        for split_file in split_files:
            result_file = agroup.homology_search(split_file, blast_out_dir, threads, identity, coverage)
            result_files.append(result_file)
        agroup.merge_files(result_file_path, result_files)
        if len(split_files) == 1:
            for file in result_files:
                os.remove(file)
        else:
            for file in split_files + result_files:
                os.remove(file)
    agroup.handle_result(result_file_path, 'rbh')  # 处理blast结果
    agroup.build_homology_graph(sub_cc_dir)  # 建立rbh


def parallel_workflow(faa_files, blast_out_dir, threads, sub_cc_dir, identity, coverage):
    total_tasks = len(faa_files)  # 添加进度条
    make_working_dir([blast_out_dir, sub_cc_dir])
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(work_flow, faa, blast_out_dir, str(threads), sub_cc_dir,
                                   str(identity), str(coverage)): faa for faa in faa_files}
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
        cc_gene = list(gen_seqs_with_headers(member).keys())  # 该cc内的基因
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
    global OUT_DIR, THREADS, redundant_num
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters.input_dir
    OUT_DIR = parameters.output_dir
    THREADS = parameters.threads

    identity = str(parameters.identity)
    coverage = str(parameters.coverage)

    # PG.py
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

    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    pp = Panproteome(input_genomes_dir)  # 初始化
    message(text=f'Start with PFAM annotation ...', label='PROCESS')
    pp.put_pfam_file(threads=THREADS, outdir=pfam_dir)  # pfam注释
    if parameters.r:
        message(text=f"Run deduplication.", label='PROCESS')
        pp = pp.remove_redundant_sequences(outdir=OUT_DIR)
    redundant_num = len(pp)
    message(text=f"After removing the redundancy, there are {len(pp)} genomes left.", label='Information')
    if parameters.pfam:
        message(text='Pfam annotate Done.', label='PROCESS')
    else:
        if not os.path.exists(os.path.join(graph_dir, 'cc_infomation.txt')):
            # 输出domain聚类的cc
            message(text=f'Start building the graph ...', label='PROCESS')
            pfam_graph = PGraph(pp, no_pfam)  # 初始化，提出需要做blast的文件
            del pp  # 删除pp释放内存

            pfam_graph.generate_graph()  # 构建网络
            pfam_graph.put_graph_file(graph_dir)  # 输出网络相关文件
            pfam_graph.la_find_partition(graph_dir, query_dir)  # 社区发现
            del pfam_graph  # 删除pfam_graph

        if parameters.pg:
            message(text='Pfam graph search Done.', label='PROCESS')
        else:
            # Pfam_cc.py
            faa_files = sorted(glob.glob(os.path.join(query_dir, '*.fa')))
            parallel_workflow(faa_files, os.path.join(homo_dir, 'blast'), THREADS, os.path.join(homo_dir, 'sub_cc'),
                              identity, coverage)  # 并行做blast

            none_dict = gen_seqs_with_headers(os.path.join(no_pfam, 'none_pfam.fa'))  # none去mapping已有的cc
            unused_sequences_file = mapping.mapping_cc_flow(os.path.join(homo_dir, 'sub_cc'), no_pfam, none_dict,
                                                            THREADS)
            # 没有mapping的结果自己做组内blast
            mapping.none_parallel(unused_sequences_file, os.path.join(homo_dir, 'blast'), THREADS,
                                  os.path.join(homo_dir, 'sub_cc'), identity, coverage)
            for infile in glob.glob(os.path.join(homo_dir, 'blast', '*.dmnd')):  # 删除diamond的db
                os.remove(infile)
            # count_pangenome.py 统计泛基因组的信息
            object_list = cc_list(os.path.join(homo_dir, 'sub_cc'), redundant_num, homo_dir)
            put_out_file(object_list, pangenome_dir)
    message(text='Analyse Done.', label='PROCESS')


if __name__ == '__main__':
    main()
