#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Time    ：2024/10/30
# @Author  ：zhaoyu
# @File    ：LSGFA.py

from modules.panproteome import *
from modules.build_graph import *
import shutil
import argparse
from modules.homology_search import *
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


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
    parser.add_argument("--no_cc_output", action="store_true",
                        help="Do not output cc file")

    args = parser.parse_args()  # general options
    return args


def make_working_dir(output_dirs):
    for dir_name in output_dirs:
        try:
            shutil.rmtree(os.path.join(OUT_DIR, dir_name))
        except FileNotFoundError:
            pass
    try:
        os.remove(os.path.join(OUT_DIR, 'orthogroups.csv'))
    except FileNotFoundError:
        pass
    # 重新创建目录
    [os.makedirs(os.path.join(OUT_DIR, dir_), exist_ok=True) for dir_ in output_dirs]


def work_flow(input_file, out_dir, threads, result_dir, identity, coverage, not_output):
    result_files = []
    agroup = DomainGroup(input_file)
    result_file_path = os.path.join(result_dir, f'{agroup.name}.txt')
    if not os.path.exists(result_file_path):
        agroup.make_db(out_dir, threads)  # make db
        split_files = agroup.split_file(out_dir)
        for split_file in split_files:
            result_file = agroup.homology_search(split_file, result_dir, threads, identity, coverage)
            result_files.append(result_file)
        agroup.merge_files(result_file_path, result_files)
        if len(split_files) == 1:
            for file in result_files:
                os.remove(file)
        else:
            for file in split_files + result_files:
                os.remove(file)
    agroup.handle_result(result_file_path)  # 处理blast结果
    agroup.build_homology_graph(out_dir, not_output)  # 建立rbh


def parallel_workflow(faa_files, out_dir, threads, identity, coverage, not_output):
    total_tasks = len(faa_files)  # 添加进度条
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(work_flow, faa, out_dir, str(threads), os.path.join(out_dir, 'blast'),
                                   str(identity), str(coverage), not_output): faa for faa in faa_files}
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


def cc_list(cc_dir, species_num):
    # 获取路径中的faa文件，处理为cc:gene_num的字典
    members = glob.glob(os.path.join(cc_dir, '*.faa'))
    object_list = []
    for member in members:
        cc_file = Fasta(member)
        cc_name = os.path.basename(member).split('.')[0]
        cc_species = set([i.split('|')[0] for i in cc_file.keys()])
        marker_ = define_marker(len(cc_species)/int(species_num))
        aGeneGroup = GeneGroup(name=cc_name, marker=marker_, species=cc_species)
        object_list.append(aGeneGroup)
        os.remove(f"{member}.flat")
        os.remove(f"{member}.gdx")
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
    global OUT_DIR, THREADS
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters.input_dir
    OUT_DIR = parameters.output_dir
    THREADS = parameters.threads

    identity = str(parameters.identity)
    coverage = str(parameters.coverage)

    # PG.py
    pfam_dir = os.path.join(OUT_DIR, 'pfam')
    graph_dir = os.path.join(OUT_DIR, 'graph')
    query_dir = os.path.join(OUT_DIR, 'query')
    no_pfam = os.path.join(OUT_DIR, 'none_pfam')

    # 判断是否进行重做
    if parameters.fg:  # 从PGraph重做
        message(text=f"Re-perform the graph search(not including annotation).", label='PROCESS')
        output_dirs = ['graph', 'query', 'none_pfam',
                       'none_pfam/fa_file',
                       'blast', 'sub_cc', 'split_file']
        make_working_dir(output_dirs)
    if parameters.fb:
        message(text=f"Re-perform the homology search(blast).", label='PROCESS')
        output_dirs = ['blast', 'sub_cc', 'split_file']
        make_working_dir(output_dirs)
    if parameters.f:
        message(text=f"Re-perform whole process(including pfam annotation).", label='PROCESS')
        output_dirs = ['graph', 'query', 'none_pfam', 'pfam',
                       'none_pfam/fa_file',
                       'blast', 'sub_cc', 'split_file']
        make_working_dir(output_dirs)

    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    pp = Panproteome(input_genomes_dir)  # 初始化
    message(text=f'Start with PFAM annotation ...', label='PROCESS')
    pp.put_pfam_file(threads=THREADS, outdir=pfam_dir)  # pfam注释
    if parameters.r:
        message(text=f"Run deduplication.", label='PROCESS')
        pp = pp.remove_redundant_sequences(outdir=OUT_DIR)
    message(text=f"After removing the redundancy, there are {len(pp)} genomes left.", label='Information')
    if parameters.pfam:
        message(text='Pfam annotate Done.', label='PROCESS')
    else:
        if not os.path.exists(os.path.join(graph_dir, 'cc_infomation.txt')):
            # 输出domain聚类的cc
            message(text=f'Start building the graph ...', label='PROCESS')
            pfam_graph = PGraph(pp, no_pfam)  # 初始化，提出需要做blast的文件
            pfam_graph.generate_graph()  # 构建网络
            pfam_graph.put_graph_file(graph_dir)  # 输出网络相关文件
            pfam_graph.la_find_partition(graph_dir, query_dir)  # 社区发现
            # pfam_graph.put_out_cc(partitions, query_dir)

        if parameters.pg:
            message(text='Pfam graph search Done.', label='PROCESS')
        else:
            # Pfam_cc.py
            faa_files = sorted(glob.glob(os.path.join(OUT_DIR, 'query', '*.fa')))
            if parameters.no_cc_output:  # 判断是否输出cc的文件
                not_output = False
            else:
                not_output = True
            parallel_workflow(faa_files, OUT_DIR, THREADS, identity, coverage, not_output)  # 并行做blast
            for infile in glob.glob(os.path.join(OUT_DIR, 'blast', '*.dmnd')):  # 删除diamond的db
                os.remove(infile)
            if parameters.no_cc_output:  # 判断是否输出cc的文件
                message(text="No sub_cc file, can't count pangenome infomation.", label='PROCESS')
            else:
                # count_pangenome.py 统计泛基因组的信息
                object_list = cc_list(os.path.join(OUT_DIR, 'sub_cc'), len(pp))
                put_out_file(object_list, OUT_DIR)
    message(text='Analyse Done.', label='PROCESS')


if __name__ == '__main__':
    main()
