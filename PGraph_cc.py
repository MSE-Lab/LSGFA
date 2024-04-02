#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2023-07-11 17:27
# @Author  : zhaoyu
# @File    : PGraph_cc.py
from modules.panproteome import *
import shutil
import argparse

global OUT_DIR, THREADS

def get_parameters():
    parser = argparse.ArgumentParser(
        description='This program can get the Pfam Network. ')
    parser.add_argument(
        '-i', '--in', type=str, dest='input_dir', default=False,
        help='The directory including all genome files')
    parser.add_argument(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f'specify a output directory default: {os.getcwd()}')
    parser.add_argument(
        '-f', '--force', type=str, dest='re_run', default=False,
        help=f'Re-perform the homology search; default: False')
    parser.add_argument(
        '-t', '--threads', type=int, dest='threads', default=8,
        help='Hmmscan threads. default: 8')
    args = parser.parse_args()  # general options
    return args


def make_working_dir(re_run):
    output_dirs = ['graph', 'query', 'none_pfam',
                   'none_pfam/fa_file']
    if re_run:  # 如果强制重做
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


@time_used(f"[{timing()}]Whole processing Done")
def main():
    global OUT_DIR, THREADS
    # general options
    parameters = get_parameters()
    input_genomes_dir = parameters.input_dir
    OUT_DIR = parameters.output_dir
    THREADS = parameters.threads
    re_run = parameters.re_run

    pfam_dir = os.path.join(OUT_DIR, 'pfam')
    graph_dir = os.path.join(OUT_DIR, 'graph')
    query_dir = os.path.join(OUT_DIR, 'query')
    no_pfam = os.path.join(OUT_DIR, 'none_pfam')

    make_working_dir(re_run)  # 创建需要的目录
    max_genome = len([file for file in os.listdir(input_genomes_dir) if file.split(".")[-1] == 'faa'])
    message(text=f'genomes Numbers: {max_genome}', label='Information')

    pp = Panproteome(input_genomes_dir)  # 初始化
    message(text=f'Start with PFAM annotation ...', label='PROCESS')
    pp.put_pfam_file(threads=100, outdir=pfam_dir)  # pfam注释

    # 输出domain聚类的cc
    message(text=f'Start building the graph ...', label='PROCESS')
    pfam_graph = PGraph(pp, no_pfam)  # 初始化，提出需要做blast的文件
    pfam_graph.generate_graph()  # 构建网络
    pfam_graph.put_graph_file(graph_dir)  # 输出网络相关文件
    partitions = pfam_graph.la_find_partition(graph_dir)  # 社区发现
    pfam_graph.put_out_cc(partitions, query_dir)


if __name__ == '__main__':
    main()
