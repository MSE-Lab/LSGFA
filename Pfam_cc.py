#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2024-03-27 20:31
# @Author  : zhaoyu
# @File    : Pfam_cc.py

import glob
import os
import shutil
import argparse
from modules.homology_search import *


def get_parameters():
    parser = argparse.ArgumentParser(
        description='This program can proceed the all to all blast, '
                    'then get the reciprocal best hit (RBH) network, '
                    'and the cc of the subgraph is output. ')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-i', '--in', type=str, dest='input_file',
                       help='Input file')
    group.add_argument('-dir', '--in_dir', type=str, dest='input_dir', default=False,
                       help="The directory including all files. "
                            "If you use this parameter, don't use the -i parameter.")
    parser.add_argument(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f'Specify a output directory default: {os.getcwd()}')
    parser.add_argument(
        '-t', '--threads', type=int, dest='threads', default=8,
        help='Homologs searching threads. default: 8')
    parser.add_argument(
        '-id', dest='identity', type=int, default=40,
        help='The identity of homology search, default = 40.')
    parser.add_argument(
        '-c', dest='coverage', type=int, default=50,
        help='The coverage of homology search, default = 50.')
    parser.add_argument("-f", "--force", action="store_true",
                        help="Force to run by overwriting existing files")
    parser.add_argument("--no-output", action="store_true",
                        help="Do not output cc file")
    args = parser.parse_args()  # general options
    return args


def make_working_dir(out_dir, force):
    output_dirs = ['blast', 'sub_cc']
    cc_list = os.path.join(out_dir, 'cc_list.txt')
    if force:
        for dir_name in output_dirs:
            try:
                shutil.rmtree(os.path.join(out_dir, dir_name))
            except FileNotFoundError:
                pass
        try:
            os.remove(cc_list)
        except FileNotFoundError:
            pass
    # 重新创建目录
    [os.makedirs(os.path.join(out_dir, dir_), exist_ok=True) for dir_ in output_dirs]


@time_used(f"[{timing()}]Whole processing Done")
def main():
    args = get_parameters()
    input_file = args.input_file
    input_dir = args.input_dir
    out_dir = args.output_dir
    threads = str(args.threads)
    force = args.force
    identity = str(args.identity)
    coverage = str(args.coverage)
    if args.no_output:
        not_output = False
    else:
        not_output = True

    make_working_dir(out_dir, force)
    result_dir = os.path.join(out_dir, 'blast')

    if input_file:  # 输入一个文件
        agroup = DomainGroup(input_file)
        agroup.homology_search(out_dir, threads, identity, coverage)  # blast
        agroup.handle_result(os.path.join(result_dir, f'{os.path.basename(agroup.name)}.txt'))  # 处理blast结果
        agroup.build_homology_graph(out_dir, not_output)  # 建立rbh
        message(text='Analyse Done.', label='PROCESS')
    else:  # 输入一个目录
        faa_files = sorted(glob.glob(os.path.join(input_dir, '*.fa')))
        for faa in faa_files:
            agroup = DomainGroup(faa)
            agroup.homology_search(out_dir, threads, identity, coverage)  # blast
            # agroup.handle_result(os.path.join(result_dir, f'{os.path.basename(agroup.name)}.txt'))  # 处理blast结果
            # agroup.build_homology_graph(out_dir, not_output)  # 建立rbh
        message(text='Analyse Done.', label='PROCESS')

if __name__ == '__main__':
    main()
