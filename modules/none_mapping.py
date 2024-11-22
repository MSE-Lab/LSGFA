#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Time    ：2024/11/15
# @Author  ：zhaoyu
# @File    ：none_mapping.py

import os
from modules.utils import gen_seqs_with_headers
from modules.homology_search import DomainGroup
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def choice_seqs(sub_cc_dir, none_dir):
    combined_sequences = ''
    # 遍历输入目录中的所有faa文件
    for fn in os.listdir(sub_cc_dir):
        if fn.endswith(".faa"):
            cc_name = os.path.basename(fn).split('.')[0]
            file_path = os.path.join(sub_cc_dir, fn)
            sequences_dict = gen_seqs_with_headers(file_path)
            # 选择最长的一条序列
            longest_id = max(sequences_dict, key=lambda k: len(sequences_dict[k]))
            id = f"{cc_name}|{longest_id}"
            seq = sequences_dict[longest_id]
            combined_sequences += f">{id}\n{seq}\n"
    # 将合成结果写入输出文件
    combined_sequences_file = os.path.join(none_dir, 'combined_sequences.fasta')
    with open(combined_sequences_file, 'w') as f:
        f.write(combined_sequences)
    return combined_sequences_file


def mapping_cc_flow(sub_cc_dir, none_dir, none_dict, threads):
    combined_sequences_file = choice_seqs(sub_cc_dir, none_dir)
    none_file = os.path.join(none_dir, 'none_pfam.fa')

    # none_pfam的部分与cc
    combinedGroup = DomainGroup(combined_sequences_file)
    combinedGroup.make_db(none_dir, threads)  # make db
    result_file = combinedGroup.homology_search(none_file, none_dir, threads, identity='99')
    combinedGroup.handle_result(result_file, 'sbh')  # 处理blast结果

    # 将单向最优的序列加到sub_cc内
    for sgenome, queries in combinedGroup.rbh.items():
        filename = os.path.join(sub_cc_dir, f'{sgenome}.faa')
        result = ''
        for seq_id in queries:
            result += f'>{seq_id}\n{none_dict[seq_id]}\n'
        with open(filename, 'a') as file:
            file.write(result)

    # 剩下的部分做组内blast
    rbh_ids = set()
    for queries in combinedGroup.rbh.values():
        rbh_ids.update(queries)  # 将每个查询添加到集合中
    unused_sequences = {seq_id: seq for seq_id, seq in none_dict.items() if seq_id not in rbh_ids}
    unused_sequences_file = os.path.join(none_dir, 'unused_sequences.faa')
    with open(unused_sequences_file, 'w') as output_file:
        for seq_id, sequence in unused_sequences.items():
            output_file.write(f'>{seq_id}\n{sequence}\n')
    return unused_sequences_file


def none_parallel(input_file, blast_out_dir, threads, sub_cc_dir, identity, coverage):
    result_files = []
    agroup = DomainGroup(input_file)  # cc文件
    result_file_path = os.path.join(blast_out_dir, f'{agroup.name}.txt')  # blast的结果
    if not os.path.exists(result_file_path):
        agroup.make_db(blast_out_dir, threads)  # make db
        split_files = agroup.split_file(blast_out_dir)

        with ThreadPoolExecutor(max_workers=threads) as executor:
            # 提交所有任务到线程池，并创建一个进度条
            futures = {executor.submit(agroup.homology_search, split_file, blast_out_dir,
                        threads, identity, coverage): split_file for split_file in split_files}
            # 使用 tqdm 创建进度条
            with tqdm(total=len(split_files)) as pbar:
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

        agroup.merge_files(result_file_path, result_files)
        if len(split_files) == 1:
            for file in result_files:
                os.remove(file)
        else:
            for file in split_files + result_files:
                os.remove(file)
    agroup.handle_result(result_file_path, 'rbh')  # 处理blast结果
    agroup.build_homology_graph(sub_cc_dir)  # 建立rbh