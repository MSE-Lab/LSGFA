import functools
import json
import multiprocessing as mp
import os.path
import subprocess as sp
import sys
import time
from tempfile import NamedTemporaryFile
import progressbar
import glob
from collections import Counter
import os
import numpy as np
import argparse


def timing():
    """
    Get current time
    return: formatted time
    """
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def message(text, label='', depart=False, sleep_time=0):
    """
    Print information on screen
    param out: exit program or not
    param text: the information will be sent.
    param label: ERROR, WARNING, '' (for prompt)
    return:
    """
    if label:
        print(f"[{timing()}] {label}: {text}")
    else:
        print(f"[{timing()}] {text}")
    if sleep_time:
        time.sleep(sleep_time)
    if depart:
        sys.exit(0)


def make_temp_file(prefix, close=False):
    temp = NamedTemporaryFile(mode='w', prefix=prefix, delete=False)
    if close:
        temp.close()
    else:
        pass
    return temp


def time_used(info=''):
    def timer(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            start = time.perf_counter()  # 使用 perf_counter 获取高精度计时
            results = function(*args, **kwargs)
            end = time.perf_counter()
            time_use = end - start
            # 输出耗时，以秒为单位，保留小数点后两位
            print(f'{info}: {time_use:.2f}s')
            # print(
            #     f'{info}: {time_use // 3600:.0f}h {(time_use % 3600) // 60:.0f}m {((time_use % 3600) % 60) % 60:.0f}s'
            # )
            return results
        return wrapper
    return timer


class FileOperator:
    def __init__(self, name: str = "", dir_: str = "", formate: str = "text", data=None):
        self.name = name
        self.dir = dir_
        self.formate = formate
        self.data = data

    def _get_full_name(self):
        if "/" in self.name:
            return self.name
        else:
            return os.path.join(self.dir, self.name)

    def read(self):
        with open(self._get_full_name()) as f:
            if self.formate == "json":
                self.data = json.load(f)
            else:
                self.data = f.readlines()

    def write(self):
        with open(self._get_full_name(), "w") as f:
            if self.formate == "json":
                f.write(json.dumps(self.data))
            else:
                f.writelines(self.data)

    def remove(self):
        os.remove(self._get_full_name())


class CmdManger:
    # 用于管理命令
    def __init__(self, process: str = "", cmd: str = "", thread: str = "1"):
        self.process = process
        self.cmd = cmd
        self.thread = str(thread)

    def homology_searching(self, query, db, out_name, id):
        if self.process == 'blastp':
            self.cmd = ' '.join(['blastp', '-query', query, '-db', db, "-outfmt 6 -evalue 1e-5", "-out", out_name])
        elif self.process == 'diamond':
            self.cmd = ' '.join([
                'diamond', 'blastp', '--more-sensitive', '-p', self.thread, '-q', query, '-d', '%s.dmnd' % db,
                '--evalue 1e-5 -f 6', '--out', out_name, '--quiet', '--query-cover', '50', '--subject-cover', '50',
                '-k', '0', '--id', id])
        elif self.process == 'mmseqs':
            self.cmd = ' '.join([
                'mmseqs', 'easy-search', query, db, out_name, '/temp', '--threads', self.thread, '-v', '1',
                '--format-mode', '0', '--remove-tmp-files', '-s', '7.5', '-e', '1e-5', ])

    def make_db(self, input_name, db):
        if self.process == 'blastp':
            self.cmd = ' '.join(['makeblastdb', '-dbtype', 'prot', '-in', input_name, '-out', db])
        elif self.process == 'diamond':
            self.cmd = ' '.join(['diamond', 'makedb', '--in', input_name, '--db', db, '--threads', self.thread])

    def mcl(self, abc_file, inflation, out):
        self.cmd = ' '.join(
            ['/media/disk4/conda_envs/UPhO/bin/mcl', abc_file, '--abc', '-I', inflation, '-o', out, '-te', self.thread,
             '-V -all'])

    def clustalo_aln(self, fa_file, out_name):
        self.cmd = ' '.join(['clustalo', '-i', fa_file, '-o', out_name])

    def fasttree(self, aln_file, out_name):
        self.cmd = ' '.join(['/opt/miniconda3/bin/fasttree', aln_file, '>', out_name])


class CallCmd:

    def __init__(self, cmd_list: list = None, process_info: str = "", threads: int = 8,
                 parallel: [False, True] = False):
        self.process_info = process_info
        self.threads = threads
        self.parallel = parallel
        self.cmd_list = cmd_list

    def call_cmd(self, cmd, queue: mp.Manager().Queue() = None):
        pro = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        pro.wait()
        if self.parallel:
            queue.put(pro.returncode)
        else:
            return pro.returncode

    def manager_queue(self, queue: mp.Manager().Queue() = None):
        task_num = len(self.cmd_list)
        task_stat_list = []
        messages = f'[{timing()}]{self.process_info:<20}|'
        progressbar_widgets_set = [messages, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
        bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=task_num)
        bar.start()
        done_num = 0
        while True:
            cmd_stat = queue.get()
            task_stat_list.append(cmd_stat)
            done_num += 1
            bar.update(done_num)
            if done_num >= task_num:
                break
        bar.finish()
        return task_stat_list

    def processing(self):
        print(f'[{timing()}]{self.process_info:.<20}')
        for cmd in self.cmd_list:
            self.call_cmd(cmd=cmd)

    def parallel_process(self):
        queue = mp.Manager().Queue()
        pool = mp.Pool(self.threads)
        for cmd in self.cmd_list:
            pool.apply_async(func=self.call_cmd, args=(cmd, queue))
        self.manager_queue(queue)
        pool.close()
        pool.join()


def gen_seqs_with_headers(fn, extract_ids=False):
    with open(fn) as f:
        fh = f.readlines()

    if not fh[0].startswith('>'):
        message(text=f"'{fn}' is not a fasta file", label='Error')
        exit(1)

    header = None
    seqs = []
    gene_seqs = dict()
    gene_ids = []
    # extract_ids这个参数的意思是，是否仅提取ID
    # False表示不仅仅提取ID，序列也要提取
    # True表示仅仅提取ID，不提取序列
    if not extract_ids:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:  # 提取整个文件
                    gene_seqs[header] = "".join(seqs)
                header = line[1:]
                seqs = []
            else:
                seqs.append(line)
        # 保存最后一个读取到的序列
        if seqs and header is not None and not extract_ids:
            gene_seqs[header] = "".join(seqs)
        return gene_seqs
    if extract_ids:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                gene_ids.append(header)
        return gene_ids


def count_da(data: dict):
    none_domain = 0
    da_list = list()
    da_pfam_num = list()
    lencov_list = list()
    for p_dict in data.values():
        da = p_dict['Domain']
        if len(da) == 0:
            none_domain += 1  # 记录有多少个没有注释到pfam的
        else:
            da_string = ','.join(sorted(da))
            da_list.append(da_string)  # 所有da的list（DA组合）
            da_pfam_num.append(len(da))  # da有多少种pfam组成
            lencov = str(sum(p_dict['LenCov']))  # da的覆盖度
            lencov_list.append(lencov)
    domain_num_dict = Counter(da_list)  # 统计每种da的数量
    da_pfam_num = Counter(da_pfam_num)
    return none_domain, domain_num_dict, lencov_list, da_pfam_num


def count_all_pfam(file_path):
    pfam_files = glob.glob(os.path.join(file_path, '*.pfam'))
    pfam_dict = dict()

    for file in pfam_files:
        with open(file, 'r') as f:
            data = json.load(f)
        pfam_dict.update(data)

    # domain的计算
    none_domain, domain_num_dict, lencov_list, da_pfam_num = count_da(pfam_dict)
    len_domain = len(domain_num_dict.keys())
    # 计算有多少序列没有pfam，每种da的数量，da的覆盖度
    len_bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    bin_counts = np.histogram(lencov_list, bins=len_bins)[0]
    total_count = len(lencov_list)

    genes_num = len(pfam_dict.keys())
    # 输出结果
    print("=" * 50)
    print(f"{'Summary of Genome Data':^50}")
    print("=" * 50)
    print(f'{"Genomes num:":<30} {len(pfam_files):>15}')  # 基因组数量
    print(f'{"Seqs num:":<30} {genes_num:>15}')  # 基因数量
    print(f'{"None pfam seqs num:":<30} {none_domain:>15}')  # 注释到None的序列数
    print(f'{"Percentage of None Pfam sequences:":<30} {none_domain / genes_num * 100:>10.2f}%')  # None的百分比
    print(f'{"DA type num:":<30} {len_domain:>15}')  # da的种类数

    print("\n" + "=" * 50)
    print(f"{'DA Composition':^50}")
    print("=" * 50)
    for key, value in da_pfam_num.items():  # DA的组成情况
        print(f"{f'DA consisting of {key} PFAMs:':<40} {value:>10}")

    print("\n" + "=" * 50)
    print(f"{'Coverage Percentage by Range':^50}")
    print("=" * 50)
    for i in range(len(bin_counts)):
        percentage = (bin_counts[i] / total_count) * 100
        print(f'{len_bins[i]} ~ {len_bins[i + 1]:<15} {percentage:.2f}%')

    print("=" * 50)


def validate_range(value, min_val, max_val, name):
    """
    验证值是否在指定范围内
    :param value: 要验证的值
    :param min_val: 最小值
    :param max_val: 最大值
    :param name: 参数名称（用于错误信息）
    :return: 如果验证通过，返回原始值
    :raises: argparse.ArgumentTypeError 如果值无效
    """
    try:
        fvalue = float(value)
    except ValueError:
        message(text="'{value}' is not a valid number for {name}", label='Error')
        exit(1)
    else:
        if not (min_val <= fvalue <= max_val):
            message(text=f"{name} must be between {min_val} and {max_val} (got {fvalue})",
                    label='Error')
            exit(1)

def validate_fasta_ids(seq_id, genome_id):
    parts = seq_id.split('|')
    if '|' not in seq_id or len(parts) != 2 or parts[0] != genome_id:
        message(text=f"'{seq_id}' is not a valid id", label='Error')
        print("Examples of valid IDs:")
        print("  >Genome1|Gene001")
        print("Please correct your FASTA files and try again")
        sys.exit(1)
