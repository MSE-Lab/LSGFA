import functools
import json
import multiprocessing as mp
import os.path
import subprocess as sp
import sys
import time
from tempfile import NamedTemporaryFile
import progressbar


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
            start = time.perf_counter() if sys.version[0] == '3' else time.clock()
            results = function(*args, **kwargs)
            end = time.perf_counter() if sys.version[0] == '3' else time.clock()
            time_use = end - start
            print(
                f'{info}: {time_use // 3600:.0f}h {(time_use % 3600) // 60:.0f}m {((time_use % 3600) % 60) % 60:.0f}s'
            )
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
