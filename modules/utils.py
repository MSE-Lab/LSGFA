import functools
import json
import os.path
import time
import sys
from tempfile import NamedTemporaryFile


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
                f'[{info}]: {time_use // 3600:.0f}h {(time_use % 3600) // 60:.0f}m '
                f'{((time_use % 3600) % 60) % 60:.0f}s')
            return results

        return wrapper

    return timer


class FileOperator:
    # 用于处理输出的文件格式
    def __init__(self, name: str = "", dir_: str = "", formate: str = "text", data=None):
        self.name = name
        self.dir = dir_
        self.formate = formate
        self.data = data

    def _get_full_name(self):
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
