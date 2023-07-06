import subprocess
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
