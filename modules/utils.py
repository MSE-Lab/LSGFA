import subprocess
import time
import sys


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


def run_command(cmd):
    """
    Test command, if stderr info get out, test result is negative
    :param cmd: command
    :return:
    """
    cap = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sdo = [x for x in cap.stdout]
    sde = [x for x in cap.stderr]
    if len(sdo) > 0 and len(sde) == 0:
        message("test run '%s' - successful" % cmd)
        return True
    else:
        message("test run '%s' - failed" % cmd)
        return False

