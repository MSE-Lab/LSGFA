#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Time    ：2025/6/25
# @Author  ：zhaoyu
# @File    ：check_software.py

import sys
import shutil
import subprocess


def check_software_installation(search='mmseqs-search', hmm_engine='pyhmmer', blast='mmseqs-search'):
    """检查所需第三方软件是否安装(按所选方法只检查实际需要的软件)"""
    # 定义需要检查的软件列表
    required_software = [
        {'name': 'hmmscan', 'command': 'hmmscan', 'test': '-h',
         'install_hint': 'conda install -c bioconda hmmer'},
        {'name': 'mmseqs', 'command': 'mmseqs', 'test': '-h',
         'install_hint': 'conda install -c bioconda mmseqs2'},
        {'name': 'diamond', 'command': 'diamond', 'test': 'help',
         'install_hint': 'conda install -c bioconda diamond'},
        {'name': 'mcl', 'command': 'mcl', 'test': '-h',
         'install_hint': 'conda install -c bioconda mcl'},
    ]

    need_hmmscan = search == 'hmmscan' and hmm_engine == 'hmmscan'
    need_diamond = blast == 'diamond'
    required_software = [
        sw for sw in required_software
        if (sw['name'] != 'hmmscan' or need_hmmscan)
        and (sw['name'] != 'diamond' or need_diamond)
    ]

    missing_software = []

    # 检查每个软件
    for sw in required_software:
        # 1. 检查命令是否存在
        if shutil.which(sw['command']) is None:
            missing_software.append(sw)
            continue

        # 2. 检查命令是否能正常运行
        try:
            result = subprocess.run(
                [sw['command'], sw['test']],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5
            )

            # 检查输出中是否包含预期的内容，不能正常运行则提示
            if result.returncode != 0 and "not found" not in result.stderr.lower():
                missing_software.append(sw)
        except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError):
            missing_software.append(sw)

    # 如果有缺失的软件，显示错误信息并退出
    if missing_software:
        print("\n" + "=" * 80)
        print("ERROR: REQUIRED SOFTWARE NOT FOUND")
        print("=" * 80)
        print("\nThis script requires the following software to be installed and available in your PATH:\n")

        for i, sw in enumerate(missing_software, 1):
            print(f"{i}. {sw['name']} ({sw['command']})")
            print(f"   Installation hint: {sw.get('install_hint', 'see README')}")

        # 提供额外的帮助信息
        print("\nTroubleshooting tips:")
        print("1. Ensure the software is installed correctly")
        print("2. Verify that the software is in your system PATH")
        print("3. For Conda users, activate your environment before running this script")
        print("4. Try running the command manually in your terminal to verify it works")

        sys.exit(1)
