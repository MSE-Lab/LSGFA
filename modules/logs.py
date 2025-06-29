#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Time    ：2025/6/26
# @Author  ：zhaoyu
# @File    ：logs.py

import logging
import os
import re
import subprocess
import sys
from datetime import datetime
import platform


def setup_logger(log_dir, log_level, log_to_console=True):
    """
    设置日志记录系统
    :param log_dir: 日志目录
    :param log_level: 日志级别 (logging.DEBUG, logging.INFO, etc.)
    :param log_to_console: 是否在控制台显示日志
    :return: 配置好的logger对象
    """
    # 确保日志目录存在
    os.makedirs(log_dir, exist_ok=True)

    # 创建唯一的日志文件名（带时间戳）
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(log_dir, f"LSGFA_{timestamp}.log")

    # 创建logger
    logger = logging.getLogger("LSGFA")
    logger.setLevel(log_level)

    # 清除现有处理器（防止重复添加）
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # 创建日志格式
    log_format = logging.Formatter(
        "%(asctime)s | %(levelname)-2s | %(module)-2s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # 文件处理器（记录所有级别日志）
    file_handler = logging.FileHandler(log_filename, encoding='utf-8')
    file_handler.setFormatter(log_format)
    file_handler.setLevel(logging.DEBUG)  # 文件记录所有级别
    logger.addHandler(file_handler)

    # 控制台处理器（可选）
    if log_to_console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(log_format)
        console_handler.setLevel(log_level)  # 控制台显示指定级别及以上
        logger.addHandler(console_handler)

    # 添加启动信息
    logger.info(f"Log file: {os.path.abspath(log_filename)}")
    logger.info(f"Log level: {logging.getLevelName(log_level)}")

    return logger


def log_command_line_args(logger, args):
    """记录命令行参数"""
    logger.info("Command Line Arguments:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg:15} = {value}")


def log_software_versions(logger):
    """记录关键软件版本"""
    logger.info("Software Versions:")
    try:
        logger.info(f"  Python: {platform.python_version()}")
        # 记录关键python库版本
        try:
            import numpy
            logger.info(f"  NumPy: {numpy.__version__}")
        except ImportError:
            pass
        try:
            import pandas
            logger.info(f"  Pandas: {pandas.__version__}")
        except ImportError:
            pass
        try:
            import igraph
            logger.info(f"  igraph: {igraph.__version__}")
        except ImportError:
            pass

        # 记录外部软件版本
        software = [
            ("hmmscan", "hmmscan -h", r"HMMER (\d+\.\d+)"),
            ("mmseqs", "mmseqs -h", r"MMseqs2 Version: (\S+)"),
            ("diamond", "diamond help", r"diamond v(\d+\.\d+\.\d+\.\d+)"),
            ("mcl", "mcl --version", r"mcl (\d+-\d+)"),
            ]

        for name, cmd, pattern in software:
            try:
                result = subprocess.run(
                    cmd.split(),
                    capture_output=True,
                    text=True,
                    check=True
                )
                # 使用正则表达式提取版本号
                match = re.search(pattern, result.stdout)
                if match:
                    logger.info(f"  {name}: {match.group(1)}")
                else:
                    logger.info(f"  {name}: version unknown")
            except Exception:
                logger.info(f"  {name}: not found")

    except Exception as e:
        logger.error(f"Error logging software versions: {str(e)}")


def log_exception(logger, exc_type, exc_value, exc_traceback):
    """记录未捕获的异常"""
    logger.error("Unhandled exception occurred:", exc_info=(exc_type, exc_value, exc_traceback))
    logger.critical("Program terminated due to unhandled exception")


