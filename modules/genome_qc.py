#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2024-01-30 16:44
# @Author  : zhaoyu
# @File    : genome_qc.py

class Seed:
    """
    用于存放挑选出来的参考基因组
    里面的功能包含调用fastani，提取ani的结果，记录ani的值
    """
    def __init__(self, name: str = "", sequence: str = ""):
        self.name = name
        self.pairs = []
        self.ani_num = []

    def get_ani_num(self): # 返回pair的ani分数
        pass

    def get_seed(self):  # 用于选择seed
        pass

    def fastani(self):  # 调用ani
        pass

    def read_ani(self): # 处理ani结果
        pass


class ANIMatrix:
    """
    用于初始化一个表格
    行为genome的name
    列为每个genome对应的ani的分数
    """
    def __init__(self, name: str = "", sequence: str = ""):
        self.name = name

    ##### 在此处添加相关需要的函数