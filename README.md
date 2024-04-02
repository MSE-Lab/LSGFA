## LSGFA简介

`LSGFA`是基于python (v3.9.16)开发的程序包，针对大规模基因组数据集进行泛基因组学研究。

本程序首先对所有蛋白进行 Pfam 结构域鉴定，再根据结构域的相似性聚类；然后在聚类序列群中构建序列相似性网络，逐个分析每个序列群的网络结构特征以判定子群是否存在；最后完成序列群的拆分。
对蛋白进行 Pfam 结构域鉴定是使用hmmscan(3.3.2)根据[Pfam数据库(36.0)](http://pfam.xfam.org/)进行鉴定，对于鉴定得到结果，保留不重叠的，长度最长的Pfam，同时，如果一条序列上的Pfam所占长度不足该序列长度的60%，则抛弃这条序列的Pfam鉴定结果。

在不同domain间进行相似性聚类时，如果两个domain有相同的Pfam，且相同的Pfam占序列长度自身序列长度的50%以上，则进行连线，同时计算两种domain间连线的权重，最后根据该网络进行社区发现，得到聚类序列群。

## 依赖包

- python-igraph (v0.10.5)
- leidenalg (v0.10.0)
- pyfasta (0.5.2)
- pandas (2.0.3)

## 软件
- blastp (v2.14.1)
- diamond (v2.1.8.162)
- hmmscan (v3.3.2)

## 数据库
- Pfam (v36.0)

## 用法

1. 基于蛋白结构域的鉴定，获得结构域相似性的聚类结果

   ```python
   python PGraph_cc.py -i input_dir -o output_dir -t threads
   ```
   如果需要重做，则添加参数 -f True
   ```python
   python PGraph_cc.py -i input_dir -o output_dir -t threads -f True
   ```
   
   根据结构域相似性聚类的结果见`query`文件夹


2. 在聚类序列群中构建RBH的网络，分析每个序列群的网络结构特征，完成序列群的拆分

   ```python
   python Pfam_cc.py -i input_file -o out_dir -t threads
   ```
   
   如果有多个需要处理的文件
   ```python
   python Pfam_cc.py -dir input_dir -o out_dir -t threads
   ```

   计算结果见指定的输出文件夹，包括：

   - `sub_cc`：子家族的文件
   - `cc_list.txt`：子家族的列表
