## LSGFA简介

`LSGFA`是基于python (v3.9.16)开发的程序包，针对大规模基因组数据集进行泛基因组学研究。

本程序首先对所有蛋白进行 Pfam 结构域鉴定，再根据结构域的相似性聚类；然后在聚类序列群中构建序列相似性网络，逐个分析序列群的网络结构特征以判定子群是否存在；最后完成序列群的拆分。

对蛋白进行 Pfam 结构域鉴定是使用hmmscan(3.3.2)根据[Pfam数据库(36.0)](http://pfam.xfam.org/)进行鉴定。对于鉴定得到结果，保留不重叠的，长度最长的Pfam，同时，如果一条序列上的Pfam所占长度不足该序列长度的60%，则抛弃这条序列的Pfam鉴定结果。

在不同domain间进行相似性聚类时，如果两个domain有相同的Pfam，且相同的Pfam占序列长度自身序列长度的50%以上，则进行连线，同时计算两种domain间连线的权重，最后根据该网络结构获取子图，得到聚类序列群。

## 依赖包

- python-igraph (v0.10.5)
- leidenalg (v0.10.0)
- pyfasta (v0.5.2)
- pandas (v2.0.3)

## 软件
- blastp (v2.14.1)
- diamond (v2.1.8.162)
- hmmscan (v3.3.2)

## 数据库
- Pfam (v36.0)

  在运行本程序前，请将 Pfam 数据库下载modules的database目录，解压后并进行本地化

​		1. 进入[数据库页面](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/)查看要下载的数据库

​		2. 将数据库下载到database目录内

​		3. 对数据库进行本地化

   ``` shell
   # 下载Pfam36.0
   wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz
   # 解压数据库文件
   gzip -d Pfam-A.hmm.gz
   # 将数据库本地化
   hmmpress Pfam-A.hmm
   ```


## 用法

1. 基于蛋白结构域的鉴定，获得结构域相似性的聚类结果

   **生成帮助页面**    快速查看软件可用的命令行选项

   ```python
   python PGraph_cc.py -h
   ```
   **一般使用方法**

   ```python
   python PGraph_cc.py -i input_dir -o output_dir -t threads
   ```
   a. 根据结构域的相似性建立的网络见 `graph` 文件夹
   - `cc_infomation.txt`：聚类结果的相关信息
   - `node_genes.txt`：网络节点属性的文件
   - `pfam_graph.gml`、`pfam_graph.txt`：网络结构的文件

   b. 没有注释到 Pfam 的序列处理结果见 `none_pfam` 文件夹

   - `none_pfam.fa`:没有注释到 Pfam 的序列

   c. 基因组的 Pfam 注释结果见 `pfam` 文件夹
   d. 根据结构域相似性聚类的结果见 `query` 文件夹

   **测试示例**

   本程序内置了两个测试数据Staphylococcus和Streptomyces

   可选择其一对软件进行测试，保证软件可以顺利运行

   ```python
   python PGraph_cc.py -i /testdata/Streptomyces -o output_dir -t threads
   ```

   

2. 在聚类序列群中构建双向最优的网络，分析每个序列群的网络结构特征，完成序列群的拆分

   这一部分用于处理结构域相似性聚类的结果，输入文件为 .fa 结尾的文件

   如果使用的是PGraph_cc.py的结果，输入目录应该是query目录

   **生成帮助页面**    快速查看软件可用的命令行选项

   ```python
   python Pfam_cc.py -h
   ```
   **一般使用方法**

   ```python
   python Pfam_cc.py -i input_file -o out_dir -t threads
   ```

   如果有多个需要处理的文件
   ```python
   python Pfam_cc.py -dir input_dir -o out_dir -t threads
   ```

   计算结果见指定的输出文件夹，包括：

   ​	a. `sub_cc`：子家族的文件

   ​	b. `blast`:聚类序列群进行blast的结果

   ​	c. `cc_list.txt`：子家族的列表
