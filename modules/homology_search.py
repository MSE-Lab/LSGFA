import os
import gc
import subprocess
import pandas as pd
from igraph import Graph
from modules.utils import *
from collections import Counter


class DomainGroup:
	"""
	存放DomainType上具有相似性的聚类序列群，并在其中构建序列相似性网络
	"""
	def __init__(self, fasta_name, graph: Graph = None):
		super().__init__()
		self.content = gen_seqs_with_headers(fasta_name)
		self.file = fasta_name
		self.name = os.path.basename(fasta_name).split('.')[0]
		self.db = None
		self.rbh = None
		self.graph = graph

	def __len__(self):
		return len(self.content)

	def homology_search(self, input_file, blast_dir, threads, method, num, identity=40, cover=50):
		blast_cmd = ''
		res = ''
		if method == 'diamond':
			# 进行ata blast
			res = os.path.join(blast_dir, f'result_{os.path.basename(input_file)}')  # 比较结果
			blast_cmd = ' '.join([
				'diamond', 'blastp', '--more-sensitive', '-p', str(threads), '-q', input_file, '-d', self.db,
				'--evalue 1e-5 -f 6', '--out', res, '--query-cover', str(cover), '--subject-cover', str(cover),
				'-k', '0', '--id', str(identity)])
		elif method == 'mmseqs-search':
			res = os.path.join(blast_dir, f'{self.name}.txt')
			blast_cmd = ' '.join([
				'mmseqs', 'easy-search', input_file, self.file, res, os.path.join(blast_dir, f'tmp_{self.name}'),
				'-s', '7.5', '-e', '1.000E-05', '--threads', str(threads), '--max-seqs', str(num),
				'--min-seq-id', str(int(identity)/100), '-c', str(int(cover)/100)])
		blast_cap = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		blast_cap.communicate()
		return res

	def make_db(self, out_dir, threads):
		# 进行ata blast
		name = self.name  # 要进行比较的fa文件的名字
		os.makedirs(out_dir, exist_ok=True)
		db = os.path.join(out_dir, name)  # db的位置
		self.db = db
		# 建db
		# message(text='Make database...', label='PROCESS')
		db_cmd = ' '.join(['diamond', 'makedb', '--in', self.file, '--db', db, '--threads', str(threads)])
		db_cap = subprocess.Popen(db_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		db_cap.communicate()

	def split_file(self, out_dir):
		result_files = []
		gene_num = len(self.content.keys())
		if gene_num <= 2000:  # 当序列数小于4000时不拆分
			return [self.file]
		else:
			# message(text='Split file...', label='PROCESS')
			split_dir = os.path.join(out_dir, 'split_file')
			os.makedirs(split_dir, exist_ok=True)
			for i in range(0, gene_num, 500):
				chunk = list(self.content.items())[i:i + 500]
				output_file = f'{self.name}_part{i//500+1}.fa'
				text_content = ''
				for key, value in chunk:
					text_content += f'>{key}\n{value}\n'
				with open(os.path.join(split_dir, output_file), 'w') as f:
					f.write(text_content)
				result_files.append(os.path.join(split_dir, output_file))
			# message(text=f'There are {len(result_files)} split files...', label='Information:')
			return result_files

	# 合并结果文件
	@staticmethod
	def merge_files(output_file, result_files):
		# message(text='Merge split files...', label='PROCESS')
		if len(result_files) == 1:
			os.rename(result_files[0], output_file)
		else:
			with open(output_file, 'w') as outfile:
				for result_file in result_files:
					with open(result_file, 'r') as infile:
						outfile.write(infile.read())
		return output_file

	def handle_result(self, result_file, tag: str):
		"""
		处理blast的结果文件
		cc自己的blast结果取双向最优
		none到cc的balst取单向最优，即none为query，cc_combine为subject
		"""
		# 读取文件内
		rbh = None
		# 仅读取以下列名的内容
		necessary_columns = ['query', 'subject', 'id', 'evalue', 'bitscore']
		data = pd.read_csv(result_file, sep='\t', header=None,
						names=['query', 'subject', 'id', 'length', 'mismatch',
								'gapopen', 'qstart', 'qend', 'sstart',
								'send', 'evalue', 'bitscore'],
						usecols=necessary_columns)
		# 提取基因组
		data[['qgenome', 'sgenome']] = data[['query', 'subject']].map(lambda x: x.split('|')[0])
		# filtered_data = data[data['qgenome'] != data['sgenome']]  # 过滤掉同组的内容
		filtered_data = data[data['query'] != data['subject']]
		del data

		if tag == 'rbh':  # cc内部的
			# result_list = filtered_data.sort_values(by=['id', 'bitscore', 'evalue'],
			# 										ascending=[False, False, True]).groupby(
			# 	['query', 'sgenome']).head(1)[['query', 'subject']]
			# result_list['pair'] = result_list.apply(lambda row: tuple(sorted([row['query'], row['subject']])), axis=1)
			# rbh = [key for key, value in Counter(result_list['pair']).items() if value == 2]
			df_sorted = filtered_data.sort_values(by=['id', 'bitscore', 'evalue'], ascending=[False, False, True])
			max_values = df_sorted.groupby(['query', 'sgenome']).agg({
				'id': 'max',
				'bitscore': 'max',
				'evalue': 'min'
			}).reset_index()
			result_list = pd.merge(df_sorted, max_values, on=['query', 'sgenome', 'id', 'bitscore', 'evalue'])
			result_list['pair'] = result_list.apply(lambda row: tuple(sorted([row['query'], row['subject']])), axis=1)
			rbh = [key for key, value in Counter(result_list['pair']).items() if value == 2]

		elif tag == 'sbh':  # none去blast的
			result_list = filtered_data.sort_values(by=['id', 'bitscore', 'evalue'],
													ascending=[False, False, True]).groupby(
				['query']).head(1)[['query', 'sgenome']]
			rbh = result_list.groupby('sgenome')['query'].apply(list).to_dict()
		del filtered_data, result_list
		gc.collect()
		# 存放双向最优的配对结果
		self.rbh = rbh

	def build_homology_graph(self, out_dir):  # 构建rbh网络
		"""
		构建序列间的双向最优匹配网络
		:param out_dir:输出目录
		:return:子图列表
		"""
		# message(text='Build homology graph...', label='PROCESS')
		vs_list = list(self.content.keys())
		cc_graph = Graph()
		cc_graph.add_vertices(vs_list)  # 添加点
		cc_graph.add_edges(self.rbh)  # 添加边
		# cc_graph.write_gml(os.path.join(out_dir, f'{self.name}.gml'))

		components_list = []  # 存放是sog的子图
		components = cc_graph.components()  # 子图
		for component in components:
			subgraph = cc_graph.subgraph(component)
			components_list.append(list(subgraph.vs['name']))
		# message(text=f'{len(components_list)} sub-Pfam were found.', label='Information')
		self.put_file(components_list, out_dir)  # 输出文件
		return

	def put_file(self, components_list, out_dir):
		cc_num = 1
		os.makedirs(out_dir, exist_ok=True)
		for cc in components_list:
			result = ''
			for seq_id in cc:
				result += f'>{seq_id}\n{self.content[seq_id]}\n'
			cc_name = FileOperator(f'{self.name}_{cc_num}.faa', out_dir, data=result)
			cc_name.write()
			cc_num += 1
