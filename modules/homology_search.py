import os
import subprocess
import pandas as pd
from igraph import Graph
from modules.utils import *
from pyfasta import Fasta
from collections import Counter


class DomainGroup:
	# 用于存放每个cc的内容
	# Ngroup是一个dict，里面存放的是gene_id:seq
	def __init__(self, fasta_name, graph: Graph = None):
		super().__init__()
		self.content = dict()
		self.file = fasta_name
		self.name = os.path.basename(fasta_name).split('.')[0]
		self.rbh = None
		self.graph = graph
		for gene_id, seq in Fasta(fasta_name).items():
			self.content[gene_id] = seq
		os.remove(f"{fasta_name}.flat")
		os.remove(f"{fasta_name}.gdx")

	def __len__(self):
		return len(self.content)

	def homology_search(self, out_dir, threads, id='40', cover='50'):
		# 进行ata blast
		name = self.name  # 要进行比较的fa文件的名字
		blast_dir = os.path.join(out_dir, 'blast')
		os.makedirs(blast_dir, exist_ok=True)
		db = os.path.join(blast_dir, name)  # db的位置
		res = os.path.join(blast_dir, f'{name}.txt')  # 比较结果
		# 建db
		message(text='Make database...', label='PROCESS')
		db_cmd = ' '.join(['diamond', 'makedb', '--in', self.file, '--db', db, '--threads', threads])
		db_cap = subprocess.Popen(db_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		db_cap.communicate()
		# blast
		message(text='Blast...', label='PROCESS')
		blast_cmd = ' '.join([
			'diamond', 'blastp', '--more-sensitive', '-p', threads, '-q', self.file, '-d', '%s.dmnd' % db,
			'--evalue 1e-5 -f 6', '--out', res, '--quiet', '--query-cover', cover, '--subject-cover', cover,
			'-k', '0', '--id', id])
		print(blast_cmd)
		blast_cap = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		blast_cap.communicate()

	def handle_result(self, result_file):
		# 读取文件内
		data = pd.read_csv(result_file, sep='\t', header=None,
						   names=['query', 'subject', 'id', 'length', 'mismatch', 'gapopen',
								  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
		# 提取基因组
		data[['qgenome', 'sgenome']] = data[['query', 'subject']].map(lambda x: x.split('|')[0])
		# filtered_data = data[data['qgenome'] != data['sgenome']]  # 过滤掉同组的内容
		filtered_data = data[data['query'] != data['subject']]
		# 排序后分组，然后去top1
		result_list = filtered_data.sort_values(by=['id', 'evalue', 'bitscore'],
												ascending=[False, True, True]).groupby(
												['query', 'sgenome']).head(1)[['query', 'subject']]

		result_list['pair'] = result_list.apply(lambda row: tuple(sorted([row['query'], row['subject']])), axis=1)
		rbh_list = [key for key, value in Counter(result_list['pair']).items() if value == 2]
		# 存放双向最优的配对结果
		self.rbh = rbh_list
		return rbh_list

	def build_homology_graph(self, out_dir, cc_file):  # 构建rbh网络
		vs_list = list(self.content.keys())
		cc_graph = Graph()
		cc_graph.add_vertices(vs_list)  # 添加点
		cc_graph.add_edges(self.rbh)  # 添加边
		cc_graph.write_gml(os.path.join(out_dir, f'{self.name}.gml'))

		components_list = []  # 存放是sog的子图
		components = cc_graph.components()  # 子图
		for component in components:
			subgraph = cc_graph.subgraph(component)
			components_list.append(list(subgraph.vs['name']))
		message(text=f'{len(components_list)} sub-Pfam were found.', label='Information')
		self.put_file(components_list, out_dir, cc_file)  # 输出文件
		return components_list

	def put_file(self, components_list, out_dir, cc_file=True):
		if cc_file:  # 输出seq
			cc_num = 1
			subcc_dir = os.path.join(out_dir, 'sub_cc')
			os.makedirs(subcc_dir, exist_ok=True)
			for cc in components_list:
				result = ''
				for seq_id in cc:
					result += f'>{seq_id}\n{self.content[seq_id]}\n'
				cc_name = FileOperator(f'{self.name}_{cc_num}.faa', subcc_dir, data=result)
				cc_name.write()
				cc_num += 1
		# 只输出gene_id不输出seq
		cc_ = 1
		cc_result = ''
		for cc in components_list:
			protein_lists = ','.join(cc)
			cc_result += f'{self.name}_{cc_}\n{protein_lists}\n'
			cc_ += 1
		with open(os.path.join(out_dir, 'cc_list.txt'), 'a') as f:
			f.write(cc_result)
