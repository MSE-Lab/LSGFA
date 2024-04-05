#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2024-01-08 15:36
# @Author  : zhaoyu
# @File    : upho.py
import glob
import re
import os
from math import fsum


# GLOBAL VARIABLES. MODIFY IF NEEDED
sep = '|'
gsep = re.escape(sep)


# CLASS DEFINITIONS
class Split:
	def __init__(self):
		self.vecs = None
		self.branch_length = None
		self.support = None
		self.name = None
		self.label = []


class Myphylo:
	"""A class for newick trees"""

	def __init__(self, newick):
		self.leaves = get_leaves(newick)
		self.splits = []
		self.ortho = []
		self.costs = {}  # Dictionary of leaf cost for inparalog evaluation'
		self.newick = newick.strip('\n')
		for leaf in self.leaves:  # Cost initialized
			self.costs[leaf] = 1.0
		split_decomposition(self)

# FUNCTION DEFINITIONS


def get_leaves(string):
	"""Find leaves names in newick files using regexp.
	Leaves names are composed of alpha numeric characters,
	underscore and a special field delimiter"""
	pattern = "(?<=[,\(])\w.*?(?=[,:;\)])"
	#    pattern = "[^\(\),;:\[\]]+%s[^\(\),;:\[\]]+" % gsep
	leaves = re.findall(pattern, string)
	return leaves


def spp_in_list(alist):
	"""Return the species from a list of sequence identifiers"""
	spp = []
	for i in alist:
		spp.append(i.Split(sep)[0])
	return spp


def complement(sub, whole):
	"""Return elements in Whole that are not present in Sub"""
	complement = set(whole) - set(sub)
	return list(complement)


def split_decomposition(tree):
	"""Add a list of splits class objects to myPhylo Class object"""
	# Part I. Where we identify matching parenthesis in the newick.
	newick = tree.newick
	leaves = tree.leaves
	parenthesis = {}  # Empty dictionary where to store parenthesis identifiers and a list of string index.
	ids = 0
	pos = 0
	closed = []
	for letter in newick:
		if letter == '(':
			ids += 1
			parenthesis[ids] = [pos]
		elif letter == ')':
			idc = ids
			while idc in closed and idc != 0:
				idc = idc - 1
			parenthesis[idc].append(pos + 1)  # ensure the string incldues the closing parenthesis
			closed.append(idc)
		pos += 1
	# Part II: Where we use string operations to identify components parts of each split.
	inspected = []
	miss_bval = 0
	for key in parenthesis.keys():
		# This extracs splits deduced from the parenthetical
		# notation ussing mappings in dictionary P.
		r_vec = newick[parenthesis[key][0]: parenthesis[key][1]]
		#        print r_vec
		vec = sorted(get_leaves(r_vec))
		covec = sorted(complement(vec, tree.leaves))
		# Complementary splits are inferred as the set of leaves not included in the parenthesis grouping.
		if vec not in inspected and covec not in inspected:
			my_splits = Split()
			my_splits.vecs = [vec, covec]
			exp = re.escape(r_vec) + r'([0-9\.]*:[0-9\.]+)'
			branch_val = re.findall(exp, tree.newick)
			try:
				my_splits.branch_length = branch_val[0].Split(':')[1]
				my_splits.support = branch_val[0].Split(':')[0]
			except:
				miss_bval += 1
			tree.splits.append(my_splits)
			inspected.append(vec)
			inspected.append(covec)
	for leaf in leaves:  # Splits leading to each terminal are included.
		vec = [leaf]
		covec = sorted(complement(vec, leaves))
		if vec not in inspected:
			inspected.append(leaf)
			my_splits = Split()
			my_splits.vecs = [vec, covec]
			exp = re.escape(leaf) + r'\:([0-9\.]+)'
			branch_val = re.findall(exp, tree.newick)
			try:
				my_splits.branch_length = branch_val[0]
			except:
				miss_bval += 1
			tree.splits.append(my_splits)


"""
remove redundancies
"""


def no_og_subsets(file, out_file):
	"""
	Takes a UPho_orthogroups.csv.
	It writes a similar formated file with one orthologroup per line but without subsets
	"""
	out = open(out_file, 'w')
	m_list = open(file).readlines()
	total_subsets = 0
	print('Master list contains %d elements' % len(m_list))
	f = open(file, 'r')
	for line in f:
		score = 0
		a = line.strip('\n').split(',')
		for b in m_list:
			b = b.strip('\n').split(',')
			# print B
			if set(a).issubset(b) and a != b:
				score += 1
				total_subsets += 1
		if score < 1:
			out.write(line)
	out.close()
	f.close()


def remove_ip(tree_list):
	gene_list = tree_list[1:]
	sp_list = []
	gene_list_new = []
	for gene in gene_list:
		sp = gene.Split('|')[0]
		if sp not in sp_list:
			gene_list_new.append(gene)
			sp_list.append(sp)
	return gene_list_new


def no_same_og_intesec(file, no_file):
	out = open(no_file, 'w')
	f = open(file, 'r')
	current = ''  # 当前正在处理的原始树
	independent = []
	for line in f:
		a = line.strip('\n').split(',')  # 读取每行内容，并处理成列表
		pattern = re.sub("_[0-9]+$", "", a[0])  # 获取原本的树的name
		if pattern == current:  # 当处理的是同一棵树时
			for i in independent:  # 当前树中独立的子树与当前的子树相比
				if a not in independent:  # 如果某棵子树没有处理过
					if len(set(a) & set(i)) > 0:  # 如果两棵子树有交集
						independent.remove(i)  # 去除Independent中与当前子树有交集的树
						winner = max([a, i], key=len)  # 保留长的那个
						independent.append(winner)
					else:
						independent.append(a)  # 如果没有交集则加入
		else:
			for i in independent:
				gene_list = remove_ip(i)
				result = f'{i[0]},{",".join(gene_list)}\n'
				out.write(result)
			current = pattern  # 当上棵树处理完了就更新Current标记
			independent = [a]
	for i in independent:  # 把最后一个写进去
		gene_list = remove_ip(i)
		result = f'{i[0]},{",".join(gene_list)}\n'
		out.write(result)
	f.close()
	out.close()


def largest_box(lo_l):
	"""Takes a list of lists (lol) and returns a lol where no list is a subset of the others, retaining only the largest"""
	nr = []
	for list_ in lo_l:
		score = 0
		for J in lo_l:
			if set(list_).issubset(J):
				score += 1
		if score < 2:
			nr.append(list_)
	return nr


def orthologs(phylo, min_taxa, b_support):
	"""This function returns populates the list of orthologs in the PhyloClass object"""
	ortho_branch = []
	# if in-paralogs are to be included, update cost value of each terminals.
	for s in phylo.splits:
		if s.support in [None, ''] or float(s.support) >= b_support:
			for i_split in s.vecs:
				otus = spp_in_list(i_split)
				if len(set(otus)) == 1 and len(otus) > 1:  # find splits representing in-paralogs and update costs
					for leaf in i_split:
						i_cost = 1.0 / len(otus)
						if i_cost < phylo.costs[leaf]:
							phylo.costs[leaf] = i_cost
						# Reduce cost value of inparlogue copies in poportion
						# to the number of inparalogs inplied by this split.
	for s in phylo.splits:
		if s.support in [None, ''] or float(s.support) >= b_support:
			for i_split in s.vecs:
				otus = spp_in_list(i_split)
				c_count = fsum(phylo.costs[i] for i in i_split)
				#                print "%s:%f" % (','.join(Otus), cCount)
				if len(set(otus)) == c_count and c_count >= min_taxa:
					if i_split not in ortho_branch:
						ortho_branch.append(i_split)
	ortho_branch = largest_box(ortho_branch)
	phylo.ortho = ortho_branch


def aggregate_splits(small, large):
	"""
	Takes two newick like splits where small is a subset of large and returns partial newick
	incluiding the two input groupings
	"""
	aggregate = large
	contents = get_leaves("(%s)" % small)
	placeholder = contents.pop()
	aggregate = aggregate.replace("%s," % placeholder, "@@@,")
	aggregate = aggregate.replace("%s)" % placeholder, "@@@)")
	placeholder = "@@@"
	for i in contents:  # remove from aggregate all leaves in small except the placeholder
		aggregate = aggregate.replace('%s,' % i, "")
		aggregate = aggregate.replace('%s)' % i, ")")
	aggregate = aggregate.replace(placeholder, small)
	return aggregate


def sub_newick(alist, my_phylo):
	"""
	This function takes a list of leaves forming a branch and a source tree, returning the newick subtree
	"""
	relevant = []
	seed = ''
	for split in my_phylo.splits:
		for vec in split.vecs:
			if set(vec).issubset(set(alist)) and len(vec) > 0:
				if len(vec) == len(alist):
					seed = "(%s)%s:%s;" % (','.join(vec), str(split.support), str(split.branch_length))
				elif len(vec) == 1:
					rep = '%s:%s' % (vec[0], str(split.branch_length))
					relevant.append(rep)
				else:
					rep = "(%s)%s:%s" % (','.join(vec), str(split.support), str(split.branch_length))
					relevant.append(rep)
	partial = seed
	relevant = sorted(relevant, key=len, reverse=True)  # order is important
	for e in relevant:
		partial = aggregate_splits(e, partial)
	partial = re.sub('None:', ':', partial)
	partial = re.sub(':None', ':1', partial)
	return partial


def upho_main(input_trees: list, max_genome, ort_list):
	"""
	Main program execution when trees are not written
	"""
	total = 0
	for tree in input_trees:
		name = os.path.basename(tree)
		count = 0
		with open(tree, 'r') as T:
			for line in T:
				p = Myphylo(line)
				orthologs(p, max_genome, 0)
				ort_num = 0
				for group in p.ortho:
					f_name = '#%s_%d,' % (name, ort_num)
					g = ','.join(group).strip(',')
					ort_list.write(f_name + g + '\n')
					count += 1
					total += 1
					ort_num += 1
	print('Total  orthogroups found: %d' % total)


"""
Get_fasta_from_Ref
"""


# Function definitions
def fasta_to_dict(input_dir):
	records = {}
	for file in glob.glob(os.path.join(input_dir, '*.fa')):
		"""Creates a dictionary of FASTA sequences in a File, with seqIs as key to the sequences."""
		with open(file, 'r') as F:
			for Line in F:
				if Line.startswith('>'):
					seqid = Line.split(' ')[0].strip('>').strip('\n')
					seq = ''
					records[seqid] = seq
				else:
					seq = records[seqid] + Line.strip('\n')
					records[seqid] = seq.upper()
	return records


def fasta_retriever(seq_id, fasta_dict):
	"""Returns a FASTA formated record from  a seqID and a fastaDict where fasta Id is key in FastaDict"""
	try:
		seq = fasta_dict[seq_id]
		return ">%s\n%s\n" % (seq_id, seq)
	except:
		print("\x1b[1;31;40mALERT: The sequence ID:  %s  was not found in the source Fasta file.\x1b[0m" % seq_id)


def gfr_main(query, outdir, prefix, input_path):
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	else:
		print('The output dir already exist!')
	counter = 0
	seq_source = fasta_to_dict(input_path)
	handle = open(query, 'r')
	for line in handle:
		if len(line) > 0:  # do not process empty lines
			line = line.replace(' ', '')  # remove white spaces
			qlist = line.strip('\n').split(',')
			qlist = [i for i in qlist if i != ""]
			if line.startswith('#'):  # means that filenames are provided in the input this being the fisrt field in the csv.
				Name = qlist.pop(0)
				og_filename = Name.strip('#') + '.fasta'
				og_outfile = open(outdir + '/' + og_filename, 'w')
			else:
				og_filename = prefix + "_" + str(counter) + ".fasta"
				og_outfile = open(outdir + '/' + og_filename, 'w')
				counter += 1
			for seqId in qlist:
				seq = fasta_retriever(seqId, seq_source)
				try:
					og_outfile.write(seq)
				except:
					print("There is a problem retrieving the seqID: {}. "
						"Verify the seqID is the exactly same in query and source files.\n".format(seqId))
					exit(1)
			og_outfile.close()
