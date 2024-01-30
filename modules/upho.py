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
class split:
	def __init__(self):
		self.vecs = None
		self.branch_length = None
		self.support = None
		self.name = None
		self.label = []


class myPhylo:
	"""A class for newick trees"""

	def __init__(self, N):
		self.leaves = get_leaves(N)
		self.splits = []
		self.ortho = []
		self.costs = {}  # Dictionary of leaf cost for inparalog evaluation'
		self.newick = N.strip('\n')
		for leaf in self.leaves:  # Cost initialized
			self.costs[leaf] = 1.0
		split_decomposition(self)

# FUNCTION DEFINITIONS


def get_leaves(String):
	"""Find leaves names in newick files using regexp.
	Leaves names are composed of alpha numeric characters,
	underscore and a special field delimiter"""
	pattern = "(?<=[,\(])\w.*?(?=[,:;\)])"
	#    pattern = "[^\(\),;:\[\]]+%s[^\(\),;:\[\]]+" % gsep
	Leaves = re.findall(pattern, String)
	return Leaves


def spp_in_list(alist):
	"""Return the species from a list of sequence identifiers"""
	spp = []
	for i in alist:
		spp.append(i.split(sep)[0])
	return spp


def complement(Sub, Whole):
	"""Return elements in Whole that are not present in Sub"""
	complement = set(Whole) - set(Sub)
	return list(complement)


def split_decomposition(Tree):
	"""Add a list of splits class objects to myPhylo Class object"""
	# Part I. Where we identify matching parenthesis in the newick.
	newick = Tree.newick
	leaves = Tree.leaves
	P = {}  # Empty dictionary where to store parenthesis identifiers and a list of string index.
	ids = 0
	idc = 0
	Pos = 0
	closed = []
	for l in newick:
		if l == '(':
			ids += 1
			P[ids] = [Pos]
		elif l == ')':
			idc = ids
			while idc in closed and idc != 0:
				idc = idc - 1
			P[idc].append(Pos + 1)  # ensure the string incldues the closing parenthesis
			closed.append(idc)
		Pos += 1
	# Part II: Where we use string operations to identify components parts of each split.
	Inspected = []
	missBval = 0
	for Key in P.keys():
		# This extracs splits deduced from the parenthetical
		# notation ussing mappings in dictionary P.
		r_vec = newick[P[Key][0]: P[Key][1]]
		#        print r_vec
		vec = sorted(get_leaves(r_vec))
		covec = sorted(complement(vec,
								  Tree.leaves))  # Complementary splits are inferred as the set of leaves not included in the parenthesis grouping.
		if vec not in Inspected and covec not in Inspected:
			mySplits = split()
			mySplits.vecs = [vec, covec]
			exp = re.escape(r_vec) + r'([0-9\.]*:[0-9\.]+)'
			BranchVal = re.findall(exp, Tree.newick)
			try:
				mySplits.branch_length = BranchVal[0].split(':')[1]
				mySplits.support = BranchVal[0].split(':')[0]
			except:
				missBval += 1
			Tree.splits.append(mySplits)
			Inspected.append(vec)
			Inspected.append(covec)
	for leaf in leaves:  # Splits leading to each terminal are included.
		vec = [leaf]
		covec = sorted(complement(vec, leaves))
		if vec not in Inspected:
			Inspected.append(leaf)
			mySplits = split()
			mySplits.vecs = [vec, covec]
			exp = re.escape(leaf) + r'\:([0-9\.]+)'
			BranchVal = re.findall(exp, Tree.newick)
			try:
				mySplits.branch_length = BranchVal[0]
			except:
				missBval += 1
			Tree.splits.append(mySplits)

"""
remove redundancies
"""

def No_OG_subsets(File, out_file):
	'''Takes a UPho_orthogroups.csv. It writes a similar formated file with one orthologroup per line but without subsets'''
	Out = open(out_file, 'w')
	M_List = open(File).readlines()
	F = open(File, 'r')
	TotalSubsets = 0
	print('Master list contains %d elements' % len(M_List))
	F = open(File, 'r')
	for Line in F:
		Score = 0
		A = Line.strip('\n').split(',')
		for B in M_List:
			B = B.strip('\n').split(',')
			# print B
			if set(A).issubset(B) and A != B:
				Score += 1
				TotalSubsets += 1
		if Score < 1:
			Out.write(Line)
	Out.close()
	F.close()

def remove_ip(tree_list):
	pat_ = tree_list[0]
	gene_list = tree_list[1:]
	sp_list = []
	for gene in gene_list:
		sp = gene.split('|')[0]
		if sp in sp_list:
			gene_list.remove(gene)
		else:
			sp_list.append(sp)
	return gene_list

def No_Same_OG_Intesec(File, no_file):
	Out = open(no_file, 'w')
	F = open(File, 'r')
	Current = ''  # 当前正在处理的原始树
	Independent = []
	for Line in F:
		A = Line.strip('\n').split(',')  # 读取每行内容，并处理成列表
		Pattern = re.sub("_[0-9]+$", "", A[0])  # 获取原本的树的name
		if Pattern == Current:  # 当处理的是同一棵树时
			for i in Independent:  # 当前树中独立的子树与当前的子树相比
				if A not in Independent:  # 如果某棵子树没有处理过
					if len(set(A) & set(i)) > 0:  # 如果两棵子树有交集
						Independent.remove(i)  # 去除Independent中与当前子树有交集的树
						Winner = max([A, i], key=len)  # 保留长的那个
						Independent.append(Winner)
					else:
						Independent.append(A)  # 如果没有交集则加入
		else:
			for i in Independent:
				gene_list = remove_ip(i)
				result = f'{i[0]},{",".join(gene_list)}\n'
				Out.write(result)
			Current = Pattern  # 当上棵树处理完了就更新Current标记
			Independent = []
			Independent.append(A)
	for i in Independent:  # 把最后一个写进去
		gene_list = remove_ip(i)
		result = f'{i[0]},{",".join(gene_list)}\n'
		Out.write(result)
	F.close()
	Out.close()


def LargestBox(LoL):
	'''Takes a list of lists (lol) and returns a lol where no list is a subset of the others, retaining only the largest'''
	NR = []
	for L in LoL:
		score = 0
		for J in LoL:
			if set(L).issubset(J):
				score += 1
		if score < 2:
			NR.append(L)
	return NR


def orthologs(Phylo, minTaxa, bsupport):
	"""This function returns populates the list of orthologs in the PhyloClass object"""
	OrthoBranch = []
	# if in-paralogs are to be included, update cost value of each terminals.
	for S in Phylo.splits:
		if S.support in [None, ''] or float(S.support) >= bsupport:
			for i_split in S.vecs:
				Otus = spp_in_list(i_split)
				if len(set(Otus)) == 1 and len(Otus) > 1:  # find splits representing in-paralogs and update costs
					for leaf in i_split:
						ICost = 1.0 / len(Otus)
						if ICost < Phylo.costs[leaf]:
							Phylo.costs[leaf] = ICost  # Reduce cost value of inparlogue copies in poportion to the number of inparalogs inplied by this split.
	for S in Phylo.splits:
		if S.support in [None, ''] or float(S.support) >= bsupport:
			for i_split in S.vecs:
				Otus = spp_in_list(i_split)
				cCount = fsum(Phylo.costs[i] for i in i_split)
				#                print "%s:%f" % (','.join(Otus), cCount)
				if len(set(Otus)) == cCount and cCount >= minTaxa:
					if i_split not in OrthoBranch:
						OrthoBranch.append(i_split)
	OrthoBranch = LargestBox(OrthoBranch)
	Phylo.ortho = OrthoBranch


def aggregate_splits(small, large):
	"""Takes two newick like splits where small is a subset of large and returns partial newick incluiding the two input groupings"""
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


def subNewick(alist, myPhylo):
	'''This function takes a list of leaves forming a branch and a source tree, returning the newick subtree'''
	relevant = []
	seed = ''
	for split in myPhylo.splits:
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

def upho_main(input_trees:list, max_genome, OrtList):
	'''Main program execution when trees are not written'''
	Total = 0
	for tree in input_trees:
		name=os.path.basename(tree)
		count = 0
		with open(tree, 'r') as T:
			for line in T:
				P = myPhylo(line)
				orthologs(P, max_genome, 0)
				ortNum = 0
				for group in P.ortho:
					FName = '#%s_%d,' % (name, ortNum)
					G = ','.join(group).strip(',')
					OrtList.write(FName + G + '\n')
					count += 1
					Total += 1
					ortNum += 1
	print('Total  orthogroups found: %d' % Total)


"""
Get_fasta_from_Ref
"""

# Function definitions
def Fasta_to_Dict(input_dir):
	Records = {}
	for file in glob.glob(os.path.join(input_dir, '*.faa')):
		"""Creates a dictionary of FASTA sequences in a File, with seqIs as key to the sequences."""
		with open(file, 'r') as F:
			for Line in F:
				if Line.startswith('>'):
					Seqid = Line.split(' ')[0].strip('>').strip('\n')
					Seq = ''
					Records[Seqid] = Seq
				else:
					Seq = Records[Seqid] + Line.strip('\n')
					Records[Seqid] = Seq.upper()
	return Records


def FastaRetriever(seqId, FastaDict):
	"""Returns a FASTA formated record from  a seqID and a fastaDict where fasta Id is key in FastaDict"""
	try:
		seq = FastaDict[seqId]
		return ">%s\n%s\n" % (seqId, seq)
	except:
		print("\x1b[1;31;40mALERT: The sequence ID:  %s  was not found in the source Fasta file.\x1b[0m" % seqId)


def gfr_main(query, outdir, prefix, input_path):
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	else:
		print('The output dir already exist!')
	Counter = 0
	seqSource = Fasta_to_Dict(input_path)
	handle = open(query, 'r')
	for line in handle:
		if len(line) > 0:  # do not process empty lines
			line = line.replace(' ', '')  # remove white spaces
			qlist = line.strip('\n').split(',')
			qlist = [i for i in qlist if i != ""]
			if line.startswith('#'):  # means that filenames are provided in the input this being the fisrt field in the csv.
				Name = qlist.pop(0)
				OG_filename = Name.strip('#') + '.fasta'
				OG_outfile = open(outdir + '/' + OG_filename, 'w')
			else:
				OG_filename = prefix + "_" + str(Counter) + ".fasta"
				OG_outfile = open(outdir + '/' + OG_filename, 'w')
				Counter += 1
			for seqId in qlist:
				seq = FastaRetriever(seqId, seqSource)
				try:
					OG_outfile.write(seq)
				except:
					print("There is a problem retrieving the seqID: %s. Verify the seqID is the exactly same in query and source files.\n" % seqId)
					exit(1)
			OG_outfile.close()
