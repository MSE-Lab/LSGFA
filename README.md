# LSGFA: Domain-based large-scale prokaryotic genomic orthologous gene inference
## Introduction
LSGFA is a pan-genome pipeline developed based on Python (v3.9.16) to perform pangenomics studies on large-scale genomic datasets.

In this pipeline, the Pfam domains of all proteins are first identified using the Pfam database (http://pfam.xfam.org/), then clustered based on the similarity of their domain architectures (DAs). A sequence similarity network (SSN) will be constructed for each clustered protein sequence group, and the network structure features of the sequence groups will be analyzed to determine the existence of subgroups. Finally, the sequence groups are fully separated. For the annotation result of Pfam domain, the non-overlapping, longest domains were retained. Additionally, for sequences not annotated to the Pfam domain or for sequences where the total coverage of the Pfam domain does not exceed 60%, the DA result is recorded as None.

For similarity clustering between different DAs, the sum of the percentage of the same domain in sequence between DAs is calculated; if the same domain percentage is more than 50%, the standard is reached. Calculate the proportion of the number of sequences in two DAs that reached the standard in their respective total number of sequences is multiplied as the weight between DAs, and finally obtain the subgraph based on this network structure to get the cluster of clustered sequences.

## Installation
Download the zip file for Linux and extract its contents to a folder of your choice to complete the installation.

```shell
# Download 
unzip lsgfa.zip -d /your_path
conda env create -f LSGFA-env.yaml
```

LSGFA has the following dependencies:

### Required dependencies
* [igraph · PyPI](https://pypi.org/project/igraph/)
* [leidenalg · PyPI](https://pypi.org/project/leidenalg/)
* [pandas · PyPI](https://pypi.org/project/pandas/)
* [tqdm · PyPI](https://pypi.org/project/tqdm/)
* [diamond](https://github.com/bbuchfink/diamond)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [HMMER](http://www.hmmer.org/)
* [mcl](https://micans.org/mcl/)

### Database
* [Pfam](http://pfam.xfam.org/)

Before running this program, please download the Pfam database.

If using **hmmscan** to annotate Pfam, download the **Pfam-A.hmm.gz.**

```shell
# Download the Pfam database, as an example for v37.1
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.1/Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

If using **MMseqs2** to annotate Pfam, download the **seed** database.

```shell
# Download the Pfam database using mmseqs2
mmseqs databases Pfam-A.seed  pfam_dir/pfam  tmp --threads 10
mmseqs createindex pfam tmp -k 5 -s 7
```

### Input file format
Input the proteome files in `fasta` format ending with `faa`.

Proteome files require sequences to be named in the format `Genome_ID|Gene_ID` .

Example:

```
>GCA_018599755.1|ORF_00001
MKLEFKKSISNKIIYTLGVLFIFLFLLGYFLPIGIDKVKSLSYSQFFFSSYTVATQLGFL
LFSFVIAYFINKEYSNKNILFYKLIGDNIFTFFYKKVAVLFFECLVFIILSITIISIIYS
DFSHYLLLIILFSLVILQYILVVGTISMVSPNILISLGISIVYWIGSVILVAINKNIFGI
VAPFEASNTMYRAVEKILNNESTFMCPTEIINTVSFFVLLFIVNTIVLLLSRKRWLKIGM
```

## Usage

The main program script is `LSGFA.py`

**Generate a help page**, a quick look at the command line options available to the software

```
python LSGFA.py -h 

optional arguments: 
	-h, --help     ashow this help message
	-i INPUT_DIR, --in INPUT_DIR
	               The directory including all genome files
	-o OUTPUT_DIR, --out OUTPUT_DIR
	               Specify a output directory. default: ./
	-f             Re-perform whole process(including pfam annotation)
	-fg            Re-perform the graph search(not including pfam annotation)
	-fb            Re-perform the homology search(blast)
	--pfam         Stop at pfam annotation
	--pg           Stop at DAgraph
	-t THREADS, --threads THREADS
	               Threads. default: 8
	-db PFAM_DB    Pfam database path.
	-search {hmmscan,mmseqs-search}
                   Select a method for blasting. default=hmmscan-ssn {1,2,3}
                   Select a method for build sequence similarity network (SSN).
                   default = 3
                       1 = Reciprocal best hit (RBH)
                       2 = Specify a reciprocal hit above the threshold (identity >= 40).
                       3 = reciprocal hits above the minimum reciprocal hit threshold.
	-blast {diamond,mmseqs-search}
                   Select a method for blasting. default=mmseqs-search
	-id IDENTITY   The identity of easy-search for the sequence with DA of None.
	               [0-100], default = 40.
	-c COVERAGE    The coverage of homology easy-search for the sequence with DA of None.
	               [0-100], default = 50.
	-e EVALUE      The coverage of Pfam domain annotation.
	               [0-1], default = 1e-5.
	-inflation INFLATION  
	               Inflation (varying this parameter affects granularity)
	               [1.2-5.0], default = 1.5.
	-partition_type {1,2,3,4,5}
                   Select a method for la.partition_type.
                   default = 4
                        1 = ModularityVertexPartition
                        2 = RBConfigurationVertexPartition
                        3 = RBERVertexPartition
                        4 = CPMVertexPartition
                        5 = SurpriseVertexPartition
	-rp RESOLUTION_PARAMETER
                   [0-1.0], default = 0.9.
                   Some methods accept resolution parameters,
                   such as RBConfigurationVertexPartition, RBERVertexPartition and CPMVertexPartition. 
                   The larger the resolution_parameter, the more subgraphs will be obtained.

```
## General method of use

```shell
python LSGFA.py -i input_dir -o output_dir -t threads -db Pfam-A.hmm
```

## Output

1. The network based on the similarity of DAs can be found in the `graph` folder.

- `cc_infomation.txt`: information about DA clustering results

- `node_genes.txt`: file of network node attributes

- `pfam_graph.gml`, `pfam_graph.txt`: files for network structure

2. The subnetworks obtained from DA clustering are shown in the `graph_cc` folder.

3. Clustering results for subnetworks are available in the `homology_search' folder.

- `blast`: blast results and ABC files for subnetworks

- `sub_cc`: subgraphs obtained after mcl clustering of subnetworks

- `sub_cc_list.txt`: list of sequences within each subgraph

4. The results of processing sequences with DA of None can be found in the `none_pfam` folder.

- `none_pfam.fa`: sequences with DA of None

5. The pangenomic taxonomic position of each OG is shown in the `pangenome` folder.

6. Pfam annotation results for the genome are available in the `pfam` folder.

## License

LSGFA is free software, licensed under GPLv3.
