# NLGenomeSweeper
Identification of NBS-LRR genes in genomic sequences
v.1.0

This was created as part of the project HealthyGrape2 at the INRA.
nicholas.toda@inra.fr (out-of-date)

## Introduction
NLGenomeSweeper is a command line bash pipeline that searches a genome for NBS-LRR (NLR) disease resistance 
genes based on the presence of the NB-ARC domain using the consensus sequence of the Pfam HMM profile (PF00931) 
and class specific consensus sequences built from Vitis vinifera. This pipeline can be used with a custom 
NB-ARC HMM consensus protein sequence(s) built for a species of interest or related species for greater power,
separately for each type of NBS-LRR (TNLs, CNLs, NLs) and combine them into a single fasta file for use. This 
pipeline shows high specificity for complete genes and structurall complete pseudogenes. However, candidate regions 
are identified but may not necessarily represent functional genes and does not itself do gene prediction. A domain 
identification step is also included and the output in gff3 format can be used for manual annotation of NLR 
genes. Therefore, it is primarily for the identification of NLR genes for a genome where either no annotation 
exists or a large number of genes are expected to be absent due to repeat masking and difficulties in annotation. 
For many genomes this may be the case. 

If you use NLGenomeSweeper for published research, please cite the paper (coming soon). Thanks!

## Obtaining NLGenomeSweeper
NLGenomeSweeper can be obtained from the [GitHub page](https://github.com/ntoda03/NLGenomeSweeper)

## Requirements
This is a bash pipeline to be run on linux machines.
The following software must be available in your path.

* Python version 3.5 or greater
* blast+
* MUSCLE aligner
* SAMtools
* bedtools
* HMMER
* InterProScan
* TransDecoder

## Using this script

### Running NLGenomeSweeper
<pre>

usage: NLGenomeSweeper [-h] -genome <fasta file> [-consensus <fasta file>]
                          [-overwrite [T/F]] [-outdir <path>] [-t <threads>]

optional arguments:

  -h, --help                      show this help message and exit

Required arguments:

  -genome <fasta file>                        A genome sequence in fasta format to search.

Input arguments (optional):

  -consensus <fasta file>                     Also search the genome using a custom NB-ARC consensus sequence(s).

Output arguments (optional):

  -overwrite [T/F]                Whether to overwrite output files if they already exist. [Default F] 
  -outdir <path>                        Path where output directory NLGenomeSweeper will be created. [Default ./] 

Program control:

  -t <threads>                             The number of threads to use.

</pre>
Search for NB-ARC domains in the genome:

*NLGenomeSweeper -genome Mrot_genome.fa*

Search using custom consensus sequences:

*NLGenomeSweeper -genome Mrot_genome.fa -consensus VvinTNL.fa*

## Output
Here is an overview of the important output files created by the script. The files will be output by default in a folder called NLGenomeSweeper which will be created in the current working directory or else the user can specify and output location. 

### Primary output files
#### NLGenomeSweeper
The primary output files of NLGenomeSweeper will be written to the output directory.

*All_candidates.bed*

This is the list of positions in the genome that were identified as probable NBS-LRR genes in bed format. Coordinates refer to the position of the NB-ARC domain and not of the potential gene.

*All_candidates.gff3*

This is the annotation gff3 file based on the results of InterProScan to predicts domains and ORFs in the at the positions of predicted NBS-LRR genes. 

