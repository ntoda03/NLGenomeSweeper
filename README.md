# NLGenomeSweeper
Identification of NBS-LRR genes in genomic sequences
v.1.0

This was created as part of the project HealthyGrape2 at the INRA.

nicholas.toda@inra.fr (out-of-date)

nicholas.toda@mnhn.fr (current)

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

If you use NLGenomeSweeper for published research, please cite the paper. 

Toda, N.; Rustenholz, C.; Baud, A.; Paslier, M.-C.L.; Amselem, J.; Merdinoglu, D.; Faivre-Rampant, P. [NLGenomeSweeper: A Tool for Genome-Wide NBS-LRR Resistance Gene Identification. Genes 2020, 11, 333.](https://www.mdpi.com/2073-4425/11/3/333)
 
Thanks!

## Obtaining NLGenomeSweeper
NLGenomeSweeper can be obtained from the [GitHub page](https://github.com/ntoda03/NLGenomeSweeper)

## Requirements & Installation
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

Quick example commands for setting up an environment to run NLGenomeSweeper using conda and common utilities. This requires approximately 90 gibibytes of free disk space.

```
# create a new environment and install needed pacakages from conda
conda create -y -n NLGenomeSweeper -c bioconda -c conda-forge \
    python=3.6 blast muscle samtools bedtools hmmer transdecoder openjdk
conda activate NLGenomeSweeper

# install interproscan with Panther
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.45-80.0/interproscan-5.45-80.0-64-bit.tar.gz
tar -xzf interproscan-5.45-80.0-64-bit.tar.gz
ln -s $(pwd)/interproscan-5.45-80.0/interproscan.sh $(dirname $(which python))/interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-14.1.tar.gz
tar -xzf panther-data-14.1.tar.gz -C interproscan-5.45-80.0/data/

# download NLGenomeSweeper and install
git clone https://github.com/ntoda03/NLGenomeSweeper.git
ln -s $(pwd)/NLGenomeSweeper/NLGenomeSweeper $(dirname $(which python))/NLGenomeSweeper
NLGenomeSweeper -h
```


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
                                  This is a protein fasta file of the domain sequence(s) only.
                                  Not recommended.

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
Warning! A potential structural classification is given but these must be verified by manual annotation.

*All_candidates.bed*

This is the list of positions in the genome that were identified as probable NBS-LRR genes in bed format. Coordinates refer to the position of the NB-ARC domain and not of the potential gene. A potential classification is given.

*Final_candidates.bed*

This is the list of positions in the genome that were identified as probable NBS-LRR genes in bed format that has passed the filter for the presence of LRRs. Coordinates refer to the position of the NB-ARC domain and not of the potential gene. A potential classification is given.

*Filtered_candidates.bed*

This is the list of positions in the genome that were identified as probable NBS-LRR genes in bed format but filtered out because they did not contain an LRR. Coordinates refer to the position of the NB-ARC domain and not of the potential gene. A potential classification is given.

*All_candidates.gff3*

This is the annotation gff3 file based on the results of InterProScan to predicts domains and ORFs in the at the positions of predicted NBS-LRR genes. 

