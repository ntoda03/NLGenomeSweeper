# NLGenomeSweeper
Identification of NBS-LRR genes in genomic sequences

## Introduction
NLGenomeSweeper is a command line pipeline that searches a genome for NBS-LRR genes based on the presence of the NB-ARC domain using the consensus sequence of the Pfam HMM profile (PF00931). It can be used with a custom NB-ARC HMM profile and consensus protein sequence(s) built for a species of interest or related species for greater power. It is highly recommended to use a custom consensus sequence and if possible to build a consensus sequence separately for each type of NBS-LRR (TNLs, CNLs, NLs) and combine them into a single fasta file for use. This pipeline shows high specificity for functional genes and mostly intact pseudogenes. However, candidate regions are identified but may not necessarily represent functional genes and does not itself do gene prediction. A domain identification step is also included and the output in gff3 format can be used for manual annotation of NBS-LRR genes. Therefore, it is primarily for the identification of NBS-LRR genes for a genome where either no annotation exists or a large number of genes are expected to be absent due to repeat masking and difficulties in annotation. For many genomes this is expected to be the case. 

If you use NLGenomeSweeper for published research, please cite the paper (coming soon). Thanks!

## Obtaining NLGenomeSweeper
NLGenomeSweeper can be obtained from the [GitHub page](https://github.com/ntoda03/NLGenomeSweeper). It is ready to use out of the box and requires no setup.

## Requirements
The following software must be installed and available in the path.

* Python version 3.5 or newer with Pandas
* blast+
* MUSCLE aligner
* SAMtools
* bedtools
* HMMER
* InterProScan
* R
* Exonerate

## Using this script

### Running CustomProfiler
<pre>

usage: CustomProfiler.py [-h] -verified <fasta file> -prefix <prefix>
                         [-overwrite [T/F]] [-outdir <path>]

optional arguments:

  -h, --help            	show this help message and exit

Required arguments:

  -verified <fasta file>	A protein fasta file of known NBS-LRR genes.
  -prefix <prefix>      	Prefix for output files. <br>

Output arguments:

  -overwrite [T/F]      	Whether to overwrite output files if they already exist. [Default F]
  -outdir <path>        	Set an output directory for program output. [Default ./NLGenomeSweeper]


</pre>
### Running NLGenomeSweeper
<pre>

usage: NLGenomeSweeper.py [-h] -genome <fasta file> [-protein <fasta file>]
                          [-gff <gff3 file>] [-consensus <fasta file>]
                          [-hmm <file>] [-overwrite [T/F]] [-outdir <path>]
                          [-t <threads>]

optional arguments:

  -h, --help            	show this help message and exit

Required arguments:

  -genome <fasta file>  	A genome sequence in fasta format to search.

Input arguments:

  -protein <fasta file>		Protein sequences in fasta format. 
  -gff <gff3 file>      	gff3 annotation of the genome. Required if searching protein sequences. 
  -consensus <fasta file>	Also search the genome using a custom NB-ARC consensus sequence(s).
  -hmm <file>           	Also search the protein file using a custom NB-ARC HMM. <br>

Output arguments:

  -overwrite [T/F]     		Whether to overwrite output files if they already exist. [Default F] 
  -outdir <path>        	Set an output directory for program output. [Default ./NLGenomeSweeper] 

Program control:

  -t <threads>          	The number of threads to use.

</pre>
### Examples
Generate a custom HMM profile and consensus sequence:

*python CustomProfiler.py -verified vvinifera_TNLgenes.fa -prefix VvinTNL*

Search for NB-ARC domains in the genome:

*python NLGenomeSweeper.py -genome Mrot_genome.fa*

Search for NB-ARC domains in proteins and the genome:

*python NLGenomeSweeper.py -genome Mrot_genome.fa -protein Mrot_protein.fa Mrot_annot.gff3*

Search using a custom HMM profile and consensus sequence:

*python NLGenomeSweeper.py -genome Mrot_genome.fa -consensus VvinTNL.fa -hmm VvinTNL.hmm*

## Output
Here is an overview of the important output files created by the script. The files will be output by default in a folder called NLGenomeSweeper which will be created in the current working directory or else the user can specify and output location. 

### Primary output files
The primary output files of CustomProfiler will be written to the folder 00_profile_creation within the output directory.

*00_profile_creation/\<prefix\>.fa*

The consensus NB-ARC domain sequence generated for the known NBS-LRR genes.

*00_profile_creation/\<prefix\>.hmm*

The NB-ARC domain HMM profile generated for the known NBS-LRR genes.


The primary output files of NLGenomeSweeper will be written to the output directory.

*Annotated_candidates.bed*

If a protein file is supplied to be searched then the script will output a list of candidate proteins identified as NBS-LRR candidate genes in bed format.

*Unannotated_candidates.bed*

This is the list of positions in the genome that were identified as probable NBS-LRR genes in bed format. Coordinates refer to the position of the NB-ARC domain and not of the potential gene.

*All_candidates_regions.annotated.gff3*

This is the annotation gff3 file based on the results of InterProScan to predicts domains and ORFs in the at the positions of predicted NBS-LRR genes. 

### Other output files of interest
If a protein file is provided to search for NBS-LRR genes then the script will also compute a species specific NB-ARC domain profile and consensus sequence using the candidates identified. These may be of interest and can be found in the following location. If a custom HMM profile is provided to the program then a second set of species specific NB-ARC profile and sequence will also be generated.

*01_candidate_identification/proteins_pfam_search/customhmmerprofile.hmm*
*01_candidate_identification/proteins_pfam_search/customhmmerprofile.fa*

Species specific HMM profile and consensus sequence generated using the pfam NB-ARC hmm.

*01_candidate_identification/proteins_custom_search/customhmmerprofile.hmm*
*01_candidate_identification/proteins_custom_search/customhmmerprofile.fa*

Species specific HMM profile and consensus sequence generated using a custom NB-ARC hmm provided by the user.
