#!/usr/bin/env python

import subprocess
import sys
import argparse
import os
import shutil
import pandas as pd

description = """

Written by Nicholas Toda
Last modified March 29, 2019
v1.0

Description:  Searches a genome for NBS-LRR genes based on the presence
of the NB-ARC domain using the consensus sequence of the Pfam HMM profile 
(PF00931) and can be used with custon domain profiles built for a species 
of interest or related species for greater power. 

For full details see: https://github.com/ntoda03/NLGenomeSweeper

"""


'''
Definition of global names for files and paths
'''
# Folder created in output dir for candidate search results
search_folder = '01_candidate_identification/'
# subfolder within the search folder for searching proteins with hmmer
protein_search_folder = 'proteins_pfam_search/'
# subfolder within the search folder for searching proteins with a custom hmm profile
protein_custom_folder = 'proteins_custom_search/'
# subfolder within the search folder for searching the genome with blast
genome_search_folder = 'genome_search/'
# the name of the species specific consensus file generated when provided with protein file
species_consensus = 'customhmmerprofile.fa'

# Folder created in output dir for domain identification results
domain_folder = '02_domain_identification/'

# Location of the Pfam NB-ARC hmm profile and consensus within the program directory
nbarc_hmm = "data/pfam_NB-ARC.hmm"
nbarc_fa = "data/pfam_NB-ARC.fa"

# Output file names
annotated_candidates = 'Annotated_candidates.bed'
unannotated_candidates = 'Unannotated_candidates.bed'
annotated_gff= 'All_candidates_regions.annotated.gff3'

# Append the contents of an input file to an output file
def append_file(input_file, output_file):
	infile = open(input_file, 'r') 
	appendFile = open(output_file,'a')
	appendFile.write('\n')
	appendFile.write(infile.read())
	appendFile.close()

##########################################################################
#
# identify_candidates_proteins
#
# Arguments:
# program_dir - directory of this program
# protein - protein file to search (fasta format)
# gff - gff3 annotation of the genome
# hmm - custom NB-ARC hmm profile to also use for search (or None)
# outdir - output directory
#
# Output:
# outdir/Annotated_candidates.bed - a bed file of identified candidates
#
# This function calls bash scripts to search a protein file in fasta
# format for NB-ARC domain containing proteins using hmmer. The search is
# always done using the pfam NB-ARC domain hmm profile and additionally
# a custom NB-ARC profile can be optimally used. 
# 
##########################################################################
def identify_candidates_proteins(program_dir, protein, gff, hmm, outdir):
	# Searching using the pfam NB-ARC hmm profile
	os.mkdir(outdir + search_folder + '/' + protein_search_folder)
	subprocess.run('{}/hmmer_candidates.sh {} {} {} {} {}'.format(program_dir, protein, gff,
		program_dir + nbarc_hmm,
		outdir + search_folder + protein_search_folder, program_dir), shell=True)
	shutil.copy(outdir + search_folder + protein_search_folder + annotated_candidates, 
		outdir + annotated_candidates)
	
	# Also searching using a custom hmm profile
	if hmm: 
		os.mkdir(outdir  + search_folder + protein_custom_folder)
		subprocess.run('{}/hmmer_candidates.sh {} {} {} {} {}'.format(program_dir,protein,gff,
			hmm,outdir + search_folder + protein_custom_folder, program_dir), shell=True)
		# Concatenate to the pfam search results
		#annotated1 = pd.read_csv(outdir + search_folder + protein_custom_folder +annotated_candidates, 
		#	sep='\t',header=None)
		#annotated2 = pd.read_csv(outdir + '/Annotated_candidates.bed', sep='\t',header=None)
		#annotated = annotated1.append(annotated2, ignore_index=True).drop_duplicates()
		#annotated.to_csv(outdir + '/Annotated_candidates.bed', sep='\t', header =False, index = False )
		append_file(outdir + search_folder + protein_custom_folder + annotated_candidates, 
			outdir + annotated_candidates)

##########################################################################
#
# identify_candidates_genome
#
# Arguments:
# program_dir - directory of this program
# genome - genome file to search (fasta format)
# consensus - consensus protein sequence(s) of NB-ARC domain in fasta format
# outdir - output directory
# cores - number of cores to use
#
# Output:
# outdir/Unannotated_candidates.bed - a bed file of identified candidates

# This function calls bash scripts to search a genome file in fasta
# format for NB-ARC domain containing sequences using tblast. The search is
# always done using the pfam NB-ARC domain consensus sequence and additionally
# a custom NB-ARC consensus protein sequence can be used. The custom file
# can contain multiple sequences and it is strongly recommended to create
# multiple consensus sequence for each type of NBS-LRR gene (TNL, CNL, NL).
# 
##########################################################################
def identify_candidates_genome(program_dir, genome, consensus, outdir, cores):
	os.mkdir(outdir + search_folder + genome_search_folder)
	shutil.copy(program_dir + nbarc_fa, 
		outdir + search_folder + genome_search_folder + 'pfam_NB-ARC.fa')
	# Also searching using a custom consensus sequence, add to query file
	if consensus:
		append_file(consensus, outdir + search_folder + genome_search_folder + 'pfam_NB-ARC.fa')
	# Also searched a protein file so a species specific NB-ARC consensus sequence was generated  
	if os.path.isfile(outdir + annotated_candidates + protein_search_folder + species_consensus):
		append_file(outdir + annotated_candidates + protein_search_folder + species_consensus, 
			outdir + search_folder + genome_search_folder + 'pfam_NB-ARC.fa')
	# Also searched a protein file using a custom hmm so a second species specific NB-ARC 
	# consensus sequence was generated 
	if os.path.isfile(outdir + annotated_candidates + protein_custom_folder + species_consensus):
		append_file(outdir + annotated_candidates + protein_custom_folder + species_consensus,
			outdir + search_folder + genome_search_folder + 'pfam_NB-ARC.fa')
	# Already ran protein based search, mask genome to avoid duplication
	if os.path.isfile(outdir + annotated_candidates):
		subprocess.run('bedtools maskfasta -fi {} -fo {} -bed {}'.format(
			genome, outdir + search_folder + genome_search_folder + 'genome.masked.fa', 
			outdir + annotated_candidates), shell=True)
		genome = outdir + search_folder + genome_search_folder + 'genome.masked.fa'
		# Remember to also include the annotated candidates for extracting sequences
		shutil.copy(outdir + annotated_candidates, 
			outdir + search_folder + genome_search_folder + 'candiates_sites.bed')
	# Run the genome based search
	subprocess.run('{}/genome_search.sh {} {} {} {}'.format(program_dir,genome, 
		outdir + search_folder + genome_search_folder, program_dir, cores), shell=True)
	shutil.copy(outdir + search_folder + genome_search_folder + unannotated_candidates, 
		outdir + unannotated_candidates)
	shutil.copy(outdir + search_folder + genome_search_folder + 'Candidate_sites.with_flanking.fa', 
		outdir + domain_folder + 'Candidate_sites.with_flanking.fa')


##########################################################################
#
# run_interproscan
#
# Arguments:
# program_dir - directory of this program
# outdir - output directory
# cores - number of cores to use
#
# Output:
# outdir/All_candidates_regions.annotated.gff3 - domain annotation of candidate sites

# This function calls a bash script to search the identified NBS-LRR genes 
# candidate sequences for ORFs and domains using interproscan. 10 kb of flanking 
# sequence on both sides are also searched. This returns an annotation file 
# in gff3 format that can be viewed directly in a genome browser that contains
# relevant information on the ORFs and domains at the candidate sites.
# 
##########################################################################
def run_interproscan(program_dir, outdir, t):
	subprocess.run('{}/region_interproscan.sh {} {}'.format(program_dir, outdir + domain_folder, t),
		shell=True)
	shutil.copy(outdir + domain_folder + annotated_gff, 
		outdir + annotated_gff)

##########################################################################
#
# main
#
##########################################################################
def main(argv):
	
	# Parse command line arguments.
	parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
	
	required_group = parser.add_argument_group('Required arguments')
	required_group.add_argument("-genome", metavar="<fasta file>", required=True, type=str,
		help="A genome sequence in fasta format to search.")
	
	input_group = parser.add_argument_group('Input arguments')
	input_group.add_argument("-protein", metavar="<fasta file>",  type=str, default = None,
		help="Protein sequences in fasta format.")
	input_group.add_argument("-gff", metavar="<gff3 file>", type=str, default = None,
		help="gff3 annotation of the genome. Required if searching protein sequences.")
	input_group.add_argument("-consensus", metavar="<fasta file>", type=str, default = None,
		help="Also search the genome using a custom NB-ARC consensus sequence(s).")
	input_group.add_argument("-hmm", metavar="<file>",  type=str, default = None,
		help="Also search the protein file using a custom NB-ARC HMM.")
	
	output_group = parser.add_argument_group('Output arguments')
	output_group.add_argument("-overwrite", metavar="[T/F]",  default='F', type=str,
		help="Whether to overwrite output files if they already exist. [Default F]")
	output_group.add_argument("-outdir", metavar="<path>",  default='./NLGenomeSweeper', type=str,
		help="Set an existing output directory for program output. [Default ./NLGenomeSweeper]")
	
	run_group = parser.add_argument_group('Program control')
	run_group.add_argument("-t", metavar="<threads>",  type=int, default = 1,
		help="The number of threads to use.")
		
	args = parser.parse_args()
	program_dir = os.path.dirname(os.path.realpath(__file__)).replace(' ', '\ ')
	
	# Error handling. 
	if not os.path.isfile(args.genome):
		parser.print_help()
		raise Exception( 'ERROR: genome file ' + args.genome + ' does not exist!')
	if args.protein and not os.path.isfile(args.protein):
		parser.print_help()
		raise Exception( 'ERROR: protein file ' + args.protein + ' does not exist!')
	if args.protein and not args.gff:
		parser.print_help()
		raise Exception( 'ERROR: gff file must be supplied when searching proteins!')
	if args.gff and not os.path.isfile(args.gff):
		parser.print_help()
		raise Exception( 'ERROR: gff file ' + args.gff + ' does not exist!')
	if args.consensus and not os.path.isfile(args.consensus):
		parser.print_help()
		raise Exception( 'ERROR: consensus file ' + args.consensus + ' does not exist!')
	if args.hmm and not os.path.isfile(args.hmm):
		parser.print_help()
		raise Exception( 'ERROR: HMM file ' + args.hmm + ' does not exist!')
	if not args.outdir.endswith('/'):
		args.outdir += '/'
	if not program_dir.endswith('/'):
		program_dir += '/'
	if os.path.exists(args.outdir + search_folder) and (args.overwrite == "F"):
		parser.print_help()
		raise Exception( 'ERROR: program output already exists and program is not set to overwrite!')
	else:
		if not os.path.exists(args.outdir):
			os.mkdir(args.outdir)
		if os.path.exists(args.outdir + search_folder):
			shutil.rmtree(args.outdir + search_folder)
		os.mkdir(args.outdir + search_folder)
		if os.path.exists(args.outdir + domain_folder):
			shutil.rmtree(args.outdir + domain_folder)
		os.mkdir(args.outdir + domain_folder)
	
	# If a protein file was given to search
	if args.protein:
		# Search protein file using hmmer
		identify_candidates_proteins(program_dir,args.protein, args.gff, args.hmm, args.outdir )
	# Search the genome using blast
	identify_candidates_genome(program_dir, args.genome, args.consensus, args.outdir, args.t)
	# Define ORFs and domains using interproscan
	run_interproscan(program_dir, args.outdir, args.t)


if __name__ == "__main__":
   main(sys.argv[1:])
