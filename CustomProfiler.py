#!/usr/bin/env python

import subprocess
import sys
import argparse
import os
import shutil
#from distutils.spawn import find_executable

description = """

Written by Nicholas Toda
Last modified March 29, 2019
v1.0

Description:  This is a helper script to create a custom HMM profile and
consensus sequence based on high quality known NBS-LRR genes. It will
isolate the NB-ARC domain from a protein fasta file. It is strongly 
suggested that if multiple classes are available (TNL, CNL, NL) that
each class is built separately and the consensus sequences aggregated into
a single fasta file for use. 

For full details see: https://github.com/ntoda03/NLGenomeSweeper

"""


'''
Definition of global names for files and paths
'''
# Folder created in output dir for results
profile_folder = '00_profile_creation'
# Location of the Pfam NB-ARC hmm profile within the program directory
nbarc_location = "data/pfam_NB-ARC.hmm"


##########################################################################
#
# main
#
##########################################################################
def main(argv):
	
	# Parse command line arguments
	parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
	
	required_group = parser.add_argument_group('Required arguments')
	required_group.add_argument("-verified", metavar="<fasta file>", required=True, type=str,
		help="A protein fasta file of known NBS-LRR genes.")
	required_group.add_argument("-prefix", metavar="<prefix>", required=True, type=str,
		help="Prefix for output files.")
		
	output_group = parser.add_argument_group('Output arguments')
	output_group.add_argument("-overwrite", metavar="[T/F]",  default='F', type=str,
		help="Whether to overwrite output files if they already exist. [Default F]")
	output_group.add_argument("-outdir", metavar="<path>",  default='./NLGenomeSweeper', type=str,
		help="Set an existing output directory for program output. [Default ./NLGenomeSweeper]")
	
	args = parser.parse_args()
	program_dir = os.path.dirname(os.path.realpath(__file__)).replace(' ', '\ ')
	

	# Error handling. 
	if '/' in args.prefix:
		parser.print_help()
		raise Exception( 'ERROR: prefix should not contain a path!')
	if not os.path.isfile(args.verified):
		parser.print_help()
		raise Exception( 'ERROR: custom profile file ' + args.verified + ' does not exist!')
	if not args.outdir.endswith('/'):
		args.outdir += '/'
	if os.path.exists(args.outdir + profile_folder) and (args.overwrite == "F"):
		parser.print_help()
		raise Exception( 'ERROR: program output already exists and program is not set to overwrite!')	
	else:
		if not os.path.exists(args.outdir):
			os.mkdir(args.outdir)
		if os.path.exists(args.outdir + profile_folder):
			shutil.rmtree(args.outdir + profile_folder)
		os.mkdir(args.outdir + profile_folder)
	
	# Run program 
	subprocess.run('{}/createProfile.sh {} {} {} {} {}'.format( program_dir, program_dir + '/' + nbarc_location, 
		args.verified, args.outdir + profile_folder + '/', program_dir, args.prefix), shell=True)

if __name__ == "__main__":
   main(sys.argv[1:])
