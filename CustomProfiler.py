#!/usr/bin/env python

import subprocess
import sys
import argparse
import os
import shutil
#from distutils.spawn import find_executable

description = """

Written by Nicholas Toda
v1.2.2

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

##########################################################################
#
# main
#
##########################################################################
def main(argv):
    program_dir = os.path.dirname(os.path.realpath(__file__)).replace(' ', '\ ')

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    required_group = parser.add_argument_group('Required arguments')
    required_group.add_argument("-verified", metavar="<fasta file>", required=True, type=str,
        help="A fasta file of protein sequences of known NBS-LRR genes.")
    required_group.add_argument("-prefix", metavar="<prefix>", required=True, type=str,
        help="Prefix for output files.")
        
    consensus_group = parser.add_argument_group('Domain arguments')
    consensus_group.add_argument("-reference_nbarc", metavar="<fasta file>",  
        default=program_dir + '/' + 'data/Vitis_vinifera_NB-ARC_consensus.fa', type=str,
        help="A reference NB-ARC domain protein sequence for comparison in fasta format. [Default data/Vitis_vinifera_NB-ARC_consensus.fa]")

    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument("-overwrite", metavar="[T/F]",  default='F', type=str,
        help="Whether to overwrite output files if they already exist. [Default F]")
    output_group.add_argument("-outdir", metavar="<path>",  default='./NLGenomeSweeper', type=str,
        help="Set an existing output directory for program output. [Default ./NLGenomeSweeper]")
    
    args = parser.parse_args()
    

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
    subprocess.run('{}/scripts/createProfile.sh {} {} {} {} {} prot'.format( program_dir, args.verified,  
         args.outdir + profile_folder + '/', program_dir, args.prefix, args.reference_nbarc, "prot"), shell=True)

if __name__ == "__main__":
   main(sys.argv[1:])
