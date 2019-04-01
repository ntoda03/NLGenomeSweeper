#!/bin/bash

#
# Create a custom HMM profile and consensus sequence of NB-ARC domain.
#
# $1: Pfam NB-ARC HMM profile
# $2: proteins to search in fasta format
# $3: output dir
# $4: directory of program
# $5: output prefix name
#

hmm=$1
proteins=$2
outputdir=$3
programdir=$4
prefix=$5

source $programdir/functions.sh

hmmsearch -E 1e-4 -o $outputdir/$prefix.hmmersearch_custom.txt $hmm $proteins
process_custom_profile $outputdir/$prefix.hmmersearch_custom.txt $proteins $outputdir $prefix

