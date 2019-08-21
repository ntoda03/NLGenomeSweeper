#!/bin/bash

#
# Search a file of proteins for NB-ARC genes using an HMM profile.
#
# $1: proteins file to search in fasta format
# $2: gff annotation
# $3: NB-ARC hmm profile
# $4: output firectory of program
# $5: directory of this program
#

proteins=$1
gff=$2
hmm=$3
outputdir=$4
programdir=$5

source $programdir/functions.sh

## Find candidates based on hmm
hmmsearch --noali -E 1e-4 --tblout $outputdir/hmmersearch.tblout \
	--domtblout $outputdir/hmmersearch.domtblout $hmm $proteins
awk '!/^#/ {print $1}' $outputdir/hmmersearch.tblout > $outputdir/hmmersearch_candidates.txt

## Create a custom species profile to find candidates
hmmsearch -E 1e-60 --domtblout $outputdir/hmmersearch_highe.txt $hmm $proteins
process_custom_profile $outputdir/hmmersearch_highe.txt $proteins $outputdir customhmmerprofile
hmmsearch --noali -E 1e-4 --tblout $outputdir/custom_hmmersearch.tblout \
	--domtblout $outputdir/custom_hmmersearch.domtblout $outputdir/customhmmerprofile.hmm $proteins
awk '!/^#/ {print $1}' $outputdir/custom_hmmersearch.tblout > $outputdir/customhmmer_candidates.txt

## Combine sets of candidates
cat $outputdir/hmmersearch_candidates.txt $outputdir/customhmmer_candidates.txt \
	| sort -k1,1 |uniq > $outputdir/Annotated_candidates.txt

## Generate output files of candidates
generate_gff $outputdir/Annotated_candidates.txt $gff $outputdir/Annotated_candidates.gff3
gff3_to_bed $outputdir/Annotated_candidates.gff3 $outputdir/Annotated_candidates.bed

