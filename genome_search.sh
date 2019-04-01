#!/bin/bash

#
# Search a genome file for NB-ARC domain genes using a consensus domain sequence.
#
# $1: genome file in fasta format to search
# $2: output directory
# $3: directory of this program
# $4: number of cores to use
#

genome=$1
outputdir=$2
programdir=$3
cores=$4

source $programdir/functions.sh

cp $genome $outputdir/genome.fa
makeblastdb -in $outputdir/genome.fa -out $outputdir/genome.fa -dbtype nucl

# Blast the consensus sequence against the genome to find candidates
tblastn -outfmt 6 -query $outputdir/pfam_NB-ARC.fa -db $outputdir/genome.fa -out $outputdir/genome_blast.txt \
	-evalue 1e-4 -num_threads $cores
	
# Merge candidates on the same strand within 500 bp
R --no-save < $programdir/mergeBlastHits.R $outputdir/genome_blast.txt $outputdir/genome_blast.merge.txt \
	500 same flanking

# Hit must be >80% the length of the NB-ARC sequence
samtools faidx $outputdir/pfam_NB-ARC.fa
join $outputdir/genome_blast.merge.txt <(awk '{print $1,$2}' $outputdir/pfam_NB-ARC.fa.fai |sort -k1,1) | \
awk '($4> (0.8*$13)) {if($9>$10){tmp=$9;$9=$10;$10=tmp};print $2,$9,$10}' | \
sed 's/ /\t/g' |bedtools sort |bedtools merge > $outputdir/Unannotated_candidates.bed

cat $outputdir/Unannotated_candidates.bed >> $outputdir/candidate_sites.bed
cat $outputdir/candidate_sites.bed |bedtools sort |bedtools merge > $outputdir/candidate_sites.bed.tmp
mv $outputdir/candidate_sites.bed.tmp $outputdir/candidate_sites.bed

# Extend candidate sites to include 10 kb of flanking sequence for domain search 
get_flanking_regions $outputdir/candidate_sites.bed $genome $outputdir/Candidate_sites