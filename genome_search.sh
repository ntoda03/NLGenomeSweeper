#!/bin/bash

#
# Search a genome file for NB-ARC domain genes using consensus domain sequence(s).
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
intron=$5

source $programdir/functions.sh

function find_genome_candidates {
	outdir=$2
	cores=$3
	format=$4
	profiles=$5
	intron=$6

	makeblastdb -in $outdir/genome.fa -out $outdir/genome.fa -dbtype nucl
	if [ $format == 'prot' ]; then

	# Blast the consensus sequence against the genome to find candidates
		tblastn -outfmt 6 -query $profiles -db $outdir/genome.fa -out $outdir/genome_blast.txt \
	        	-evalue 1e-4 -num_threads $cores
	elif [ $format == 'nucl' ]; then
		blastn -outfmt 6 -query $profiles -db $outdir/genome.fa -out $outdir/genome_blast.txt \
                        -evalue 1e-4 -num_threads $cores
	fi

	# Merge candidates within 10 bp and then combine across introns (default <500bp)
	awk '{if ($9>$10){tmp=$9;$9=$10;$10=tmp} print $2,$9,$10}' $outdir/genome_blast.txt |sed 's/ /\t/g'  \
	        | bedtools sort | bedtools merge -d 10 | awk '{print $1,$2,$3,$3-$2}' |sed 's/ /\t/g' > $outdir/genome_blast.merge.txt
	R --no-save < $programdir/mergeExons.R $outdir/genome_blast.merge.txt $outdir/genome_blast.merge2.txt $intron
	awk '{printf "%s:%s-%s %s:%s-%s,len,%s\n",$1,$2,$3,$1,$2,$3,$4}' $outdir/genome_blast.merge2.txt > $outdir/genome_blast.merge_pos.txt
	extract_seq $outdir/genome_blast.merge_pos.txt $outdir/genome.fa $outdir/genome_blast.merge.fa

	# Hit must be >80% the length of the closest NB-ARC sequence
	samtools faidx $profiles
	makeblastdb -in $profiles -out $profiles -dbtype $format
        if [ $format == 'prot' ]; then
		blastx -outfmt 6 -query $outdir/genome_blast.merge.fa -db $profiles -out $outdir/genome_blast2.txt \
	        	-evalue 1e-4 -num_threads $cores -max_target_seqs 1
        elif [ $format == 'nucl' ]; then
                blastn -outfmt 6 -query $outdir/genome_blast.merge.fa -db $profiles -out $outdir/genome_blast2.txt \
                        -evalue 1e-4 -num_threads $cores -max_target_seqs 1
	fi
	join <(awk '{print $2,$1}' $outdir/genome_blast2.txt |sort -k1,1 |sed 's/,len,/ /g') <(awk '{print $1,$2}' $profiles.fai |sort -k1,1) \
	        |awk '((0.8*$4*3) < $3) {print $2}' |sed 's/:\|-/\t/g' |bedtools sort > $outdir/Unannotated_candidates.bed
	awk '{printf "%s:%s-%s\n",$1,$2,$3}' $outdir/Unannotated_candidates.bed > $outdir/Unannotated_candidates.pos.txt
	extract_seq $outdir/Unannotated_candidates.pos.txt $outdir/genome.fa $outdir/Unannotated_candidates.fa
}

# First pass identification of candidates based on consensus sequences
mkdir -p $outputdir/first_pass
cp $outputdir/pfam_NB-ARC.fa $outputdir/first_pass/
cp $genome $outputdir/first_pass/
find_genome_candidates $genome $outputdir/first_pass $cores prot $outputdir/first_pass/pfam_NB-ARC.fa $intron

# Use candidates to create species specific consensus sequences
mkdir -p $outputdir/profiler
$programdir/createProfile.sh $outputdir/first_pass/Unannotated_candidates.fa $outputdir/profiler $programdir \
	species_specific $outputdir/pfam_NB-ARC.fa nucl

# Second pass, use the species specific consensus sequences to do a second pass of candidate identification
mkdir -p $outputdir/second_pass
cp $genome $outputdir/second_pass/
find_genome_candidates $genome $outputdir/second_pass $cores prot $outputdir/profiler/species_specific_domains.fa $intron

# Combine results from first and second pass
cat $outputdir/first_pass/Unannotated_candidates.bed $outputdir/second_pass/Unannotated_candidates.bed | \
	bedtools sort | bedtools merge -d 10 > $outputdir/Unannotated_candidates.bed
awk '{printf "%s:%s-%s\n",$1,$2,$3}' $outputdir/Unannotated_candidates.bed > $outputdir/Unannotated_candidates.pos.txt
extract_seq $outputdir/Unannotated_candidates.pos.txt $genome $outputdir/Unannotated_candidates.fa
cp $outputdir/Unannotated_candidates.bed $outputdir/candidate_sites.bed

# If a protein search was run then there are already some candidate sites	
#cat $outputdir/Unannotated_candidates.bed >> $outputdir/candidate_sites.bed
#cat $outputdir/candidate_sites.bed |bedtools sort |bedtools merge > $outputdir/candidate_sites.bed.tmp
#mv $outputdir/candidate_sites.bed.tmp $outputdir/candidate_sites.bed

# Extend candidate sites to include 10 kb of flanking sequence for domain search 
get_flanking_regions $outputdir/candidate_sites.bed $genome $outputdir/Candidate_sites

