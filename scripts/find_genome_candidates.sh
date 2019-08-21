#!/bin/bash

#
#

genome=$1
outputdir=$2
cores=$3
format=$4
profiles=$5
intron=$6
programdir=$7

outname="All_candidates"

source $programdir/scripts/functions.sh

# Blast the consensus sequence against the genome to find candidates
if [ $format == 'prot' ]; then
        tblastn -outfmt 6 -query $profiles -db $genome -out $outputdir/genome_blast.txt \
       		-evalue 1e-4 -num_threads $cores
elif [ $format == 'nucl' ]; then
        blastn -outfmt 6 -query $profiles -db $genome -out $outputdir/genome_blast.txt \
                -evalue 1e-4 -num_threads $cores
fi

# Merge candidates within 10 bp and then combine across introns (default <500bp)
awk '{if ($9>$10){tmp=$9;$9=$10;$10=tmp} print $2,$9,$10}' $outputdir/genome_blast.txt |sed 's/ /\t/g'  \
        | bedtools sort | bedtools merge -d 10 | awk '{print $1,$2,$3,$3-$2}' |sed 's/ /\t/g' > $outputdir/genome_blast.merge.txt
merge_exons $outputdir/genome_blast.merge.txt $outputdir/genome_blast.merge2.txt $intron
awk '{printf "%s:%s-%s %s:%s-%s,len,%s\n",$1,$2,$3,$1,$2,$3,$4}' $outputdir/genome_blast.merge2.txt > $outputdir/genome_blast.merge_pos.txt
extract_seq $outputdir/genome_blast.merge_pos.txt $genome $outputdir/genome_blast.merge.fa

# Hit must be >80% the length of the closest NB-ARC sequence
samtools faidx $profiles
makeblastdb -in $profiles -out $profiles -dbtype $format
if [ $format == 'prot' ]; then
        blastx -outfmt 6 -query $outputdir/genome_blast.merge.fa -db $profiles -out $outputdir/genome_blast2.txt \
        	-evalue 1e-4 -num_threads $cores -max_target_seqs 1
elif [ $format == 'nucl' ]; then
        blastn -outfmt 6 -query $outputdir/genome_blast.merge.fa -db $profiles -out $outputdir/genome_blast2.txt \
                -evalue 1e-4 -num_threads $cores -max_target_seqs 1
fi
samtools faidx $profiles
join <(awk '{print $2,$1}' $outputdir/genome_blast2.txt |sort -k1,1 |sed 's/,len,/ /g' |uniq) <(awk '{print $1,$2}' $profiles.fai |sort -k1,1) \
        |awk '((0.8*$4*3) < $3) {print $2}' |sed 's/:\|-/\t/g' |bedtools sort |uniq > $outputdir/$outname.bed
awk '{printf "%s:%s-%s\n",$1,$2,$3}' $outputdir/$outname.bed > $outputdir/$outname.pos.txt
extract_seq $outputdir/$outname.pos.txt $genome $outputdir/$outname.fa

