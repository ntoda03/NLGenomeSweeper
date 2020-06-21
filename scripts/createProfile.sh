#!/bin/bash

#
# Create custom consensus sequences based on domains desired.
#
# This program takes potential domains identified in the genome and compares
# them to known high quality consnesus sequences of that domain. For each
# potential domain the closest reference sequence if found and sequences are
# grouped by their closest reference. The grouped sequences are aligned
# and a protein consensus sequence is created for each group.
# 
# Inputs:
# $1: Fasta file of sequences containing domain extracted from the genome
# $2: output directory
# $3: directory of this program
# $4: prefix of output file
# $5: Consensus sequences in fasta format, protein or nucleotide
# $6: Format of the consensus sequences, 'prot' or 'nucl'
#
# Outputs:
# $outputdir/$prefix.fa: Fasta file of species specific consensus sequences, proteins
# 

#hmm=$1
sequences=$1
outputdir=$2
programdir=$3
prefix=$4
consensus=$5
format=$6

source $programdir/scripts/functions.sh

# Blast the consensus domain sequences against the sequences to search and extract the overlapped region 
cp $sequences $outputdir/sequences_cluster_search.fa
sed -i 's/:\|-/./g' $outputdir/sequences_cluster_search.fa
makeblastdb -in $outputdir/sequences_cluster_search.fa -out $outputdir/sequences_cluster_search.fa -dbtype $format
if [ $format == 'prot' ]; then
	blastp -outfmt 6 -query $consensus -db $outputdir/sequences_cluster_search.fa -out $outputdir/consensus_blast.txt -evalue 1e-4
elif [ $format == 'nucl' ]; then
	tblastn -outfmt 6 -query $consensus -db $outputdir/sequences_cluster_search.fa -out $outputdir/consensus_blast.txt -evalue 1e-4
fi
awk '{if ($9>$10){tmp=$9;$9=$10;$10=tmp} print $2,$9,$10}' $outputdir/consensus_blast.txt |sed 's/ /\t/g' | \
        bedtools sort | bedtools merge -d 50 | \
	awk '$3 - $2 > 200 {printf "%s:%s-%s\n",$1,$2,$3}'  > $outputdir/consensus.list.txt
extract_seq $outputdir/consensus.list.txt $outputdir/sequences_cluster_search.fa $outputdir/consensus.fa
sed -i 's/:[^ ]*//g' $outputdir/consensus.fa

# Blast the extracted sequence against the consensus to cluster based on similarity to them
cp $consensus $outputdir/consensus_cluster_search.fa
makeblastdb -in $outputdir/consensus_cluster_search.fa -out $outputdir/consensus_cluster_search.fa -dbtype prot
if [ $format == 'prot' ]; then
	blastp -outfmt 6 -query $outputdir/consensus.fa -db $outputdir/consensus_cluster_search.fa -out \
		$outputdir/consensus_blast_clustering.txt -evalue 1e-4 -max_target_seqs 1
elif [ $format == 'nucl' ]; then
	blastx -outfmt 6 -query $outputdir/consensus.fa -db $outputdir/consensus_cluster_search.fa -out \
                $outputdir/consensus_blast_clustering.txt -evalue 1e-4 -max_target_seqs 1
fi
awk -v output=$outputdir '{print $1 >> ""output"/consensus_clusters."$2".txt" }' $outputdir/consensus_blast_clustering.txt

# For each cluster, create a consensus sequence and profile (must be >200 aa or >600 bp)
for file in $outputdir/consensus_clusters.*.txt
do
    basename=${file##*/}
    extract_seq $file $outputdir/consensus.fa $file.fa
	if [ $format == 'nucl' ]; then
		(cd $outputdir && TransDecoder.LongOrfs -m 50 -t $basename.fa)
		(cd $outputdir && TransDecoder.Predict --no_refine_starts --single_best_only -t $basename.fa)	
		create_custom_profile $file.fa.transdecoder.pep $outputdir $basename.profile
	else
		create_custom_profile $file.fa $outputdir $basename.profile
	fi
done
cat $outputdir/consensus_clusters.*.profile.fa > $outputdir/consensus_clusters.fa
remove_smalls $outputdir/consensus_clusters.fa $outputdir/$prefix.fa 200

