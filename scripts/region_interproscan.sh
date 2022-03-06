#!/bin/bash

#
# Predict ORFs and domains of a fasta nucleotide sequence file.
#
# Interproscan is run on sequences containing NBS-LRR domains and
# flanking sequences.
#
# Inputs:
# $1: Output directory
# $2: Number of cores to use
# $3: Directory of the program
#
# Outputs:
# $outputdir/../All_candidates.gff3: Annotation file of domains for potential NBS-LRR genes
#

# Inputs
outputdir=$1
cores=$2
programdir=$3

source $programdir/scripts/functions.sh
source $programdir/scripts/interproscan_functions.sh
set -e

sdout=$outputdir/NLGenomeSweeper.out
errout=$outputdir/NLGenomeSweeper.err

function exit_error {
        echo "Error: There was a problem. Check "$outputdir"/NLGenomeSweeper.err for information"
        exit 1
}

# Interproscan needs at least 3MB RAM to run
free_mem_MB=$(check_mem_usage)
if [ $free_mem_MB -lt 3000 ]; then
    echo "Not enough free memory to run interproscan. 3G mimimum is required."
    echo "Only $free_mem_MB MB of free memery available." 
    exit 1
fi

## Split the file because interproscan has problems with too many sequences
seq_num=$(grep ">" $outputdir/Candidate_sites.with_flanking.fa |wc -l)
seq_num=$(($seq_num-1))

if [ "$seq_num" -eq "-1" ]; then
	echo "No candidates found." | tee -a $sdout
    touch $outputdir/../All_candidates.gff3 $outputdir/../Final_candidates.bed $outputdir/../Filtered_candidates.bed
	exit
fi

split_fasta $outputdir/Candidate_sites.with_flanking.fa $outputdir/Candidate_sites.with_flanking.fa_chunk_0000

echo -e "\nRunning interproscan domain and ORF identification..."  | tee -a $sdout
echo -e "This may take a while."  | tee -a $sdout
#echo -e "Current status can be found in $outputdir/NLGenomeSweeper.out"  | tee -a $sdout
run_interproscan $outputdir $seq_num $cores

# Interproscan can hang or crash randomly, check that it was successful
if [ "$flag" -eq 1 ]; then
    echo -e "\n\nDomain identification complete."  | tee -a $sdout
else
    echo -e "There was a problem with interproscan. Check $errout for more information." | tee -a $sdout
    exit 1
fi

echo "Converting names to human readable..." | tee -a $sdout
convert_interpro_names $outputdir $seq_num

echo "Classifying structure..." | tee -a $sdout
classify_interpro_structure $outputdir 

echo "Filtering for presence of LRRs..." | tee -a $sdout
# Remove candidates with no LRR
grep "Filtered" -v $outputdir/All_candidates.classified.bed > $outputdir/Final_candidates.bed
grep "Filtered" $outputdir/All_candidates.classified.bed >  $outputdir/Filtered_candidates.bed

num_found=$(wc -l < $outputdir/Final_candidates.bed)
echo -e "\nCandidate filtering for LRRs. $num_found final candidates found after filtering." | tee -a $sdout

echo "Outputing final candidates and annotations" | tee -a $sdout
cp $outputdir/Final_candidates.bed $outputdir/../Final_candidates.bed
cp $outputdir/Filtered_candidates.bed $outputdir/../Filtered_candidates.bed
cp $outputdir/All_candidates.classified.bed $outputdir/../All_candidates.bed
cat $outputdir/annotation.*.gff3 > $outputdir/All_candidates.gff3
sed -i '1s/^/##gff-version 3\n/' $outputdir/All_candidates.gff3
cp $outputdir/All_candidates.gff3 $outputdir/../All_candidates.gff3

echo "

Run complete!

These results are meant to be used with RNA-Seq data for manual annotation of the genes of interest. 

" | tee -a $sdout

