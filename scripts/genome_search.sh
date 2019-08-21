#!/bin/bash

#
# Search a genome file for NBS-LRR genes using consensus domain sequence(s).
#
# $1: genome file in fasta format to search
# $2: output directory
# $3: directory of this program
# $4: number of cores to use
# $5: The maximum size (bp) of introns
#

genome=$1
outputdir=$2
programdir=$3
cores=$4
intron=$5

source $programdir/scripts/functions.sh

sdout=$outputdir/NLGenomeSweeper.out
errout=$outputdir/NLGenomeSweeper.err
outname="All_candidates"

function exit_error {
	echo "Error: There was a problem. Check "$outputdir"/NLGenomeSweeper.err for information"
	exit 1
}

function format_error {
	file=$1
	format=$2

	echo "Error: The file "$file" is not in the correct "$format" format." | tee -a $errout
	exit_error
}

#{
echo "

Running NLGenomeSweeper

Any problems can be submitted via the GitHub page https://github.com/ntoda03/NLGenomeSweeper

" | tee -a $sdout


## Check dependencies
echo "Checking software requirements..." | tee -a $sdout
DEPENDENCIES=(samtools bedtools blastp blastx TransDecoder.LongOrfs TransDecoder.Predict muscle interproscan hmmbuild)
dep_flag=1
for dependency in "${DEPENDENCIES[@]}"
do
	if ! which $dependency &> /dev/null; then
    		echo "The dependency" $dependency " could not be found. Ensure that is installed and in your path." | tee -a $errout
		depflag=0
	fi
done
if [ $dep_flag -eq 0 ]; then
	exit_error
fi
echo -e "All required software was found.\n" | tee -a $sdout

cp $genome $outputdir/genome.fa
makeblastdb -in $outputdir/genome.fa -out $outputdir/genome.fa -dbtype nucl >> $sdout 2>> $errout || format_error $genome "nucleotide fasta"

#################### First pass identification of candidates based on consensus sequences #######################
echo "Beginning first pass of candidate finding..." | tee -a $sdout
mkdir -p $outputdir/first_pass
cp $outputdir/pfam_NB-ARC.fa $outputdir/first_pass/
makeblastdb -in $outputdir/first_pass/pfam_NB-ARC.fa -out $outputdir/first_pass/pfam_NB-ARC.fa \
	-dbtype prot > $sdout 2> $errout || format_error "custom consensus" "protein fasta"
$programdir/scripts/find_genome_candidates.sh $outputdir/genome.fa $outputdir/first_pass $cores prot \
	$outputdir/first_pass/pfam_NB-ARC.fa $intron $programdir >> $sdout 2>> $errout
if [ -s "$outputdir/first_pass/$outname.fa" ]
then 
   echo -e "First pass complete.\n" | tee -a $sdout
else
   exit_error
fi


######################## Use candidates to create species specific consensus sequences ##########################
echo "Creating species and class specific custom sequences..." | tee -a $sdout
mkdir -p $outputdir/profiler
$programdir/scripts/createProfile.sh $outputdir/first_pass/$outname.fa $outputdir/profiler $programdir \
	species_specific $outputdir/pfam_NB-ARC.fa nucl >> $sdout 2>> $errout
if [ -s "$outputdir/profiler/species_specific_domains.fa" ]
then
   echo -e "Custom profiles complete.\n" | tee -a $sdout
else
   exit_error
fi

### Second pass, use the species specific consensus sequences to do a second pass of candidate identification ###
echo "Beginning second pass of candidate finding..." | tee -a $sdout
mkdir -p $outputdir/second_pass
cat $outputdir/profiler/species_specific_domains.fa | sed 's/consensus_clusters.//g' |sed 's/.txt//g' > $outputdir/Species_specific_consensus.fa
cat $outputdir/first_pass/pfam_NB-ARC.fa $outputdir/profiler/species_specific_domains.fa > $outputdir/second_pass/pfam_NB-ARC.fa
$programdir/scripts/find_genome_candidates.sh $outputdir/genome.fa $outputdir/second_pass $cores prot \
	$outputdir/second_pass/pfam_NB-ARC.fa $intron $programdir >> $sdout 2>> $errout
if [ -s "$outputdir/second_pass/$outname.bed" ]
then
   echo -e "Second pass complete.\n" | tee -a $sdout
else
   exit_error
fi

################################# Combine results from first and second pass #####################################
echo "Extracting sequences..." | tee -a $sdout
cat $outputdir/first_pass/$outname.bed $outputdir/second_pass/$outname.bed | \
	bedtools sort | bedtools merge -d 10 > $outputdir/$outname.bed
#cp $outputdir/second_pass/$outname.bed |bedtools sort > $outputdir/$outname.bed
awk '{printf "%s:%s-%s\n",$1,$2,$3}' $outputdir/$outname.bed > $outputdir/$outname.pos.txt
extract_seq $outputdir/$outname.pos.txt $genome $outputdir/$outname.fa
cp $outputdir/$outname.bed $outputdir/candidate_sites.bed
# If a protein search was run then there are already some candidate sites	
#cat $outputdir/Unannotated_candidates.bed >> $outputdir/candidate_sites.bed
#cat $outputdir/candidate_sites.bed |bedtools sort |bedtools merge > $outputdir/candidate_sites.bed.tmp
#mv $outputdir/candidate_sites.bed.tmp $outputdir/candidate_sites.bed
# Extend candidate sites to include 10 kb of flanking sequence for domain search 
get_flanking_regions $outputdir/candidate_sites.bed $genome $outputdir/Candidate_sites
echo -e "Extraction complete." | tee -a $sdout

num_found=$(wc -l < $outputdir/$outname.bed)
echo -e "\nCandidate search complete. $num_found candidates found." | tee -a $sdout


