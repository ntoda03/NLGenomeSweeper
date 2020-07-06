#!/bin/bash

# Get directory of the script
programdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $programdir/../scripts/functions.sh
source $programdir/../scripts/interproscan_functions.sh

## Check dependencies
echo  "Checking software requirements..."  | tee -a $sdout
DEPENDENCIES=(samtools bedtools blastp blastx TransDecoder.LongOrfs TransDecoder.Predict muscle interproscan hmmbuild)
dep_flag=1
for dependency in "${DEPENDENCIES[@]}"
do
    if ! which $dependency &> /dev/null; then
        echo "The dependency" $dependency " could not be found. Ensure that is installed and in your path." 
        dep_flag=0
    fi
done
if [ $dep_flag -eq 1 ]; then
    echo -e "All required software was found.\n"
else
    echo "Some dependencies not found. Exiting."
    exit
fi

#echo "################### Running unit tests ###################"

outputdir=$programdir/testing_dir/01_candidate_identification/genome_search
outname="All_candidates"
genome=$programdir/data/TAIR10_chr1.sample.fa
rm -R $outputdir 2> /dev/null
mkdir -p $outputdir
cp $programdir/data/01_candidate_identification/genome_search/$outname.pos.txt $outputdir/

test_extract_seq(){
    extract_seq $outputdir/$outname.pos.txt $genome $outputdir/$outname.fa
    assertTrue "Error: extract_seq failed" "[ -r '$outputdir/$outname.fa' ]"
    test1=$(cat $outputdir/$outname.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/$outname.fa)
    assertEquals "Warning: extract_seq output no what expected" "$test1" "$test2" 
}

cp $programdir/data/01_candidate_identification/genome_search/candidate_sites.bed $outputdir/
test_get_flanking_regions(){
    get_flanking_regions $outputdir/candidate_sites.bed $genome $outputdir/Candidate_sites
    assertTrue "Error: get_flanking_regions failed" "[ -r '$outputdir/Candidate_sites.with_flanking.fa' ]"
    assertTrue "Error: get_flanking_regions failed" "[ -r '$outputdir/Candidate_sites.wflanking.txt' ]"
    test1=$(cat $outputdir/Candidate_sites.with_flanking.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/Candidate_sites.with_flanking.fa)
    assertEquals "Warning: get_flanking_regions output no what expected" "$test1" "$test2" 
    test1=$(cat $outputdir/Candidate_sites.wflanking.txt)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/Candidate_sites.wflanking.txt)
    assertEquals "Warning: get_flanking_regions output no what expected" "$test1" "$test2" 
}

cp $programdir/data/01_candidate_identification/genome_search/genome_blast.merge.txt $outputdir/
test_merge_exons(){
    merge_exons $outputdir/genome_blast.merge.txt $outputdir/genome_blast.merge2.txt 1000
    assertTrue "Error: merge_exons failed" "[ -r '$outputdir/genome_blast.merge2.txt' ]"
    test1=$(cat $outputdir/genome_blast.merge2.txt)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/genome_blast.merge2.txt)
    assertEquals "Warning: merge_exons output no what expected" "$test1" "$test2" 
}

outputdir2=$programdir/testing_dir/02_domain_identification
rm -R $outputdir2 2> /dev/null
mkdir -p $outputdir2
cp $programdir/data/02_domain_identification/Candidate_sites_interpro.*.gff3 $outputdir2/
test_convert_interpro_names(){
    convert_interpro_names $outputdir2 8
    for i in 0 1 2 3 4 5 6 7 8; 
    do
        assertTrue "Error: convert_interpro_names failed" "[ -r '$outputdir2/annotation.$i.gff3' ]"
        test1=$(cat $outputdir2/annotation.$i.gff3)
        test2=$(cat $programdir/data/02_domain_identification/annotation.$i.gff3)
        assertEquals "Warning: convert_interpro_names output no what expected" "$test1" "$test2" 
    done
}

cp $programdir/data/02_domain_identification/All_candidates.bed $outputdir2/
test_classify_interpro_structure(){
    classify_interpro_structure $outputdir2
    assertTrue "Error: classify_interpro_structure failed" "[ -r '$outputdir2/All_candidates.bed' ]"
    test1=$(cat $outputdir2/All_candidates.classified.bed)
    test2=$(cat $programdir/data/02_domain_identification/All_candidates.classified.bed)
    assertEquals "Warning: classify_interpro_structure output no what expected" "$test1" "$test2" 
}

#echo "################### Running Integration tests ###################"

outputdir_firstpass=$programdir/testing_dir/01_candidate_identification/genome_search/first_pass
rm -R $outputdir_firstpass 2> /dev/null
mkdir -p $outputdir_firstpass
cp $genome $outputdir_firstpass/genome.fa
cp $programdir/data/01_candidate_identification/genome_search/pfam_NB-ARC.fa $outputdir_firstpass/
makeblastdb -in $outputdir_firstpass/genome.fa -out $outputdir_firstpass/genome.fa -dbtype nucl  > /dev/null
makeblastdb -in $outputdir_firstpass/pfam_NB-ARC.fa -out $outputdir_firstpass/pfam_NB-ARC.fa \
    -dbtype prot > /dev/null
test_find_genome_candidates(){
    $programdir/../scripts/find_genome_candidates.sh $outputdir_firstpass/genome.fa $outputdir_firstpass 12 prot \
        $outputdir_firstpass/pfam_NB-ARC.fa 1000 $programdir/../ > /dev/null
    assertTrue "Error: find_genome_candidates failed" "[ -r '$outputdir_firstpass/$outname.fa' ]"
    test1=$(cat $outputdir_firstpass/$outname.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/first_pass/$outname.fa)
    assertEquals "Warning: find_genome_candidates output not what expected" "$test1" "$test2" 
}


outputdir_profile=$programdir/testing_dir/01_candidate_identification/genome_search/profiler
rm -R $outputdir_profile 2> /dev/null
mkdir -p $outputdir_profile
cp $programdir/data/01_candidate_identification/genome_search/first_pass/$outname.fa $outputdir_profile/
cp $programdir/data/01_candidate_identification/genome_search/pfam_NB-ARC.fa $outputdir_profile/
test_createProfile(){
    $programdir/../scripts/createProfile.sh $outputdir_profile/$outname.fa $outputdir_profile $programdir/../ \
        species_specific_domains $outputdir_profile/pfam_NB-ARC.fa nucl > /dev/null 2> /dev/null
    assertTrue "Error: find_genome_candidates failed" "[ -r '$outputdir_profile/species_specific_domains.fa' ]"
    test1=$(cat $outputdir_profile/species_specific_domains.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/profiler/species_specific_domains.fa)
    assertEquals "Warning: find_genome_candidates output no what expected" "$test1" "$test2" 
}

#$programdir/../scripts/find_genome_candidates.sh $outputdir/genome.fa $outputdir/second_pass $cores prot \
#    $outputdir/second_pass/pfam_NB-ARC.fa 1000 $programdir >> $sdout 2>> $errout
#$outputdir/second_pass/$outname.bed


#echo "################### Running complete Program ###################"

#outputdir=$programdir/testing_dir/NLGenomeSweeper/01_candidate_identification/genome_search
. shunit2
