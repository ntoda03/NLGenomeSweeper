#!/bin/bash

#
# testing_harness.sh
#
# This runs some basic unit and integration tests of for the pipeline.
# Also runs the complete pipeline on test data.
#
# To run: 
# (PATH_TO_DIR/testing_harness.sh)
#
# Inputs:
# None
#
# Temporary files:
# In the current working directory a new directory testing_dir will be created for temp files.
#
# Outputs:
# Prints whether requirements are found and whether tests are passed.
#

# Get directory of the script
programdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
outdir=./
mkdir -p $outdir/testing_dir

genome=$programdir/data/TAIR10_chr1.sample.fa

source $programdir/../scripts/functions.sh
source $programdir/../scripts/interproscan_functions.sh

## Check dependencies
echo  "Checking software requirements..."  
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


test_extract_seq(){
    outputdir=$outdir/testing_dir/01_candidate_identification/genome_search
    outname="All_candidates"
    rm -R $outputdir 2> /dev/null
    mkdir -p $outputdir
    cp $programdir/data/01_candidate_identification/genome_search/$outname.pos.txt $outputdir/
    extract_seq $outputdir/$outname.pos.txt $genome $outputdir/$outname.fa
    assertTrue "Error: extract_seq failed" "[ -r '$outputdir/$outname.fa' ]"
    test1=$(cat $outputdir/$outname.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/$outname.fa)
    assertEquals "Warning: extract_seq output no what expected" "$test1" "$test2" 
    rm -R $outputdir 2> /dev/null
}

test_get_flanking_regions(){
    outputdir=$outdir/testing_dir/01_candidate_identification/genome_search
    outname="All_candidates"
    rm -R $outputdir 2> /dev/null
    mkdir -p $outputdir
    cp $programdir/data/01_candidate_identification/genome_search/candidate_sites.bed $outputdir/
    get_flanking_regions $outputdir/candidate_sites.bed $genome $outputdir/Candidate_sites
    assertTrue "Error: get_flanking_regions failed" "[ -r '$outputdir/Candidate_sites.with_flanking.fa' ]"
    assertTrue "Error: get_flanking_regions failed" "[ -r '$outputdir/Candidate_sites.wflanking.txt' ]"
    test1=$(cat $outputdir/Candidate_sites.with_flanking.fa)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/Candidate_sites.with_flanking.fa)
    assertEquals "Warning: get_flanking_regions output no what expected" "$test1" "$test2" 
    test1=$(cat $outputdir/Candidate_sites.wflanking.txt)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/Candidate_sites.wflanking.txt)
    assertEquals "Warning: get_flanking_regions output no what expected" "$test1" "$test2" 
    rm -R $outputdir 2> /dev/null
}

test_merge_exons(){
    outputdir=$outdir/testing_dir/01_candidate_identification/genome_search
    outname="All_candidates"
    rm -R $outputdir 2> /dev/null
    mkdir -p $outputdir
    cp $programdir/data/01_candidate_identification/genome_search/genome_blast.merge.txt $outputdir/
    merge_exons $outputdir/genome_blast.merge.txt $outputdir/genome_blast.merge2.txt 1000
    assertTrue "Error: merge_exons failed" "[ -r '$outputdir/genome_blast.merge2.txt' ]"
    test1=$(cat $outputdir/genome_blast.merge2.txt)
    test2=$(cat $programdir/data/01_candidate_identification/genome_search/genome_blast.merge2.txt)
    assertEquals "Warning: merge_exons output no what expected" "$test1" "$test2" 
    rm -R $outputdir 2> /dev/null
}

test_convert_interpro_names(){
    outputdir2=$outdir/testing_dir/02_domain_identification
    rm -R $outputdir2 2> /dev/null
    mkdir -p $outputdir2
    cp $programdir/data/02_domain_identification/Candidate_sites_interpro.*.gff3 $outputdir2/
    convert_interpro_names $outputdir2 8
    for i in 0 1 2 3 4 5 6 7 8; 
    do
        assertTrue "Error: convert_interpro_names failed" "[ -r '$outputdir2/annotation.$i.gff3' ]"
        test1=$(cat $outputdir2/annotation.$i.gff3)
        test2=$(cat $programdir/data/02_domain_identification/annotation.$i.gff3)
        assertEquals "Warning: convert_interpro_names output no what expected" "$test1" "$test2" 
    done
    rm -R $outputdir2 2> /dev/null
}

test_classify_interpro_structure(){
    outputdir2=$outdir/testing_dir/02_domain_identification
    rm -R $outputdir2 2> /dev/null
    mkdir -p $outputdir2
    cp $programdir/data/02_domain_identification/All_candidates.bed $outputdir2/
    cp $programdir/data/02_domain_identification/annotation.*.gff3 $outputdir2/
    classify_interpro_structure $outputdir2
    assertTrue "Error: classify_interpro_structure failed" "[ -r '$outputdir2/All_candidates.bed' ]"
    test1=$(cat $outputdir2/All_candidates.classified.bed)
    test2=$(cat $programdir/data/02_domain_identification/All_candidates.classified.bed)
    assertEquals "Warning: classify_interpro_structure output no what expected" "$test1" "$test2" 
    rm -R $outputdir2 2> /dev/null
}

#echo "################### Running Integration tests ###################"

test_find_genome_candidates(){
    outputdir_firstpass=$outdir/testing_dir/01_candidate_identification/genome_search/first_pass
    outname="All_candidates"
    rm -R $outputdir_firstpass 2> /dev/null
    mkdir -p $outputdir_firstpass
    cp $genome $outputdir_firstpass/genome.fa
    cp $programdir/data/01_candidate_identification/genome_search/pfam_NB-ARC.fa $outputdir_firstpass/
    makeblastdb -in $outputdir_firstpass/genome.fa -out $outputdir_firstpass/genome.fa -dbtype nucl  > /dev/null
    makeblastdb -in $outputdir_firstpass/pfam_NB-ARC.fa -out $outputdir_firstpass/pfam_NB-ARC.fa \
        -dbtype prot > /dev/null
    $programdir/../scripts/find_genome_candidates.sh $outputdir_firstpass/genome.fa $outputdir_firstpass 12 prot \
        $outputdir_firstpass/pfam_NB-ARC.fa 1000 $programdir/../ > /dev/null
    assertTrue "Error: find_genome_candidates failed" "[ -r '$outputdir_firstpass/$outname.fa' ]"
    #test1=$(cat $outputdir_firstpass/$outname.fa)
    #test2=$(cat $programdir/data/01_candidate_identification/genome_search/first_pass/$outname.fa)
    #assertEquals "Warning: find_genome_candidates output not what expected" "$test1" "$test2" 
    rm -R $outputdir_firstpass 2> /dev/null
}


test_createProfile(){
    outputdir_profile=$outdir/testing_dir/01_candidate_identification/genome_search/profiler
    rm -R $outputdir_profile 2> /dev/null
    mkdir -p $outputdir_profile
    cp $programdir/data/01_candidate_identification/genome_search/first_pass/$outname.fa $outputdir_profile/
    cp $programdir/data/01_candidate_identification/genome_search/pfam_NB-ARC.fa $outputdir_profile/
    $programdir/../scripts/createProfile.sh $outputdir_profile/$outname.fa $outputdir_profile $programdir/../ \
        species_specific_domains $outputdir_profile/pfam_NB-ARC.fa nucl > /dev/null 2> /dev/null
    assertTrue "Error: find_genome_candidates failed" "[ -r '$outputdir_profile/species_specific_domains.fa' ]"
    #test1=$(cat $outputdir_profile/species_specific_domains.fa)
    #test2=$(cat $programdir/data/01_candidate_identification/genome_search/profiler/species_specific_domains.fa)
    #assertEquals "Warning: find_genome_candidates output no what expected" "$test1" "$test2" 
    rm -R $outputdir_profile 2> /dev/null
}


test_region_interproscan_empty(){
    outputdir_interpro=$outdir/testing_dir/02_domain_identification
    rm -R $outputdir_interpro 2> /dev/null
    mkdir -p $outputdir_interpro
    echo "" > $outputdir_interpro/Candidate_sites.with_flanking.fa
    echo "" > $outputdir_interpro/All_candidates.bed
    $programdir/../scripts/region_interproscan.sh  $outputdir_interpro 2 $programdir/../ > /dev/null 2> /dev/null
    assertTrue "Error: region_interproscan failed" "[ -r '$outputdir_interpro/../All_candidates.gff3' ]"
    #test1=$(cat $outputdir_interpro/All_candidates.gff3)
    #test2=$(cat $programdir/data/02_domain_identification/All_candidates.gff3)
    #assertEquals "Warning: region_interproscan output no what expected" "$test1" "$test2" 
    rm -R $outputdir_interpro 2> /dev/null
}


test_region_interproscan(){
    outputdir_interpro=$outdir/testing_dir/02_domain_identification
    rm -R $outputdir_interpro 2> /dev/null
    mkdir -p $outputdir_interpro
    cp $programdir/data/02_domain_identification/All_candidates.bed $outputdir_interpro/
    cp $programdir/data/02_domain_identification/Candidate_sites.with_flanking.fa $outputdir_interpro/Candidate_sites.with_flanking.fa
    $programdir/../scripts/region_interproscan.sh  $outputdir_interpro 2 $programdir/../ > /dev/null 2> /dev/null
    assertTrue "Error: region_interproscan failed" "[ -r '$outputdir_interpro/All_candidates.gff3' ]"
    #test1=$(cat $outputdir_interpro/All_candidates.gff3)
    #test2=$(cat $programdir/data/02_domain_identification/All_candidates.gff3)
    #assertEquals "Warning: region_interproscan output no what expected" "$test1" "$test2" 
    rm -R $outputdir_interpro 2> /dev/null
}


#echo "################### Running complete Program ###################"

test_NLRGenomeSweeper(){
    outputdir_NLG=$outdir/testing_dir/NLG
    genome=$programdir/data/TAIR10_chr1.sample.fa
    rm -R $outputdir_NLG 2> /dev/null
    mkdir -p $outputdir_NLG
    $programdir/../NLGenomeSweeper -genome $genome -outdir $outputdir_NLG -t 2 > /dev/null 2> /dev/null
    assertTrue "Error: NLGenomeSweeper failed" "[ -r '$outputdir_NLG/NLGenomeSweeper/All_candidates.gff3' ]"
    #test1=$(cat $outputdir_NLG/NLGenomeSweeper/All_candidates.gff3)
    #test2=$(cat $programdir/data/02_domain_identification/All_candidates.gff3)
    #assertEquals "Warning: NLGenomeSweeper output no what expected" "$test1" "$test2" 
    rm -R $outputdir_NLG 2> /dev/null
}

test_CustomProfiler(){
    outputdir_profiler=$outdir/testing_dir/profiler
    genes=$programdir/data/TAIR_TNLs.fa
    rm -R $outputdir_profiler 2> /dev/null
    mkdir -p $outputdir_profiler
    $programdir/../CustomProfiler.py -verified $genes -prefix tair_tnls -outdir $outputdir_profiler > /dev/null 2> /dev/null
    assertTrue "Error: CustomProfiler failed" "[ -r '$outputdir_profiler/00_profile_creation/tair_tnls.fa' ]"
    test1=$(cat $outputdir_profiler/00_profile_creation/tair_tnls.fa)
    test2=$(cat $programdir/data/00_profile_creation/tair_tnls.fa)
    assertEquals "Warning: CustomProfiler output no what expected" "$test1" "$test2" 
    rm -R $outputdir_profiler 2> /dev/null
}

. shunit2
