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
#set -o history -o histexpand
#set -x
set -e

sdout=$outputdir/NLGenomeSweeper.out
errout=$outputdir/NLGenomeSweeper.err

function exit_error {
        echo "Error: There was a problem. Check "$outputdir"/NLGenomeSweeper.err for information"
        exit 1
}


#trap 'rv=\$?; echo "Error \$rv" >> $sdout' ERR
#trap 'rv=\$?; echo "Exit \$rv" >> $sdout'  EXIT

#
# Check the current amount of free memory in MB.
#
# Input: None
#
# Output: The amount of free memory in MB, returned to stdout
#
# DO NOT PRINT ANYTHING TO STDOUT! Value is returned via stdout.
#
function check_mem_usage {
    local free_mem_MB=$(free --mega | awk 'NR == 2 {print $4}')
    echo $free_mem_MB
}

#
# Handle a request to submit an interproscan job.
#
# Job will be run in the background. The job will only be run if there is no
# existing output file and the input fasta file exists. THe job is only run
# if there is at least 3G of RAM. This does not checks limits on threads
# number so that should be checked before submitting.
#
# Input: 
# $1: output directory
# $2: The chunk number for the current fasta file being submitted
#
# Output: The ID of the Interproscan job being run, returned to stdout
#
# DO NOT PRINT ANYTHING TO STDOUT! Value is returned via stdout.
# Log info is printed to STDERR
#
function submit_interproscan_job {
    local outputdir=$1
    local i=$2

    # The return Interproscan job ID
    local running_job=0

    # File hasn't been run yet
    if [ ! -f $outputdir/Candidate_sites_interpro.$i.gff3 ]; then
        # The input file exists
        if [ -f $outputdir/Candidate_sites.with_flanking.fa_chunk_0000${i} ]; then
            # Only submit new jobs if there is at meast 3G of ram free, otherwise wait
            free_mem_MB=$(check_mem_usage)
            if [ $free_mem_MB -lt 3000 ]; then
                echo "High memory usage. Waiting..." > >(tee -a $sdout >&2) 
                while [ $free_mem_MB -lt 3000 ]; do
                    sleep 10
                    free_mem_MB=$(check_mem_usage)
                done
            fi
            # Checks are passed, submit interproscan job
            echo -e "Launching $i/$file_num" > >(tee -a $sdout >&2) 
            interproscan -appli PANTHER,Gene3D,Pfam,SMART,Coils -i $outputdir/Candidate_sites.with_flanking.fa_chunk_0000${i} \
                -t n -o $outputdir/Candidate_sites_interpro.$i.gff3 -f GFF3 > $sdout.interproscan_$i 2> $errout.interproscan_$i &
            running_job=$!
            sleep 20
        fi
    fi

    echo $running_job
}

#
# Handle job being run in the background.
#
# This handles a job that have been submitted to run in the background.
# If the job has finished or run for longer than 6400 a 1 is returned to
# indicate that it should be removed from the list of active jobs. If it
# Has been running for too long it will be killed.
#
# Input: 
# $1: PC ID of the background job being run
#
# Output: A flag indicating if a process has finished/been killed (1) or not (0), returned to stdout
#
# DO NOT PRINT ANYTHING TO STDOUT! Value is returned via stdout.
# Log info is printed to STDERR
#
function handle_bg_jobs {
    local pc_id=$1

    local remove_process=0

    curr_status=$(ps aux | awk -v pc_id=$pc_id '$2==pc_id {print $11}')
    curr_runtime=$(ps -o etime= -p $pc_id | tr '-' ':' | \
        awk -F: '{ total=0; m=1; } { for (i=0; i < NF; i++) {total += $(NF-i)*m; m *= i >= 2 ? 24 : 60 }} {print total}')

    # Process is no longer running, remove from list of current jobs
    if [ -z "$curr_status" ]; then
        echo -e "Job complete" > >(tee -a $sdout >&2) 
        remove_process=1
    # Process has been running for too long, possible problem with interproscan, kill job
    elif [ "$curr_runtime" -gt 6400 ]; then
        kill $pc_id
        echo -e "Process with status $curr_status killed for hanging" > >(tee -a $sdout >&2) 
        remove_process=1
    fi

    echo $remove_process
}


#
# Verify whether interproscan ran succesfully.
#
# Check whether an output gff file was created for every input fasta file
#  to be run with Interproscan.
#
# Input: 
# $1: output directory
# $2: total number of sequences
#
# Output: A flag indicating whether all files ran successfully (1) or not (0), returned to stdout
#
# DO NOT PRINT ANYTHING TO STDOUT! Value is returned via stdout.
# Log info is printed to STDERR
#
function verify_interproscan_run {
    local outputdir=$1
    local seq_num=$2

    local flag=1

    # Check every sequence chunk to see whether it ran
    for i in $(eval echo "{0..$seq_num}");
    do
        # Output does not exist for chunk
        if [ ! -f $outputdir/Candidate_sites_interpro.$i.gff3 ]; then
            # Input exists for chunk
            if [ -f $outputdir/Candidate_sites.with_flanking.fa_chunk_0000${i} ]; then
                # Did ot run successfully
                flag=0
                echo -e "Problem: Candidate_sites.with_flanking.fa_chunk_0000${i} did not run." > >(tee -a $sdout >&2) 
            fi
        fi
    done
    if [ "$flag" -eq 1 ]; then
        echo -e "Check complete. No problems." > >(tee -a $sdout >&2) 
    else
        echo -e "Some sequences did not execute. Retrying..." > >(tee -a $sdout >&2) 
    fi

    echo $flag
}


#
# Run interproscan jobs.
#
# This handles the running and maintenance of Interproscan jobs.
# Number of jobs running at a time is limited by threads and
# each job must have 3G RAM free before launching.
#
# Input: 
# $1: output directory
# $2: total number of sequences
#
# Output: None
#
function run_interproscan {
    local outputdir=$1
    local seq_num=$2
    local cores=$3

    flag=0
    run_cycle=0
    file_num=$seq_num
    declare -a bg_processes=()

    # Potentially run multiple times if some Interproscan jobs failed
    while [[ "$run_cycle" -lt 10 && "$flag" -eq 0 ]]
    do
        ## Interproscan is run in the background, number of cores limits number of instances running
        ## Although it has a multithread option, that does not take advantage of the cpus
        ## Memory usage is really the limiting factor.
        for i in $(eval echo "{0..$seq_num}");
        do
            # Try to submit job (function checks memory usage)
            running_job=$(submit_interproscan_job $outputdir $i)
            if [ "$running_job" -ne 0 ]; then
                bg_processes+=($running_job)
            fi
            # Handle running background processes, wait if reached max number of threads
            while [ "${#bg_processes[@]}" -eq "$cores" ]; do
                echo "Reached max number of jobs. Waiting..." | tee -a $sdout
                index_run=0
                for pc_id in "${bg_processes[@]}"; do
                    remove_process=$(handle_bg_jobs $pc_id)
                    if [ $remove_process -eq 1 ]; then
                        unset bg_processes[$index_run]
                        bg_processes=( ${bg_processes[@]} )
                    fi
                    index_run=$((index_run+1))
                done
                if [ "${#bg_processes[@]}" -eq "$cores" ]; then
                    sleep 20
                fi
            done
        done
        echo "All jobs submitted. Waiting..." | tee -a $sdout
        wait

        ## Check to make sure interproscan run on all files, this can be a problem sometimes
        echo "Checking output files." | tee -a $sdout
        flag=$(verify_interproscan_run $outputdir $seq_num)
        run_cycle=$((run_cycle+1))
    done
}


# Convert gff files and filter for presence of Leucine-rich_repeats
function convert_interpro_names {
    local outputdir=$1
    local seq_num=$2

    for i in $(eval echo "{0..$seq_num}"); do
        # Combine the files, renaming domains of interest to their common names
        cat $outputdir/Candidate_sites_interpro.$i.gff3 | sed 's/ /__/g' | \
            awk 'NF>3 {
                split($1,a,":");split(a[2],b,"-");$1=a[1];
                if ($3=="ORF"){start=$4;stop=$5;strand=$7;$4=$4+b[1];$5=$5+b[1];} 
                else if( $3=="protein_match" || $3=="polypeptide" ){ 
                    if(strand=="+"){$4=$4*3+b[1]+start;$5=$5*3+b[1]+start;}
                    else if(strand=="-"){$7=strand;tmp1=b[1]+stop-$5*3+1;tmp2=b[1]+stop-$4*3+1;$4=tmp1;$5=tmp2}
                }
                if($3!="nucleic_acid"){print $0} }' | sort -k1,1 -k 4,4 -k 5,5 |uniq |sed 's/ /\t/g' | \
            sed 's/[^ \t]*Name=\|[^ \t]*ID=/Name=/g' | sed 's/;.*//g' | uniq  | \
            sed 's/PF07725\|SM00367\|PF13855\|SM00369\|G3DSA:3.80.10.10\|Leucine-rich_repeats,_typical_\(most_populated\)_subfamily/LRR/g' | \
            sed 's/PTHR23155[^ \t]*\|PTHR11017[^ \t]*/LRR_protein/g' | \
            sed 's/SM00255\|PF01582\|TIR_domain\|TIR_domain_profile\.\|Toll_-_interleukin_1_-_resistance\|G3DSA:3.40.50.10140\|SSF52200|/TIR/g' | \
            sed 's/PTHR36766\|PF05659/RPW8/g' | sed 's/PF00931/NB-ARC/g' | uniq | sort -k 1,1 -k 6,6 -k 7,7 | \
            awk '{printf "%s\t%s\tDomain\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$4,$5,$6,$7,$8,$9}' > $outputdir/annotation.$i.gff3
    done
}


# Try to classify the structure of the NBS-LRR gens based on upstream and downstream domains
function classify_interpro_structure {
    local outputdir=$1

    >$outputdir/to_delete.txt
    >$outputdir/All_candidates.classified.bed

    # Get a potential structure classification for candidates
    i=0
    while read chr start stop
    do
        # Check that candidates have an LRR
        lrr_count=$(grep "Name=LRR$" $outputdir/annotation.$i.gff3 | wc -l)

        # Try to determine what strand domain is on to look for flanking domains
        strand=$(awk -v start=$((start-100)) -v stop=$((stop+100)) \
            '(($5 >= start && $5 <= stop) || ($4 >= start && $4 <= stop)) && $9 == "Name=NB-ARC" {print $7}' $outputdir/annotation.$i.gff3 \
            |uniq |head -n 1)
        if [ -z "$strand" ]; then strand="+"; fi
        awk '$2!="getorf" {print $0}' $outputdir/annotation.$i.gff3 |sort -k5,5 > $outputdir/processing.tmp
        sed -i '/LRR_protein/d' $outputdir/processing.tmp
        if [ $strand == "-" ]; then
            awk -v stop=$stop '$4>stop {print $0}' $outputdir/processing.tmp > $outputdir/before.tmp
            awk -v start=$start '$5<start {print $0}' $outputdir/processing.tmp > $outputdir/after.tmp
            before=$(head -n 1 $outputdir/before.tmp |awk -v strand=$strand '{if(strand!=$7){$9="Name=None"} print $9}' |sed 's/.*Name=//g')
            after=$(tail -n 2 $outputdir/after.tmp |awk -v strand=$strand '{if(strand==$7 && $9=="Name=LRR"){print $9}}' |sed 's/.*Name=//g' |uniq)
        else
            awk -v start=$start '$5<start {print $0}' $outputdir/processing.tmp > $outputdir/before.tmp
            awk -v stop=$stop '$4>stop {print $0}' $outputdir/processing.tmp > $outputdir/after.tmp
            before=$(tail -n 1 $outputdir/before.tmp |awk -v strand=$strand '{if(strand!=$7){$9="Name=None"} print $9}' |sed 's/.*Name=//g')
            after=$(head -n 2 $outputdir/after.tmp |awk -v strand=$strand '{if(strand==$7 && $9=="Name=LRR"){print $9}}' |sed 's/.*Name=//g' |uniq)        
        fi

        # Start classification, all have N and add structure based on flanking interproscan domains
        class="N"
        if [ ! -z "$before" ]; then
            if [ "$before" == "TIR" ]; then class="T$class"; fi
            if [ "$before" == "Coil" ]; then class="C$class"; fi
            if [ "$before" == "RPW8" ]; then class="R$class"; fi
        fi
        if [ ! -z "$after" ]; then
            if [ "$after" == "LRR" ]; then class+="L"; fi
        fi

        # Candidates don't have any LRRs downstream, mark to remove
        if [ "$lrr_count" -eq 0 ]; then
            class+=";Filtered"
            echo -e "$chr\t$start\t$stop\t$class" >> $outputdir/to_delete.txt
        fi
        echo -e "$chr\t$start\t$stop\t$class" >> $outputdir/All_candidates.classified.bed
        i=$((i+1))
    done < $outputdir/All_candidates.bed
}


####################################################################################
# 
# Main
#
####################################################################################


# Interproscan needs at least 3MB RAM to run
free_mem_MB=$(check_mem_usage)
if [ $free_mem_MB -lt 3000 ]; then
    echo "Not enough free memory to run interproscan. 3G mimimum is required." 
    exit 1
fi

## Split the file because interproscan has problems with too many sequences
seq_num=$(grep ">" $outputdir/Candidate_sites.with_flanking.fa |wc -l)
seq_num=$(($seq_num-1))
split_fasta $outputdir/Candidate_sites.with_flanking.fa $outputdir/Candidate_sites.with_flanking.fa_chunk_0000

echo -e "\nRunning interproscan domain and ORF identification..."  | tee -a $sdout
echo -e "This may take several hours."  | tee -a $sdout
#echo -e "Current status can be found in $outputdir/NLGenomeSweeper.out"  | tee -a $sdout
run_interproscan $outputdir $seq_num $cores

# Interproscan can hang or crash randomly, check taht it was successful
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

