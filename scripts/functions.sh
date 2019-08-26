#!/bin/bash

#
# Create a new gff file containing only the genes contained in the input
# file. Warning: uses grep so gene names should be formatted accordingly.
#
# $1: file containing genes to get
# $2: gff annotation file
# $3: output file in gff format
#
function generate_gff {
    genes=$1
    gff=$2
    output=$3
    
    >$output
    while read gene
    do
        grep $gene $gff | awk -v gene=$gene '{$9=gene; print $0}' >> $output
    done < $genes
}

#
# Retreive genes from a gff3 file and transform into bed format.
#
# $1: input file of gff3 annotation
# $2: output file in bed format
#
function gff3_to_bed {
    gff=$1
    output=$2
    
    awk '$3=="gene" {print $0;}' $gff | awk '{printf "%s\t%s\t%s\t%s\n",$1,$4,$5,$9}' | bedtools sort > $output
    if [ ! -s $output ]; then
        awk '$3=="mRNA" {print $0;}' $gff | awk '{printf "%s\t%s\t%s\t%s\n",$1,$4,$5,$9}' | bedtools sort > $output
    fi
}

#
# Merge positions from a file in format 'chr\tstart\tstop\tlength' where
# stop is always greater than start. Posiitons will be merged if their 
# distance is less than or equal to the max_size. Lengths will be added
# when positions are merged.
#
# $1: input file, format 'chr\tstart\tstop\tlength'
# $2: output file
# $3: maximum gap size
#
function merge_exons {
    infile=$1
    outfile=$2
    max_size=$3

    awk -v max_size=$max_size '{
       if ( NR == 1 ){
            choice=1
            last_chr=$1; last_start=$2; last_stop=$3; last_len=$4;}
       else if ( $1==last_chr && ($2-last_stop)<=max_size ){
            choice=2
            last_stop=$3; last_len += $4 }
       else { 
            choice=3
            print last_chr, last_start, last_stop, last_len;
            last_chr=$1; last_start=$2; last_stop=$3; last_len=$4;
       }
    }
    END {print last_chr, last_start, last_stop, last_len;} ' $infile | sed 's/ /\t/g' > $outfile

}


#
# Extract sequences from a fasta file based on a text file of positions of the
# format 'id' or 'id:start-stop' where id is the name of an element in the fasta file
#
# $1: file containing sequences to extract, position optional
# $2: genome file
# $3: output file
#
function extract_seq {
    extract=$1
    genome=$2
    output=$3
    
    samtools faidx $genome
    >$output
    while IFS=' ' read col1 col2
    do
        if [ ! -z ${col2} ]; then samtools faidx $genome $col1 | sed "s/>.*/>$col2/g" >> $output
        else samtools faidx $genome $col1 >> $output
        fi
    done < $extract
}

#
# Remove short sequences from a fasta file
#
# $1: fasta file
# $2: output file
# $3: minimum sequence length (inclusive)
#
function remove_smalls {
    infile=$1
    outfile=$2
    cutoff=$3

    samtools faidx $infile

    >$infile.extract_list
    while read chr length start l1 l2
    do
        if [[ "$length" -ge "$cutoff" ]]; then
            echo $chr >> $infile.extract_list
        fi
    done < $infile.fai

    extract_seq $infile.extract_list $infile $outfile
}

#
# Delete a line from a file based on a file of strings.
#
# $1: file containing strings to look for
# $2: file to delete from 
#
function delete_seq {
    substrings=$1
    file=$2
    
    while read id
    do
          sed -i '/$id/d' $file
    done < $substrings
    sed -i '/^\s*$/d' $file
}

#
# Extract blast hit sequences from a fasta file.
#
# $1: input file, blast output in format 6
# $2: sequence fasta file to extract from
# $3: output file
#
function extractBlastSeq {
    input=$1
    sequence=$2
    output=$3
    
    awk '{printf "%s\t%s\t%s\n",$2,$9,$10}' $input | sort -k1,1| uniq | \
    awk '{printf "%s:%s-%s\n", $1,$2,$3}' > $input.extract.bed
    extract_seq $input.extract.bed $sequence $output
}

#
# Extract blast hit sequences from a fasta file.
#
# $1: input file, blast output in format 6
# $2: sequence fasta file to extract from
# $3: output file
#
function extractBlastedScaff {
    input=$1
    sequence=$2
    output=$3
    
    awk '{print $1}' $input |sort -k1,1 |uniq > $input.extractScaff.bed
    extract_seq $input.extractScaff.bed $sequence $output
}

#
# Processing steps to do profile creation from hmmer output file.
#
# $1: input file, output from hmmer
# $2: proteins to search in fasta format
# $3: output directory
# $4: name prefix
#
function process_custom_profile {
    input=$1
    proteins=$2
    outputdir=$3
    name=$4

    awk '!/^#/ {printf "%s:%s-%s\n",$1,$18,$19}' $input > $input.hmmersearch.list.txt
    extract_seq $input.hmmersearch.list.txt $proteins $input.hmmersearch.fa
    create_custom_profile $input.hmmersearch.fa $outputdir $name
}

#
# Creation of hmmer output custom profile.
#
# $1: input file, output from hmmer
# $2: output directory
# $3: name prefix
# $4: format, nucl or prot
#
function create_custom_profile {
    input=$1
    outputdir=$2
    name=$3
    
    muscle -clw -in $input -out $input.muscle.clw
    hmmbuild -n $name $outputdir/$name.hmm $input.muscle.clw
    hmmemit -c $outputdir/$name.hmm > $outputdir/$name.fa
}
 
#
# Create a custom HMM profile and consensus sequence of NB-ARC domain.
#
# $1: Pfam NB-ARC HMM profile
# $2: proteins to search in fasta format
# $3: output prefix
#
#function create_custom_profile {
#    hmm=$1
#    proteins=$2
#    output=$3
#    
#    hmmsearch -E 1e-4 -o $output.hmmersearch_custom.txt $hmm $proteins
#    process_custom_profile $output.hmmersearch_custom.txt $proteins $output
#}


#
# Extend a bed file to 10 kb in either direction.
#
# $1: bed file of positions of interest
# $2: genome file
# $3: output prefix
#
function get_flanking_regions {
    sites=$1
    genome=$2
    output=$3
    
    samtools faidx $genome
    join <( sort -k1,1 $sites) <( awk '{print $1,$2}' $genome.fai |sort -k1,1) \
        | awk '{if($2-10000 < 0){$2=0} else{ $2=$2-10000};if($3+10000 >$4){$3=$4} else{ $3=$3+10000};printf "%s\t%s\t%s\n",$1,$2,$3}' \
        |bedtools sort | awk '{printf "%s:%s-%s\n", $1,$2,$3}'  > $output.wflanking.txt
    extract_seq $output.wflanking.txt $genome $output.with_flanking.fa
}


#
# Split a fasta file into multiple fasta files, one sequence per file
#
# $1: input fasta file
# $2: output prefix for split files
#
function split_fasta {
    fasta=$1
    prefix=$2

    samtools faidx $fasta
    i=0
    while read pos start stop l1 l2; do
        samtools faidx $fasta $pos >> $prefix$i
        i=$((i+1))
    done < $fasta.fai
}

#
# Split a fasta file into multiple fasta files (n files)
#
# $1: input fasta file
# $2: n, number of files to split into (may create fewer files if n > number of lines in file)
# $3: output prefix for split files
#
function split_fasta_n {
    fasta=$1
    n_split=$2
    prefix=$3
    
    samtools faidx $fasta
    awk '{print $1}' $fasta.fai > $fasta.fai.count
    total_lines=$(wc -l $fasta.fai.count |awk '{print $1}')
    ((lines_per_file = ($total_lines + $n_split - 1) / $n_split))
    split --lines=$lines_per_file $fasta.fai.count $fasta.fai_chunk_ -d -a 3
    for i in $(eval echo "{0..$n_split}"); do
        if [ -f $fasta.fai_chunk_$i ]; then
            extract_seq $fasta.fai_chunk_$i $fasta ${prefix}$i
        fi
    done
}
