#!/bin/bash

#
# Predict ORFs and domains of a fasta nucleotide sequence file.
#
# $1: output directory
# $2: number of cores to use
#

outputdir=$1
cores=$2

# Split the file because interproscan hangs with too many sequences
fastasplit -f $outputdir/Candidate_sites.with_flanking.fa -o $outputdir/ -c 500

# interproscan is run in the background, number of cores limits number of instances running
for i in {000..499};
do
	echo "Running $i/500"
	if [ ! -f $outputdir/Candidate_sites_interpro.$i.gff3 ]; then
		if [ -f $outputdir/Candidate_sites.with_flanking.fa_chunk_0000${i} ]; then
			interproscan -appli PANTHER,Gene3D,Pfam,SMART,Coils -i $outputdir/Candidate_sites.with_flanking.fa_chunk_0000${i} \
				-t n -o $outputdir/Candidate_sites_interpro.$i.gff3 -f GFF3 &
		fi
	fi
	if [ $(($((10#$i)) % $cores)) -eq $(($cores - 1)) ]; then
		wait
	fi
done
wait

# Combine the files, renaming domains of interest to their common names
cat $outputdir/Candidate_sites_interpro.[0-9][0-9][0-9].gff3 | sed 's/ /__/g' | \
	awk 'NF>3 {
		split($1,a,":");split(a[2],b,"-");$1=a[1];
		if ($3=="ORF"){start=$4;stop=$5;strand=$7;$4=$4+b[1];$5=$5+b[1];} 
		else if( $3=="protein_match" || $3=="polypeptide" ){ 
			if(strand=="+"){$4=$4*3+b[1]+start;$5=$5*3+b[1]+start;}
			else if(strand=="-"){$7=strand;tmp1=b[1]+stop-$5*3+1;tmp2=b[1]+stop-$4*3+1;$4=tmp1;$5=tmp2}
		}
		if($3!="nucleic_acid"){print $0} }' | sort -k1,1 -k 4,4 -k 5,5 |uniq |sed 's/ /\t/g' |
    sed 's/[^ \t]*Name=\|[^ \t]*ID=/Name=/g' | sed 's/;.*//g' | uniq  | \
    sed 's/PF13855\|SM00369\|G3DSA:3.80.10.10\|Leucine-rich_repeats,_typical_\(most_populated\)_subfamily/LRR/g' | \
    sed 's/PTHR23155[^ \t]*\|PTHR11017[^ \t]*/LRR_protein/g' | \
    sed 's/PF01582\|TIR_domain\|TIR_domain_profile\.\|Toll_-_interleukin_1_-_resistance\|G3DSA:3.40.50.10140\|SSF52200|/TIR/g' | \
    sed 's/PTHR36766\|PF05659/RPW8/g' | sed 's/PF00931/NB-ARC/g' | \
    uniq | sort -k 1,1 -k 6,6 -k 7,7 |awk '{printf "%s\t%s\tDomain\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$4,$5,$6,$7,$8,$9}' > \
    $outputdir//All_candidates_regions.annotated.gff3

sed -i '1s/^/##gff-version 3\n/' $outputdir/All_candidates_regions.annotated.gff3
