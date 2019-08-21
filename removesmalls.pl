## removesmalls.pl
#!/usr/bin/perl
use strict;
use warnings;

my $ifile; my $ofile;
my $cutoff=$ARGV[2];
my $flag = 0;
my $seq; my $name; my $read;

open($ifile, $ARGV[0] );
open($ofile, ">$ARGV[1]");

while( $read = <$ifile> ){
        chomp $read;
	if( $read =~ />/){
		if($flag == 0){
			$flag = 1;
			$name = $read;
			$seq = "";
		} elsif ( length($seq) > $cutoff ){
			print $ofile $name."\n".$seq."\n";
			$name = $read;
			$seq = "";
		} else{
                        $name = $read;
                        $seq = "";
		}
	} else{
		$seq .= $read;
	}
}
if ( length($seq) > $cutoff ){
	print $ofile $name."\n".$seq."\n";}

