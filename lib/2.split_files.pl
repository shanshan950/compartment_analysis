#!/usr/bin/perl
use strict;
my ($file) = @ARGV;
open(INFILE, $file) or die "Can't open file";
while(my $line = <INFILE>){
	chomp $line;
	my($bin1,$bin2,$count)=split "\t", $line;
	my $chr_one= ((split /_/,$bin1)[0]);
	close OUTFILE;
	my $new_file=$chr_one;
	$new_file.=".matrix";
	open(OUTFILE,">>$new_file") or die "Can't open: $new_file $!";
	print OUTFILE "$line\n";
	}
close OUTFILE;
close INFILE;
exit;

