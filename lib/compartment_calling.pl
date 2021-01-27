#!/usr/bin/perl
use strict;
my $usage = "Usage: ./compartment calling.pl <ChIPSeq_file> <component_file>\n";
my ($ChIPSeq_file,$component_file) = @ARGV;
if(not defined $ChIPSeq_file){
        die($usage);
}
my $chip_hash;
open(IN, $ChIPSeq_file);
while(my $line= <IN>){
	chomp $line;
	my ($chr,$from,$to,$bin,@rest)= split "\t", $line;
	my $str = join("_",$chr,$bin);
	if($chip_hash->{$str}){
		$chip_hash->{$str} = $chip_hash->{$str} + 1;             
	}else{
		$chip_hash->{$str} = 1;
	}
}
close(IN);
open(IN, $component_file) || die("Error: Cannot open file $component_file!\n");
while(my $line = <IN>){
	my ($chr_bin,$compo1,$compo2,$compo3) = split "\t", $line;
	my $count=0;
	if($chip_hash->{$chr_bin}){
		$count = $chip_hash->{$chr_bin};
	}
	print join("\t",$chr_bin,$count,$compo1,$compo2,$compo3);	
}
close(IN);
exit;
