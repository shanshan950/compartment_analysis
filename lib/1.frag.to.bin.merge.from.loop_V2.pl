#!/usr/bin/perl
#my usage = "Usage: ./new_test.pl <frag_to_bin_file> <loop_file> \n";
use strict;
use warnings;
my ($frag_to_bin_file,$loop_file) = @ARGV;
my $frag_hash;
open(IN, $frag_to_bin_file);
while(my $line=  <IN>){
	chomp $line;
	my ($chr,$frag,$bin) = split "\t", $line;
	my $str= join ("_",$chr,$bin);
	$frag_hash->{$frag} = $str;
}
close(IN);

open(IN, $loop_file) || die("Error: Cannot open file $loop_file!\n");
my $data_hash;
while(my $loop_line = <IN>){
	chomp $loop_line;
        my ($frag1,$frag2,$count) = split "\t", $loop_line;
	my $binID1 = $frag_hash->{$frag1};
	my $binID2 = $frag_hash->{$frag2};
	if(not defined $data_hash->{$binID1}){
		foreach my $id1 (keys %{$data_hash}){
			foreach my $id2 (keys %{$data_hash->{$id1}}){
				print join("\t", $id1, $id2, $data_hash->{$id1}->{$id2})."\n";
			}
		delete $data_hash->{$id1};
		}
	}
	$data_hash->{$binID1}->{$binID2} += $count;
}
close(IN);
exit;
