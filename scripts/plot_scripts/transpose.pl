#!/usr/bin/perl
use strict;

my @AoA=&load_file_as_AoA($ARGV[0]);

for (my $j=0; $j<scalar@{$AoA[$j]}; $j++){
	for (my $i=0; $i<scalar@AoA; $i++){
		print $AoA[$i][$j];
		print "\t" if $i <(scalar@AoA-1);
		};
		print "\n";
}


sub load_file_as_array {
		my @self=@_;
		my $name_file=$self[0]; 
open(FILE, '<', "$name_file") || die "non riesco ad aprire il file $name_file\n";;
		my @array=(<FILE>);
close FILE;

		return @array;
					}
########################################################################			
sub make_array_from_scalar {
	my @self=@_;
	my $scalar=$self[0];
	my $separator=$self[1] || "\t";
	chomp $scalar;
	my @array=split("$separator", "$scalar");
	return @array;
						}
########################################################################
sub load_file_as_AoA{
	my @self=@_;
	my $file_name=$self[0]; 
	my $separator=$self[1] || "\t";
	my @array=&load_file_as_array("$file_name");
	my @AoA;
	
	for (my $j=0; $j<scalar@array; $j++){
		chomp $array[$j];
		my @tmp=split ("$separator", $array[$j]);
		push (@AoA, [@tmp]);
							}
	return @AoA;
					}
					
