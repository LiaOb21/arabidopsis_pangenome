#!/usr/bin/perl
use List::Util qw/shuffle/;
use strict;

use Getopt::Long;
my $file;
my $iteration;
my $verbose;
my $output;
my $sample;
my $mode;
GetOptions ("file=s" => \$file,
			"iterations=i" => \$iteration, 
			"campionamenti=s" =>\$sample,
			"mode=s" =>\$mode, 
            "output=s"   => \$output)   
or die("Error in command line arguments: lo script si lancia con perl shuffle_count --file nome_file --iterations 10 --sample 5-10-15.. --classificazione 93-75-23-1 --mode REI \n");



open(STAT, '>', "$output");

#apre il file con la matrice
my @AoA=&load_file_as_AoA("$file");

#matrice
my @accessions=(1..scalar@AoA-1);
my @num=split("-", $sample);

for (my $k=0; $k<scalar@num; $k++){
	#definisco alcune variabili;
	my $a=0;my @conteggio=(); my @classification=(); my @av=(); my @stv=();
	
	while ($a<$iteration) {
		$a++;
		my @acc=&shuffled_rei(\@accessions, $num[$k]) if $mode eq 'REI';
		my @acc=&shuffled(\@accessions, $num[$k]) if $mode eq 'NOT_REI';
		@accessions=@{$acc[0]};
		my @shuffled=@{$acc[1]};

#stampo il sunset della matrice in un file tmp
&print_shuffled(\@AoA, \@shuffled);

#apro il file con i dati ricampionati
	my @AoA_sh=&load_file_as_AoA("tmp.txt");

#conto il numero totale di geni	
	my $conteggio=&count_genes(\@AoA_sh); 
	push(@conteggio, $conteggio);
#classifico i geni 
 my @ty=&classify_genes(\@AoA_sh, $num[$k]);
	push(@classification, [@ty]);
}
#print scalar@classification;
#vado alle stampe
my @classification2=&transpose_matrix(\@classification);
for (my $i=0; $i<scalar@classification2; $i++){
	$av[$i]=&average(\@{$classification2[$i]});
	$stv[$i]=&stdev(\@{$classification2[$i]});
	}
my $av=&average(\@conteggio);unshift (@av, $av);
my $stv=&stdev(\@conteggio); unshift (@stv, $stv);
	
	
	
print STAT $num[$k]."\t";
for (my $l=0; $l<scalar@av; $l++){
print STAT $av[$l]."\t".$stv[$l]."\t";}
print STAT "\n";
}

sub average{
        my @self = @_;
        my @data=@{$self[0]};
       
        my $total = 0;
        for (my $i=0; $i<scalar@data; $i++) {
                $total += $data[$i];
        }
        my $average = $total / scalar@data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub shuffled{
my @self=@_;
my @accessions=@{$self[0]};
my $num=$self[1];
my $count=0;
while ($count<50){
@accessions=shuffle(@accessions);
$count++;
}
my @shuffled=@accessions[0..$num-1];
my @return=([@accessions], [@shuffled]);
return @return;
}

sub make_array_from_scalar {
	my @self=@_;
	my $scalar=$self[0];
	my $separator=$self[1] || "\t";
	chomp $scalar;
	my @array=split("$separator", "$scalar");
	return @array;
						}
											
sub load_file_as_array {
		my @self=@_;
		my $name_file=$self[0];
open(FILE, '<', "$name_file");
		my @array=(<FILE>);
close FILE;

		return @array;
					}
									
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
					
sub print_shuffled {

my @self=@_;
my @data=@{$self[0]};
my @shuffled=@{$self[1]};

open(output, ">", "tmp.txt");
for(my $i=0; $i<scalar@shuffled; $i++){
	my @temp=@{$data[$shuffled[$i]]};
	my $stampo=join("\t", @temp);
	print output $stampo."\n";
}
close output;
}

sub count_genes{
	my @self=@_;
	my @selected=@{$self[0]};
	my $COUNT=0;
	for (my $i=1; $i<scalar@{$selected[0]}; $i++){
	my $count=0;
		for (my $j=0; $j<scalar@selected; $j++){
			
			if ($selected[$j][$i] >0){
				
				$COUNT++;
				last}
			}}
return $COUNT;			
	}

sub shuffled_rei{
my @self=@_;
my @accessions=@{$self[0]};
my $num=$self[1];
my $count=0;
my @shuffled=();

while ($count<$num){
@accessions=shuffle(@accessions);
my $el=$accessions[1];
push (@shuffled, $el);
$count++;
}

my @return=([@accessions], [@shuffled]);
return @return;
}

sub classify_genes {
	my @self=@_;
	my @array=@{$self[0]};
	my $max=$self[1];
	my @thresholds=($max+1, $max, int($max*0.8), 2, 1);
	my @count=();
	my @return=();
	
		for (my $i=1; $i<scalar@{$array[0]}; $i++){
			for (my $j=0; $j<scalar@array; $j++){
				$count[$i]=$count[$i]+1 if $array[$j][$i] >0 ;
				$count[$i]=$count[$i]+0 if $array[$j][$i] <1 ;  
								}
							}
				
		
		for(my $i=1; $i<scalar@thresholds; $i++){
			my $con=0;
			for(my $j=0; $j<scalar@count; $j++){
				
				$con++ if (($count[$j]< $thresholds[$i-1]) && ($count[$j] >=$thresholds[$i])); 
								}
					push(@return, $con);							}
							
		return @return;
	}

sub transpose_matrix{
	my @self=@_;
	my @AoA=@{$self[0]};
	my @transposed=();
	for(my $j=0; $j<scalar@{$AoA[0]}; $j++){	
		for (my $i=0; $i<scalar@AoA; $i++){
			$transposed[$j][$i]=$AoA[$i][$j];
										}
											}
									
	return @transposed;
					}		
