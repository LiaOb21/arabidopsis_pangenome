#!/usr/bin/perl
use strict;

my %hash=&load_table_as_hash("$ARGV[0]", "1");


open(INPUT, '<', "$ARGV[1]");
my $linea=<INPUT>;
my @array=split("\t", $linea); 
shift(@array); 

while ($linea=<INPUT>){
	chomp $linea; 
		my $path_name=&take_name($linea);
		my @nodes=&take_nodes($linea);
		my @counts=&count_node_types(\@nodes, \%hash, \@array);
		my $stampo=join("\t", @counts);
		print $path_name."\t".$stampo."\n";
		&print_R_ready("data_set_R.txt", $path_name, \@counts);
				}
						

sub take_name{
my @self=@_;
my @tmp=split("\t", $self[0]);
return $tmp[0];
			}
			
sub take_nodes{
my @self=@_;
my @tmp=split("\t", $self[0]);
shift(@tmp);
return @tmp;
			}
			
sub count_node_types{
my @self=@_;
my @nodes=@{$self[0]};
my %hash=%{$self[1]};
my @loci=@{$self[2]};

my %count=();
$count{'core'}=0;$count{'softcore'}=0;$count{'dispensable'}=0;$count{'private'}=0;$count{'missing'}=0;
	for (my $i=0; $i<scalar@nodes; $i++){
	my @tmp=@{$hash{$loci[$i]}} if exists $hash{$loci[$i]};
	$count{$tmp[1]}++ if $nodes[$i]>0 && exists $hash{$loci[$i]};
	#$count{"missing"}++ if $nodes[$i]<1 && exists $hash{$loci[$i]};
	}

my @composizione=($count{'core'}, $count{'softcore'}, $count{'dispensable'}, $count{'private'}); #, $count{'missing'});
return @composizione;
}

sub load_table_as_hash {
	my @self=@_;
	my $file_dati=$self[0];
	my @AoA=load_file_as_AoA("$file_dati");
	my %hash=();
	for (my $j=0; $j<scalar@AoA; $j++){
		@{$hash{$AoA[$j][0]}}=@{$AoA[$j]};
	}
										
return %hash;
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
sub print_R_ready {
my @self=@_;
my $file=$self[0];
my $path_name=$self[1];
my @count=@{$self[2]};
my @return;
open(OUTPUT, '>>', $file);
my @types=("core", "softcore", "dispensable", "private", "missing");
for (my $i=0; $i<scalar@count; $i++){
push (@return, $path_name."\t".$types[$i]."\t".$count[$i]."\n");
}
print OUTPUT @return;
}

