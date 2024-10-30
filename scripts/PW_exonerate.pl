#!/usr/bin/perl
use strict;
use Getopt::Long;

# Function to display help message
sub show_help {
    print "Usage: $0 --protein_fasta PROTEIN_FASTA --genome_fasta GENOME_FASTA --diamond_results DIAMOND_RESULTS --suffix SUFFIX\n";
    exit;
}

# Get command-line options
my ($protein_fasta, $genome_fasta, $diamond_results, $suffix);
GetOptions(
    'protein_fasta=s' => \$protein_fasta,
    'genome_fasta=s' => \$genome_fasta,
    'diamond_results=s' => \$diamond_results,
    'suffix=s' => \$suffix,
) or show_help();

# Check if all required options are provided
if (!defined $protein_fasta || !defined $genome_fasta || !defined $diamond_results || !defined $suffix) {
    show_help();
}

# Load protein and genome sequences into hashes
my %prot = &load_fasta_as_hash($protein_fasta);
my %gen = &load_fasta_as_hash($genome_fasta);

# Load the diamond results into an array of arrays (AoA)
my @AoA = load_file_as_AoA($diamond_results);

for (my $i = 0; $i < scalar @AoA; $i++) {
    my $protein_file = "Protein_$suffix.txt";   # Unique Protein file name
    my $genome_file = "Genome_$suffix.txt";     # Unique Genome file name
    my $output_file = "exonerate_results_$suffix.gff";  # Unique output file name

    # Write the current protein sequence to the unique file
    open(PROT, '>', $protein_file) or die "Cannot open file $protein_file: $!";
    print PROT ">" . $AoA[$i][1] . "\n" . $prot{$AoA[$i][1]};
    close PROT;

    # Write the current genome sequence to the unique file
    open(GEN, '>', $genome_file) or die "Cannot open file $genome_file: $!";
    print GEN ">" . $AoA[$i][0] . "\n" . $gen{$AoA[$i][0]};
    close GEN;

    # Run exonerate with the unique file names and output the results to a unique file
    my $cmd = "exonerate --query $protein_file --target $genome_file --model protein2genome --showquerygff no --showtargetgff yes --maxintron 5000 --showvulgar yes --ryo \"%qi\t%ti\t%pi\t%ql\t%em\t%ps\t%qab\t%qae\t%tab\t%tae\t%r\t%s\n\" >> $output_file";
    system($cmd) == 0 or die "Failed to run exonerate\n";
}

sub load_fasta_as_hash {
    my @self = @_;
    my $nome_file_fasta = $self[0];
    my $type = $self[1] || "String";

    my @array = &load_file_as_array("$nome_file_fasta");

    my %hash_fasta;
    my $seq;
    my $nome;

    for (my $i = 0; $i < scalar @array; $i++) {
        chomp $array[$i];

        if ($array[$i] =~ /^>([A-Za-z0-9#.:_-]*)/) {
            if ((length $seq) > 0) {
                if ($type eq 'String') {
                    $hash_fasta{$nome} = $seq;
                    undef $seq;
                } else {
                    my @seq_ = split("", $seq);
                    @{$hash_fasta{$nome}} = @seq_;
                    undef $seq;
                }
            }
            $nome = $1;
        } elsif ($array[$i] =~ /^[A-Za-z]/) {
            $seq = $seq . $array[$i];
        }

        if ($i == (scalar @array - 1)) {
            if ($type eq 'String') {
                $hash_fasta{$nome} = $seq;
                undef $seq;
            } else {
                my @seq_ = split("", $seq);
                @{$hash_fasta{$nome}} = @seq_;
                undef $seq;
            }
        }
    }
    return %hash_fasta;
}

sub load_file_as_AoA {
    my @self = @_;
    my $file_name = $self[0];
    my $separator = $self[1] || "\t";
    my @array = &load_file_as_array("$file_name");
    my @AoA;

    for (my $j = 0; $j < scalar @array; $j++) {
        chomp $array[$j];
        my @tmp = split("$separator", $array[$j]);
        push(@AoA, [@tmp]);
    }
    return @AoA;
}

sub load_file_as_array {
    my @self = @_;
    my $name_file = $self[0];
    open(FILE, '<', "$name_file") || die "non riesco ad aprire il file $name_file\n";
    my @array = (<FILE>);
    close FILE;

    return @array;
}

sub make_array_from_scalar {
    my @self = @_;
    my $scalar = $self[0];
    my $separator = $self[1] || "\t";
    chomp $scalar;
    my @array = split("$separator", "$scalar");
    return @array;
}
