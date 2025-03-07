#!/usr/bin/perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

# Input: List of contig names
my @contigs = <STDIN>;

# Make hash of contig list
my %hash;
foreach my $element (@contigs) {
    chomp $element;
    $hash{$element} = 1;
}

# Load synteny file (e.g. proteins_raw_gtdbtax_cog_synteny.tsv.gz)
my $filename = $ARGV[0];
# Unzip file
my $fh = IO::Uncompress::Gunzip->new($filename)
    or die "Cannot open '$filename': $GunzipError\n";

# Search column
my $col_num = $ARGV[1] - 1;

# Read the file line by line
while (my $line = <$fh>) {
    my @array = split /\t/, $line;
    # Print hit line
    if (exists $hash{$array[$col_num]} && defined $hash{$array[$col_num]} && $hash{$array[$col_num]} > 0) {
        print $line;
    }
}

# Close file hundle
$fh->close();
