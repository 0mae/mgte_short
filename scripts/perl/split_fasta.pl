#!/usr/bin/perl

use strict;

## USAGE: perl split_fasta.pl seqs fasta out_dir digits
# seqs: the number of sequences the split fasta have
# digits: the number of digits of sufix
# fasta: input fasta file. The split files are written into the fasta.split directory

## Feature: This doesn't need large memory unlike seqkit split command

my $seqs = 0;
my $n = 0;
my $digits = "%0" . "$ARGV[3]" . "d";
my @chunk;

my $str = $ARGV[1];
$str =~ s/.*\///;
my $prefix = "$ARGV[2]/$str";

open IN, $ARGV[1] or die "Cannot open $ARGV[1]: $!";
while(<IN>){
	if (m/^>/ && $seqs == $ARGV[0] - 1){
		$seqs++;
		my $k = sprintf($digits,$n);
		my $filename = "$prefix.$k";
		open OUT, "> $filename";
		print OUT @chunk;
		close OUT;
		$seqs = 0;
		@chunk = ();
		push @chunk, $_;
		$n++;
	} elsif (m/^>/ && $seqs < $ARGV[0] - 1){
		$seqs++;
		push @chunk, $_;
	} elsif (eof(IN)){
		push @chunk, $_;
		my $k = sprintf($digits,$n);
		my $filename = "$prefix.$k";
		open OUT, "> $filename";
		print OUT @chunk;
		close OUT;
	} else {
		push @chunk, $_;
	}
}
close IN;
