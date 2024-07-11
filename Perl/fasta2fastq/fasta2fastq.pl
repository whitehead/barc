#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

#######################################################################
#
# Copyright(c) 2014-2018 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Inmaculada Barrasa
#
# Version: 1.0
# Version: 1.1: Update to work with multi-line fasta input and print to stdout -- GB 5 Oct 2018
# 
#######################################################################

# Take a fasta file of short reads and convert it to a fastq file Sanger Qualities format. 
# Since there are no qualities use I as quality for all the bases

my %opts;
my $debug = 0;
my $fileIn;
my $ID;
my $seq = "NOTASEQ";
my $qual;

# GB - 12/6/16
my $dummyQuality = "I";

getopts('f:d:l:',\%opts);
if ( $opts{'f'} ) {
	$fileIn= $opts{'f'};
}

my $info = "\nConvert a fasta file into fastq, 
using '$dummyQuality' as the quality for all bases.
\nUSAGE: $0 -f fastaFile.fa > fastqFile.fq\n\n";

if (!$fileIn) {
	print "$info";
	exit;
}

open (FH, "$fileIn") || die "can not open $fileIn\n";

while (<FH>) {
	
	chomp;
	s/^\s+//;
	my $lane = $_;
	if ($lane =~ /\>/) {
		$lane =~ s/\>//;
		#PRINT prev seq is differnet than"-":
		if (!($seq =~ /NOTASEQ/)) {
		#make a fake quality controls string
			$qual = $seq;
			$qual =~ s/\w/$dummyQuality/g;
			print "@".$ID."\n".$seq."\n+\n".$qual."\n";
		}
		$ID = $lane;
		$seq = "";
	
	} else {
		$seq = $seq . $_;
		$qual = $seq;
		$qual =~ s/\w/$dummyQuality/g;
	}
}

print "@".$ID."\n".$seq."\n+\n".$qual."\n";

print STDERR "\nAll done!\n\n";
