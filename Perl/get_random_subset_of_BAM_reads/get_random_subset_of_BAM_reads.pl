#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2011 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George W. Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 1.0
#
###########################################################################

# Basic ideas:
#  1 - Take off header rows (and save for printing)
#  2 - Shuffle remaining row number
#  3 - Select first n rows of shuffled rows (ir n reads are requested)
#  4 - Print original header and selected rows (in same order as in original file)

# Note that each row is treated individually, so multiple-mapped reads are counted as >1 read
# !!! This code doesn't work for paired end reads.

# Sample command: 
#  ./get_subset_of_SAM_reads.pl accepted_hits.sam 100 > 100.sam

if (! $ARGV[1])
{
	print STDERR "\nGet a random subset of mapped reads in a SAM file\n";
	print STDERR "USAGE: $0 SAMfile numReadsWanted > out.sam\n\n";
	exit;
}

$SAM_file = $ARGV[0];
$numReadsWanted = $ARGV[1];

$header = "";

if ($SAM_file =~ /\.SAM/i)
{
	print STDERR "Reading SAM file to get number of non-header reads ...\n";
	open (SAM, $SAM_file) || die "Cannot open $SAM_file) for reading: $!";
}
elsif ($SAM_file =~ /\.BAM/i)
{
	print STDERR "Reading BAM file to get number of non-header reads ...\n";
	open (SAM, "samtools view -h $SAM_file |") || die "Cannot open $SAM_file) for reading: $!";
}

$lineNum = 0;
while (<SAM>)
{
	if (! /^@/)	# Save these line numbers
	{
		push @lineNumbers, $.;
		$numMappedReads++;
	}
	else	# Save header to print later
	{
		$header .= $_;
	}
	
	if ($. % 1000000 == 0)	# Print dots every 1M rows
	{
		print STDERR ".";
	}
}
close (SAM);

if ($numReadsWanted > $numMappedReads)
{
	print STDERR "\nHey dude and/or dudette!\n";
	print STDERR "$SAM_file has $numMappedReads rows of reads\n"; 
	print STDERR "  but you want $numReadsWanted rows.\n";
	print STDERR "Choose a subset of mapped reads that is less than $numMappedReads, OK?!?\n\n";
	exit;
}

print STDERR "I counted $numMappedReads rows in $SAM_file\n";

for ($i = 1; $i <= $numMappedReads; $i++)
{
	push @lineNumbers, $i;
}

# Do the shuffle
print STDERR "Shuffling reads by line number....\n";
fisher_yates_shuffle( \@lineNumbers );    # permutes @array in place
print STDERR "Done shuffling reads by line number.\n";

# Get the first $numReadsWanted rows
for ($i = 0; $i < $numReadsWanted; $i++)
{
	$lineNumbersWanted{$lineNumbers[$i]} = 1;
}

print STDERR "Going through $SAM_file and selecting random set of $numReadsWanted reads.";

# Start new file with header
print $header;

# Read SAM and print wanted rows
if ($SAM_file =~ /\.SAM/i)
{
	print STDERR "Going through the SAM file ($SAM_file) and selecting random set of $numReadsWanted reads ....\n";
	open (SAM, $SAM_file) || die "Cannot open $SAM_file) for reading: $!";
}
elsif ($SAM_file =~ /\.BAM/i)
{
	print STDERR "Going through the BAM file ($SAM_file) and selecting random set of $numReadsWanted reads ....\n";
	open (SAM, "samtools view -h $SAM_file |") || die "Cannot open $SAM_file) for reading: $!";
}
$lineNum = 0;
while (<SAM>)
{
	$lineNum++;
	
	if ($lineNumbersWanted{$lineNum})
	{
		if (! /^@/)	# Don't print header lines
		{
			print $_;
		}
	}
	if ($. % 1000000 == 0)	# Print dots every 1M rows
	{
		print STDERR ".";
	}
}
print STDERR "\n";

#####################

# From Perl Cookbook - 4.17. Randomizing an Array
# http://biocomputing5.wi.mit.edu/perl_books/cookbook/ch04_18.htm

sub fisher_yates_shuffle 
{
	my $array = shift;
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}
