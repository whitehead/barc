#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2017 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Count CRISPR guides in a pooled screen (one sample) 
# Guides are expected to occur at beginning of a read sequence (and need to perfectly match).
# 
# Input files:  guidesFile (see My_guides_list.txt) and gzipped fastq file (see Reads_sample.fq.gz).
#
# Author: George Bell
#
# Version 1.0 -- 8 September 2017
# Version 1.1 -- 24 Nov. 2020 BY
#
#######################################################################


if (! $ARGV[2])
{
	print STDERR "\nRead compressed fastq file and count guides\n"; 
	print STDERR "USAGE $0 guidesFile Reads.fq.gz outputFile output\n\n";
	exit;
}

$guidesFile = $ARGV[0];
$fastqGz = $ARGV[1];
$output = $ARGV[2];
$output_guide = $output . ".txt";
$output_not_guide = $output . "_other_abundent_guides.txt";

# Expect that the guide sequence is found in this column of the guildes file
$acceptedGuidesSeqField = 2;
# After counting our guides, we'll print other sequences found that weren't in the guides list
# This may be helpful for quality control
# Print any not-on-list guide-length sequences with at least $minSeqCountsToPrint counts
$minSeqCountsToPrint = 1000;

###
###  Count first @crisprGuideLengths nts of each read sequence
###

print STDERR "\nGetting accepted guides from $guidesFile.\n";

# Read guides file to get sequence length(s) that we'll need to match
open(GUIDES, $guidesFile) || die "Cannot open $guidesFile for reading: $!";
while (<GUIDES>)
{
	chomp;
	if ($. > 1)	# Skip header
	{
		@f = split /\t/, $_;
		$isGuideLength{ length $f[$acceptedGuidesSeqField - 1] } = 1;
	}
}
@crisprGuideLengths = sort {$a <=> $b} keys %isGuideLength;
close(GUIDES);

print STDERR "Reading $fastqGz and counting guides (first @crisprGuideLengths nts) ....\n";
open(IN, "gunzip -c $fastqGz |") || die "Can’t open pipe to $fastqGz";
while (<IN>)
{
	chomp;
	if (($. -2) % 4 == 0)
	{
		foreach $crisprGuideLength (@crisprGuideLengths)
		{
			$guideSeq = uc(substr($_, 0, $crisprGuideLength));
			# print "$. => $guideSeq\n";
			$guideToCounts{$guideSeq}++;
		}
	}
}

print STDERR "Printing output ....\n";
print "The number of reads for each guide is printing to $output_guide .\n";

open(OUT, ">$output_guide") || die "Cannot write to $output_guide: $!";

open(GUIDES, $guidesFile) || die "Cannot open $guidesFile for reading: $!";
while (<GUIDES>)
{
	# ID	Guide
	# Gene_A_1	GCCTACTACAGAGAAGCTG
	# Gene_A_2	AGGGCAGAATCCACAAACTG

	chomp;
	if ($. == 1)
	{
		print OUT "$_\t$fastqGz\n";
	}
	else
	{
		@f = split /\t/, $_;
		$guideSeq = uc($f[$acceptedGuidesSeqField - 1]);
		$isAcceptedGuide{$guideSeq} = 1;
		if (! $guideToCounts{$guideSeq})
		{
			$guideToCounts{$guideSeq} = 0;
		}
		print OUT "$_\t$guideToCounts{$guideSeq}\n";
	}
}

close(OUT);

#print "\n==================== Other apparent guides not in expected list =================\n\n";

open(OUT2, ">$output_not_guide") || die "Cannot write to file: $!";

foreach $guide (sort { $guideToCounts{$b} <=> $guideToCounts{$a} } keys %guideToCounts)
{
	if ($guideToCounts{$guide} && ! $isAcceptedGuide{$guide}) # Skip 0s
	{
		if ($guideToCounts{$guide} >= $minSeqCountsToPrint)
		{
			print OUT2 "$guide\t$guideToCounts{$guide}\n";
		}
		else
		{
			last;
		}
	}
}
close (OUT2);

print "Other apparent guides not in expected list, but with more than $minSeqCountsToPrint is printed in $output_not_guide\n";
#print "Other apparent guides with less than $minSeqCountsToPrint counts have not been printed in $output_not_guide .\n";
print STDERR "\nAll done!\n\n";
