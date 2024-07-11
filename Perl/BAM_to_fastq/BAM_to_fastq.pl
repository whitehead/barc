#! /usr/bin/env perl
use warnings;

##################################################################
#
# Copyright(c) 2012 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#
# Convert a single-read SAM file into a fastq file
# Version 1.1: Pay attention to strand of original sequence encoded in flag
#              and reverse complement if necessary
# Version 1.2 -- 18 October 2017 -- Improve handling of "read reverse strand" flag
# Version 1.3 -- 19 January 2017 -- Process both BAM and SAM
#
# Alternatives to this script: bamToFastq, samtools bam2fq, bam2fastx, bamtools convert
#
##################################################################


if (! $ARGV[0])
{
	print STDERR "\nConvert SAM to fastq\n";
	print STDERR "USAGE: $0 in.sam > out.fastq\n";
	print STDERR "\nNotes: Reads are given numbers as IDs.\n";
	print STDERR "       Reads with a '16' flag are reverse-complemented.\n\n";
	exit;
}

$sam = $ARGV[0];

# Use this as the sequence ID
$seqNum = 0;

if ($sam =~ /\.SAM/i)
{
	open (SAM, $sam) || die "Cannot open $sam for reading: $!";
} 
elsif ($sam =~ /\.BAM/i)
{
	open (SAM, "samtools view $sam |") || die "Cannot open $sam for reading: $!";
}
while (<SAM>)
{
	# SRR015084.905423	73	2L	50799	0	37M	*	0	0	CGGAAAGTCTCTGTCGTCAGGTTGAAGAATACTTGGA	IAHIII?IIAIIFIII:IIIIIIIIIIIIIIIIIIII	NM:i:0


	if (! /^\@/)	# Skip header rows
	{
		$seqNum++;
	
		@f = split /\t/, $_;
		
		$flag = $f[1];
		$seq = $f[9];
		$qual = $f[10];
		
		# if ($flag >= 16)	# Reverse strand; need to reverse complement
		if ($flag & hex("0x10"))	# Reverse strand; need to reverse complement
		{
			$seq = reverse($seq);
    		$seq =~ tr/ACGTacgt/TGCAtgca/;
		
			# Reverse the quality scores too
			$qual = reverse($qual);
		}
		
		
		print "@" . $seqNum . "\n" . $seq . "\n";
		print "+" . $seqNum . "\n" . $qual . "\n";
	}
		
}