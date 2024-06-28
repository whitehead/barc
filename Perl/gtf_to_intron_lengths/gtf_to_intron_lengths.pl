#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2013 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Make a list of intron lengths for all transcripts of a genome from a GTF file
# George W. Bell, BaRC
# Version 1.0: 19 July 2013
# Version 1.1: 12 August 2013 -- Modify code to work with GTF files that include exon number
# Version 1.2: 22 Nov 2022 -- Also accept gtf.gz files as input [GB]
#
#######################################################################

if (! $ARGV[0])
{
	print STDERR "\nGet intron lengths from a GTF file\n";
	print STDERR "USAGE: $0 genes.gtf > intron_lengths.txt\n\n";
	exit;
}

$gtfFile = $ARGV[0];
if ($gtfFile =~ /\.gz/i)
{
	open (GTF, "zmore $gtfFile |") || die "Cannot open $gtfFile for reading: $!";
}
else
{
	open (GTF, "$gtfFile") || die "Cannot open $gtfFile for reading: $!";
}

while (<GTF>)
{
	chomp;
	
	@f = split /\t/, $_;
	
	$featureType = $f[2];
	$start = $f[3];
	$stop =  $f[4];
	# $strand = $f[6];
	$gene_transcript = $f[8];
	
	if ($gene_transcript =~ /gene_id (\S+)/)
	{
		$geneID = $1;
		$geneID =~ s/\"//g;
		$geneID =~ s/;$//g; 
	}

	if ($gene_transcript =~ /transcript_id (\S+)/)
	{
		$transcriptID = $1;
		$transcriptID =~ s/\"//g;
		$transcriptID =~ s/;$//g; 
	}
	
	if ($geneID && $transcriptID)
	{
		$gene_transcript = "$geneID\t$transcriptID";
	}
	else
	{
		print STDERR "WARNING: Missing transcript or gene id ($_)\n";
		# exit;
	}
	
	# Only look at exons
	if ($featureType eq "exon")
	{
		push(@{$transcript2starts{$gene_transcript}}, $start);
		$transcriptStarts2ends{$gene_transcript}{$start} = $stop;
		
		# Each transcript should have only one strand (and only one chr, so let's skip that)
		# $transcript2strand{$gene_transcript} = $strand;
		
		# print "$gene_transcript => $start to $stop\n";
	}
}

foreach $gene_transcript (sort keys %transcriptStarts2ends)
{
	@startsSorted = sort { $a<=>$b } @{$transcript2starts{$gene_transcript}};
	
	($geneID, $transcriptID) = split /\t/, $gene_transcript;

	for ($i = 0; $i < $#startsSorted; $i++)
	{
		$thisEnd = $transcriptStarts2ends{$gene_transcript}{$startsSorted[$i]};
		$nextStart = $startsSorted[$i+1];
		$intronLength = $nextStart - $thisEnd;
		# print "$gene_transcript\t$nextStart - $thisEnd = $intronLength\n";
		
		print "$geneID\t$transcriptID\t$intronLength\n";
	}
}
