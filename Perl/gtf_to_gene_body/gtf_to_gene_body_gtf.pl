#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2018 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# George Bell, BaRC
# Version 1.0: 30 Nov 2018
# Version 1.1: 3 Jan 2018 -- Fix OBOB (start coord)
# Version 1.2: 4 Jan 2018 -- Require choice of gene body by gene or by transcript
# Version 1.3: 22 Nov 2022 -- Also accept gtf.gz files as input [GB]
#
#######################################################################

use List::Util qw(max min);

if (! $ARGV[1])
{
	print STDERR "\nMake a gene body GTF file from a GTF file,\n";
	print STDERR "  creating one gene body per transcript\n";
	print STDERR "  or one gene body per gene (from the union of transcript exons).\n";
	print STDERR "USAGE: $0 genes.gtf transcript|gene > gene_bodies.gtf\n\n";
	exit;
}

# Needs to be "transcript" or "gene"
$featureType = $ARGV[1];

if ($featureType eq "transcript")
{
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

		# Skip header lines
		if (! /^#/)
		{
			@f = split /\t/, $_;

			$chr = $f[0];
			$start = $f[3];
			$end = $f[4];
			$strand = $f[6];
			$featureDescription = $f[8];

			if ($featureDescription =~ /gene_id (\S+)/)
			{
				$geneID = $1;
				$geneID =~ s/\"//g;
				$geneID =~ s/;$//g; 
			}

			if ($featureDescription =~ /transcript_id (\S+)/)
			{
				$transcriptID = $1;
				$transcriptID =~ s/\"//g;
				$transcriptID =~ s/;$//g; 
			}
			# $transcript2gene{$transcriptID} = $geneID;

			if ($f[2] eq "exon")
			{
				# Skip any gene that has no "exon" line (GB -- 25 August 2016)
				$transcript2gene{$transcriptID} = $geneID;

				# Get source tag
				$transcript2source{$transcriptID} = $f[1];
				# Get score
				$transcript2score{$transcriptID} = $f[5];

				# print "Gene ($geneID) <= transcript ($transcriptID) ==> $chr\n";

				# Check that we're on the same chr and strand
				if ($transcript2chr{$transcriptID} && $transcript2chr{$transcriptID} ne $chr)
				{
					print STDERR "WARNING: Is $transcriptID ($geneID) on $transcript2chr{$transcriptID} or $chr ?\n";
				}
				elsif (! $transcript2chr{$transcriptID})
				{
					$transcript2chr{$transcriptID} = $chr;
				}

				if ($transcript2strand{$transcriptID} && $transcript2strand{$transcriptID} ne $strand)
				{
					print STDERR "WARNING: Is $transcriptID ($geneID) on $transcript2strand{$transcriptID} or $strand?\n";
				}
				elsif (! $transcript2strand{$transcriptID})
				{
					$transcript2strand{$transcriptID} = $strand;
				}

				if (! $transcript2start{$transcriptID})
				{
					$transcript2start{$transcriptID} = $start;
					$transcript2end{$transcriptID} = $end;
				}
				else
				{
					$transcript2start{$transcriptID} = min ($start, $transcript2start{$transcriptID});
					$transcript2end{$transcriptID} = max ($end, $transcript2end{$transcriptID});

					# print STDERR "For $transcriptID looking at $start of $start and end of 
				}
			}
		}
	}

	foreach $transcriptID (sort keys %transcript2gene)
	{
		# Prevent OBOB when converting from gtf to bed
		# $transcript2start{$transcriptID}--;

		# print "$transcript2chr{$transcriptID}\t$transcript2start{$transcriptID}\t$transcript2end{$transcriptID}\tgene_id \"$transcript2gene{$transcriptID}\"; transcript_id \"$transcriptID\";\t.\t$transcript2strand{$transcriptID}\n";
		print "$transcript2chr{$transcriptID}\t$transcript2source{$transcriptID}\t$featureType\t$transcript2start{$transcriptID}\t$transcript2end{$transcriptID}\t$transcript2score{$transcriptID}\t$transcript2strand{$transcriptID}\t.\tgene_id \"$transcript2gene{$transcriptID}\"; transcript_id \"$transcriptID\";\n";
	}
}

elsif ($featureType eq "gene")
{
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

		# Skip header lines
		if (! /^#/)
		{
			@f = split /\t/, $_;

			$chr = $f[0];
			$start = $f[3];
			$end = $f[4];
			$strand = $f[6];
			$featureDescription = $f[8];

			if ($featureDescription =~ /gene_id (\S+)/)
			{
				$geneID = $1;
				$geneID =~ s/\"//g;
				$geneID =~ s/;$//g; 
			}

			if ($featureDescription =~ /transcript_id (\S+)/)
			{
				$transcriptID = $1;
				$transcriptID =~ s/\"//g;
				$transcriptID =~ s/;$//g; 
			}
			# $transcript2gene{$transcriptID} = $geneID;

			if ($f[2] eq "exon")
			{			
				# Get source tag
				$gene2source{$geneID} = $f[1];
				# Get score
				$gene2score{$geneID} = $f[5];

				# print "Gene ($geneID) <= transcript ($transcriptID) ==> $chr\n";

				# Check that we're on the same chr and strand
				if ($gene2chr{$geneID} && $gene2chr{$geneID} ne $chr)
				{
					print STDERR "WARNING: Is $geneID on $gene2chr{$geneID} or $chr ?\n";
				}
				elsif (! $gene2chr{$geneID})
				{
					$gene2chr{$geneID} = $chr;
				}

				if ($gene2strand{$geneID} && $gene2strand{$geneID} ne $strand)
				{
					print STDERR "WARNING: Is $geneID on $gene2strand{$geneID} or $strand?\n";
				}
				elsif (! $gene2strand{$geneID})
				{
					$gene2strand{$geneID} = $strand;
				}

				if (! $gene2start{$geneID})
				{
					$gene2start{$geneID} = $start;
					$gene2end{$geneID} = $end;
				}
				else
				{
					$gene2start{$geneID} = min ($start, $gene2start{$geneID});
					$gene2end{$geneID} = max ($end, $gene2end{$geneID});

					# print STDERR "For $transcriptID looking at $start of $start and end of 
				}
			}
		}
	}

	foreach $geneID (sort keys %gene2chr)
	{
		# Prevent OBOB when converting from gtf to bed
		# $gene2start{$geneID}--;

		# print "$gene2chr{$geneID}\t$gene2start{$geneID}\t$gene2end{$geneID}\tgene_id \"$gene2gene{$geneID}\"; gene_id \"$geneID\";\t.\t$gene2strand{$geneID}\n";
		print "$gene2chr{$geneID}\t$gene2source{$geneID}\tgene\t$gene2start{$geneID}\t$gene2end{$geneID}\t$gene2score{$geneID}\t$gene2strand{$geneID}\t.\tgene_id \"$geneID\"; transcript_id \"$geneID\";\n";
	}
}

else
{
	print STDERR "\nERROR: Gene body type (argument 2) needs to be *transcript* or *gene*.\n\n";
	exit;
}
