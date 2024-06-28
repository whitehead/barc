#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2014 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# George Bell, BaRC
# Version 1.0: 11 Dec 2014
# Version 1.1: Change output file to be sorted by gene ID and add variable resets;
#              Add gene lengths files to outputs (14 June 2017, GB)
# Version 1.2: 22 Nov 2022 -- Also accept gtf.gz files as input [GB]
#
#######################################################################

use Statistics::Lite qw(:all);
use File::Basename;
 
if (! $ARGV[0])
{
	print STDERR "\nGet transcript and gene lengths from a GTF file\n";
	print STDERR "USAGE: $0 genes.gtf > transcript_lengths.txt\n\n";
	exit;
}

$gtfFilePath = $ARGV[0];

# Make output file names for gene lengths files
$gtfFile = basename($gtfFilePath);
$MeanGeneLengthsFile = "$gtfFile.meanGeneLength.txt";
$MedianGeneLengthsFile = "$gtfFile.medianGeneLength.txt";

if ($gtfFile =~ /\.gz/i)
{
	open (GTF, "zmore $gtfFilePath |") || die "Cannot open $gtfFile for reading: $!";
}
else
{
	open (GTF, "$gtfFilePath") || die "Cannot open $gtfFile for reading: $!";
}

while (<GTF>)
{
	if (! /^#/)
	{
		chomp;
		$geneID = "";
		$transcriptID = "";

		@f = split /\t/, $_;

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

		if ($geneID && $transcriptID)
		{	
			$exonWidth = $f[4] - $f[3] + 1;

			$transcript2gene{$transcriptID} = $geneID;
			push @{ $gene2transcripts{$geneID} }, $transcriptID;
			# print STDERR "$geneID => $transcriptID\n";

			if ($f[2] eq "exon")
			{
				$transcript2length{$transcriptID} += $exonWidth;
				# print STDERR "Adding $exonWidth to $transcriptID\n";
			}
		}
	}
}

# Process one gene at a time
for $geneID ( keys %gene2transcripts ) 
{
    @transcriptsThisGene = @{ $gene2transcripts{$geneID} };
    %seen = ( );
	@transcriptsThisGene_nr = ( );
    foreach $item (@transcriptsThisGene) 
    {
	    unless ($seen{$item}) 
	    {
	        # if we get here, we have not seen it before
	        $seen{$item} = 1;
	        push(@transcriptsThisGene_nr, $item);
	    }
	}
    
    @transcriptLengths = ();
    foreach $transcriptID (@transcriptsThisGene_nr)
    {
    	push @transcriptLengths, $transcript2length{$transcriptID};
    }
    $gene2meanLength{$geneID}   = sprintf("%.0f", mean @transcriptLengths);
    $gene2medianLength{$geneID} = sprintf("%.0f", median @transcriptLengths);
    # print "$geneID\t@transcriptsThisGene_nr\t@transcriptLengths\n";
}

$MeanGeneLengthsFile = "$gtfFile.meanGeneLength.txt";
$MedianGeneLengthsFile = "$gtfFile.medianGeneLength.txt";

open (MEAN_GENE, ">$MeanGeneLengthsFile") || die "Cannot open $MeanGeneLengthsFile for writing: $!";
open (MEDIAN_GENE, ">$MedianGeneLengthsFile") || die "Cannot open $MedianGeneLengthsFile for writing: $!";

# Sort genes by ID and print them in that order
@transcriptsSortedByGene = sort { $transcript2gene{$a} cmp $transcript2gene{$b} } keys %transcript2gene;
# Print header on main output file
print "Gene ID\tTranscript ID\tTranscript length\n";
foreach $transcriptID (@transcriptsSortedByGene)
{
	$geneID = $transcript2gene{$transcriptID};
	print "$geneID\t$transcriptID\t$transcript2length{$transcriptID}\n";
	
	if (! $printedThisGene{$geneID})
	{
		print MEAN_GENE "$geneID\t$gene2meanLength{$geneID}\n";
		print MEDIAN_GENE "$geneID\t$gene2medianLength{$geneID}\n";
	}
	
	# For MEAN_GENE and MEDIAN_GENE files
	$printedThisGene{$geneID} = 1;
}

print STDERR "\nAll done!  For gene lengths, see files
 - $MeanGeneLengthsFile
 - $MedianGeneLengthsFile
Note that gene lengths are calculated as the mean or median of all transcripts lengths of that gene.\n\n";

