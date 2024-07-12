#! /usr/bin/env perl
use warnings;


#######################################################################
# Copyright(c) 2011 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#
# Version: 1.2: 2 Dec 2010
# Version: 1.3: 3 Jan 2013
# Version: 2.0: 7 Jul 2020 -- Process BED of defined format into standard GTF file with different values for gene and transcript IDs.

# GFF format
# 1. seqname - The name of the sequence. Must be a chromosome or scaffold.
# 2. source - The program that generated this feature.
# 3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
# 4. start - The starting position of the feature in the sequence. The first base is numbered 1.
# 5. end - The ending position of the feature (inclusive).
# 6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
# 7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
# 8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
# 9. group - All lines with the same group are linked together into a single item. 

###
###  BED start and end coords are 0-based ? (In the UCSC browser, describe the first 100 nt as 0-100)
###  GFF start coords are 1-based, but end coords are 0-based ? (In the UCSC browser, describe the first 100 nt as 1-100)
###
###  See http://genome.ucsc.edu/FAQ/FAQformat.html for more info

########################################################################

if (! $ARGV[2])
{
	print STDERR "\nConvert BED into GTF.\n\n";
	print STDERR "USAGE: $0 bedFile source feature > gffFile\n";
	print STDERR "Ex: $0 foo.bed WIBR exon > foo.gff\n\n";
	print STDERR "Note that if the 4th BED field includes two semicolon-delimited values,\n  they will be used for gene and transcript IDs.\n";
	print STDERR "Otherwise the gene and transcript IDs will be set as identical values (from the 4th BED field).\n\n";
	exit;
}

$bedFile = $ARGV[0];
$SOURCE = $ARGV[1];
$FEATURE = $ARGV[2];

open (BED, $bedFile) || die "Cannot open $bedFile for reading: $!";
while (<BED>)
{
	if ($_ =~ /^type/){next;}	# Add on 3 Jan 2013

	chomp;
	($chr, $start, $end, $name, $score, $strand) = split (/\t/, $_);
	$gene_id = $name;
	$transcript_id = $name;
	
	# If the 4th field contains two "words" separated by a semicolon,
	# the first value will be used as gene ID and the second as transcript ID.
	if ($name =~ /(\S+); ?(\S+)/)
	{
		$gene_id = $1;
		$transcript_id = $2;
	}
	
	if (! $score) { $score = 1; }
	if (! $strand) { $strand = "."; }
	
	# GFF and BED use different coordinate systems !
	$start++;
	
	print "$chr\t$SOURCE\t$FEATURE\t$start\t$end\t$score\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\";\n";
}

