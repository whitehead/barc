#! /usr/bin/env perl
use warnings;

#######################################################################
# Copyright(c) 2011 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Parsing BLAST reports with BioPerl's Bio::SearchIO module
# WI Unix and Programming Skills for Biologists - class 3
#
# See documentation at http://www.bioperl.org/HOWTOs/html/SearchIO.html
# Last updated by George Bell -- March 2007
#
###########################################################################

use Bio::SearchIO;

# Prompt the user for the filename if it's not given as an argument
if (! $ARGV[0])
{
	print STDERR "\nParse a BLAST file into tab-delimited fields\n";
	print STDERR "USAGE: $0 blastFile > parsed_output\n\n";
	exit;
}
else
{
	$inFile = $ARGV[0];
}

# Create a SearchIO "object" from the BLAST report
# This report may contain multiple concatenated reports
$report = new Bio::SearchIO( -file => "$inFile", -format => "blast");

# Print out one line of all desired fields, delimited by tabs:
# QUERY NAME, HIT NUM, HIT DESCRIPTION, HIT E-VALUE, FRACTION IDENTITY   [and ending with a newline]

# 1

print "queryName\t";
print "queryAcc\t";
print "queryLength\t";
print "queryDesc\t";
print "dbName\t";
print "numSeqsInDb\t";
print "numHits\t";
print "hitNumber\t";

print "hitName\t";
print "hitAcc\t";
print "hitDesc\t";
print "hitEvalue\t";
print "hitBits\t";
print "hitLength\t";
print "numHsps\t";

print "hspEvalue\t";
print "fracIdentical\t";
print "fracConserved\t";
print "numGaps\t";
print "hspLength\t";
print "hspQuerySeq\t";
print "hspHitSeq\t";
print "hspConsensus\t";
print "hspRank\t";

print "queryStrand\t";
print "hitStrand\t";
print "queryStart\t";
print "queryEnd\t";
print "hspStart\t";
print "hspEnd\t";

print "hspScore\t";

# GB - 4 Feb 2010
# print "hspBits\n";
print "hspBits\t";
print "hspQueryFrame[1,2,3]\n";

# Go through BLAST reports one by one		     
while($result = $report->next_result)
{
	# Reset hit number for each query  
	$hitNumber = 0;

	# Go through each each matching sequence
	while($hit = $result->next_hit)
	{
		$hitNumber++;

		# Go through each each HSP for this sequence
		while ($hsp = $hit->next_hsp)
		{
			# Print some tab-delimited data about this HSP
			$queryName = $result->query_name;
			$queryAcc = $result->query_accession;
			$queryLength = $result->query_length;
			$queryDesc = $result->query_description;
			$dbName = $result->database_name;
			$numSeqsInDb = $result->database_entries;
			$numHits = $result->num_hits;

			$hitName = $hit->name;
			$hitAcc = $hit->accession;
			$hitDesc = $hit->description;
			$hitEvalue = $hit->significance;
			$hitBits = $hit->bits;
			$hitLength = $hit->length;
			$numHsps = $hit->num_hsps;

			$hspEvalue = $hsp->evalue;
			$fracIdentical = $hsp->frac_identical;
			$fracConserved = $hsp->frac_conserved;
			$numGaps = $hsp->gaps;
			$hspLength = $hsp->hsp_length;
			$hspQuerySeq = $hsp->query_string;
			$hspHitSeq = $hsp->hit_string;
			$hspConsensus = $hsp->homology_string;
			$hspLength = $hsp->hsp_length;
			$hspRank = $hsp->rank;

			$queryStrand = $hsp->strand('query');
			$hitStrand = $hsp->strand('hit');
			$queryStart = $hsp->start('query');
			$queryEnd = $hsp->end('query');
			$hspStart = $hsp->start('hit');
			$hspEnd = $hsp->end('hit');

			$hspScore = $hsp->score;
			$hspBits = $hsp->bits;
			
			# GB - 4 Feb 2010
			$hspQueryFrame = $hsp->query->frame + 1;
			
			print "$queryName\t";
			print "$queryAcc\t";
			print "$queryLength\t";
			print "$queryDesc\t";
			print "$dbName\t";
			print "$numSeqsInDb\t";
			print "$numHits\t";
			
			print "$hitNumber\t";

			print "$hitName\t";
			print "$hitAcc\t";
			print "$hitDesc\t";
			print "$hitEvalue\t";
			print "$hitBits\t";
			print "$hitLength\t";
			print "$numHsps\t";

			print "$hspEvalue\t";
			print "$fracIdentical\t";
			print "$fracConserved\t";
			print "$numGaps\t";
			print "$hspLength\t";
			print "$hspQuerySeq\t";
			print "$hspHitSeq\t";
			print "$hspConsensus\t";
			print "$hspRank\t";

			print "$queryStrand\t";
			print "$hitStrand\t";
			print "$queryStart\t";
			print "$queryEnd\t";
			print "$hspStart\t";
			print "$hspEnd\t";

			print "$hspScore\t";
			
			# GB - 4 Feb 2010
			# print "$hspBits\n";
			print "$hspBits\t";
			print "$hspQueryFrame\n";
		}
	}
}

