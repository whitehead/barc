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
   print STDERR "\nParse a BLAST file into tab-delimited fields like with BLAST's -m 8 output\n";
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

print "Query id	Subject id\t% identity\talignment length\tmismatches\tgap openings\tq. start\tq. end\ts. start\ts. end\te-value\tbit score\n";


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
			# $queryName = $result->query_name;
			$queryAcc = $result->query_accession;
			$hitAcc = $hit->accession;
			$hspLength = $hsp->hsp_length;			
			$hspEvalue = $hsp->evalue;
			$fracIdentical = $hsp->frac_identical;
			$numGaps = $hsp->gaps;
			$hspStart = $hsp->start('hit');
			$hspEnd = $hsp->end('hit');
			$hspBits = $hsp->bits;

			# Need to check this to reverse coords if necessary 			
			$queryStrand = $hsp->strand('query');
			
			$hspNumIdentical = $hsp->num_identical;
			$numMismatches = $hspLength - $hspNumIdentical;
			
			if ($queryStrand eq "1")
			{
				$queryStart = $hsp->start('query');
				$queryEnd = $hsp->end('query');
			}
			else
			{
				$queryStart = $hsp->end('query');
				$queryEnd = $hsp->start('query');			
			}

			print "$queryAcc\t";
			print "$hitAcc\t";
			print sprintf("%.2f", $fracIdentical * 100), "\t";
			print "$hspLength\t";

			print "$numMismatches\t";
			print "$numGaps\t";
			print "$queryStart\t";
			print "$queryEnd\t";
			print "$hspStart\t";
			print "$hspEnd\t";
			print "$hspEvalue\t";
			print "$hspBits\n";
		}
	}
}

