#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2011-2012 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# For each query sequence, sort the BLAT hits to get the best one
#   by score (matches - 3 * mismatches).
# In the case of >1 hits with the same score, all are selected.
# What should be done where the best hit is questionable (with a weak score)?
# If you choose to filter, you require at least $FILTER_THRESHOLD % (set to 50%)
#   of the query sequence to be aligned to the genome.  If that's not true,
#   that alignment is ignored. 
#
# George Bell, BaRC
#
########################################################################

######  ??? Start by getting rid of all hits to random chromosomes  ######

# If filtering is selected, ignore sequences having a best hit that aligns 
#   less than $FILTER_THRESHOLD % of its length
$FILTER_THRESHOLD = 50;
$BAD_BLAT_EXT = "poor.blat";

if (! $ARGV[1] || $ARGV[1] eq ">")
{
	print STDERR "\nUSAGE:  $0 BLAT_file sorted_BLAT_results Filter[T/F]\n";
	print STDERR "Sorts a BLAT file by choosing the hit with the best score for each query.\n";
	print STDERR "Score = matches - 3 * mismatches\n";
	print STDERR "In the case of >1 hits with the same score, all are selected.\n";
	print STDERR "If filter is set as 'T', aligned bases must be at least $FILTER_THRESHOLD\% of length of query.\n\n";
	
	exit(0);
}

$blat_output = $ARGV[0];
$sorted_BLAT_results = $ARGV[1];
$filter = $ARGV[2];
if ($filter =~ /[Tt1]/)
{
	$filter = 1;
}
else
{
	$filter = 0;
}

$best_score = $last_query = 0;
$lineNum = 0;

# Open temp file (that will be made non-redundant at the end)
open (SORTED_BLAT, ">$sorted_BLAT_results.tmp") || die "cannot open $sorted_BLAT_results for writing: $!";

if ($filter) # 13 Jan 2011 (GB)
{
	# Open file for poor BLAT hits
	open (BAD_BLAT, ">$sorted_BLAT_results.$BAD_BLAT_EXT") || die "cannot open $sorted_BLAT_results.$BAD_BLAT_EXT writing: $!";
}

open (RAW_BLAT, $blat_output) || die "cannot open $blat_output for reading: $!";
while(<RAW_BLAT>)
{
	chomp($_);
	
	if (/^\d/)	# Skip header lines
	{
		@fields = split(/\t/, $_);
		$line = $_;
		$query = $fields[9];
		$score = $fields[0] - 3 * $fields[1];

		# If we're looking at data for the same query sequence again
		if ($query eq $last_query)
		{
			if ($score > $best_score)
			{
				@best_fields = @fields;
				$best_score = $score;

				# Make format like usual BLAT "psl" output
				$best_line = $line;

				# Make format like refFlat.txt UCSC annotation file
				# getRefFlatLine(); 
			}
			elsif ($score == $best_score)	# Two hits with the same score
			{
				$best_line .= "\n$line";	# Add this hit to the other
			}
		}
		else	# If we're looking at a query sequence for the first time
		{
			if ($best_line)
			{
				# Filtering
			
				if (! $filter)
				{
					print SORTED_BLAT "$best_line\n";
				}
				else
				{
					@best_lines = split (/\n/, $best_line);
					{
						foreach $one_of_best_lines (@best_lines)
						{
							@best_fields = split (/\t/, $one_of_best_lines);
							$queryLength = $best_fields[10];
							$matchLength = $best_fields[0];
							{
								if ((100 * $matchLength / $queryLength) > $FILTER_THRESHOLD)
								{
									print SORTED_BLAT "$one_of_best_lines\n";
								}
								else
								{
									print BAD_BLAT "$one_of_best_lines\n";
								}
							}
						}
					}
				}
			}
			@best_fields = @fields;
			$best_score = $score;
			$best_line = $line;

			# Make format like refFlat.txt UCSC annotation file
			# getRefFlatLine(); 
		}
		$last_query = $query;
	}
	else
	{
		if ($lineNum < 5)	# header lines
		{
			if (! $lineNum)	# first header line: Should say "psLayout version 3"
			{
				print SORTED_BLAT "$_\tSorted BLAT output - best hit(s) for each query sequence\n";
			}
			else
			{
				print SORTED_BLAT "$_\n";
			}
		}
	}
	$lineNum++;
}
print SORTED_BLAT "$best_line\n";

close (SORTED_BLAT);

# Start with $sorted_BLAT_results.tmp and print non-redundant lines to $sorted_BLAT_results
makeNonRedundant();

#########################  Subroutines  #########################

sub makeNonRedundant
{
	$numRemoved = 0;

	open (SORTED_BLAT, "$sorted_BLAT_results.tmp") || die "cannot open $sorted_BLAT_results.tmp for reading: $!";	
	open (FINAL_SORTED_BLAT, ">$sorted_BLAT_results") || die "cannot open $sorted_BLAT_results for writing: $!";	
	while (<SORTED_BLAT>)
	{		
		if (! $gotThisLine{$_})
		{
			print FINAL_SORTED_BLAT $_;
		}
		else
		{
			# print STDERR "Already got this hit:\n$_";
			$numRemoved++;
		}
		$gotThisLine{$_} = 1;
	}
	
	print STDERR "Sorted (and removed $numRemoved redundancies)\nFinal non-redundant sorted output is $sorted_BLAT_results\n";
	
	if ($filter)
	{ 
		print STDERR "Filtering performed:  \n";
		print STDERR "	Entries that have been filtered out are in $sorted_BLAT_results.$BAD_BLAT_EXT\n"; 
	}
	else
	{
		print STDERR "Filtering NOT performed.\n";
	}
	
	unlink("$sorted_BLAT_results.tmp");
}

sub getRefFlatLine 
{
	##################################################################
			
	# Make format like refFlat.txt UCSC annotation file
			
	$xm = $fields[9];
	if  ($xm =~ /(XM_\d*)/)
	{
		$xm = $1;
	}
	$chr = $fields[13];
	$strand = $fields[8];
	$txStart = $fields[15];
	$txEnd = $fields[16];
	$cdsStart = "?";
	$cdsEnd = "?";
	$exonCount = $fields[17];
	$exonStarts = $fields[20];
	@exonStarts = split (/,/, $exonStarts);
	@exonSizes = split (/,/, $fields[18]);
			
	$exonEnds = "";
	for ($i = 0; $i <= $#exonStarts; $i++)
	{
		$thisEnd = $exonStarts[$i] + $exonSizes[$i];
		$exonEnds .= "$thisEnd,";
	}

	$best_line = "\t$xm\t$chr\t$strand\t$txStart\t$txEnd\t$cdsStart\t$cdsEnd\t$exonCount\t$exonStarts\t$exonEnds\n";
			
	##################################################################
}
