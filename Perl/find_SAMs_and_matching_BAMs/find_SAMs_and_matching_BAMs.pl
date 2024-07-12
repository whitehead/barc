#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2020 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Recursively search a desired directory for SAM files and matching BAM files
# Possible matches: Foo.sam => Foo.bam or Foo.sorted.bam
#
# Author: George W Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 1.0 -- 11 June 2020
#
###########################################################################

use File::Find;

if (! $ARGV[0])
{
	print STDERR "\nLook for SAM and matching BAM files (recursively) in a selected path.\n";
	print STDERR "Also print commands to delete (if BAM is present) or convert (if BAM is missing) SAM.\n\n";
	print STDERR "USAGE: $0 /path/to/search\n\n";
	exit;
}

$path = $ARGV[0];

# Commands that user can act on (but NOT executed by this script)
@samToBamCommands = ();
@deleteSamCommands = ();

print "SAM files found (and matching BAM(s), if also found):\n\n";

find( \&lookAtThisFile, $path);

print "\n# If matching BAM file exists, delete SAM file.\n";
print join "\n", @deleteSamCommands;

print "\n\n# If matching BAM file doesn't exist, create it.\n";
print "# Note that these commands will fail if the SAM file doesn't have a header.\n";
print join "\n", @samToBamCommands;
print "\n";

sub lookAtThisFile 
{
	$fileName = $File::Find::name;
	if ($fileName =~ /\.SAM$/i)
 	{
 		@matchingBAM = ();
 		
		$matchingBAMname_1 = $fileName;
		$matchingBAMname_2 = $fileName;
		$matchingBAMname_1 =~ s/\.SAM$/.BAM/i;
		$matchingBAMname_2 =~ s/\.SAM$/.bam/i;
 		
 		if (-e $matchingBAMname_1)
 		{
 			push @matchingBAM, $matchingBAMname_1;
 		}
 		if (-e $matchingBAMname_2)
 		{
 			push @matchingBAM, $matchingBAMname_2;
 		}
 		
		$matchingBAMname_1 = $fileName;
		$matchingBAMname_2 = $fileName;
		$matchingBAMname_1 =~ s/\.SAM$/.sorted.BAM/i;
		$matchingBAMname_2 =~ s/\.SAM$/.sorted.bam/i;

		if (-e $matchingBAMname_1)
		{
			push @matchingBAM, $matchingBAMname_1;
		}
		elsif (-e $matchingBAMname_2)
		{
			push @matchingBAM, $matchingBAMname_2;		 	
		}
 		
 		if (@matchingBAM)
 		{
 			print "SAM+BAM: $fileName <==> @matchingBAM\n";
 			$deleteSamCommand = "rm -f $fileName";
 			push @deleteSamCommands, $deleteSamCommand;
 		}
 		else
 		{
 			print "SAM only: $fileName\n";
 			
 			$matchingBAMname = $fileName;
 			$matchingBAMname =~ s/\.SAM$/.bam/i;
 			$samToBamCommand = "samtools view -bS $fileName > $matchingBAMname";
			push @samToBamCommands, $samToBamCommand;
 		}
 	}
	return;
}
