#! /usr/bin/env perl
use warnings;


###############################################################################
#
# Copyright(c) 2008 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Prathapan Thiru
# Version: 1.0
# History: Last modified Sep 11, 2008
#
# Comment: 
#   This script converts user-listed text files into one Excel file
#   with each sheet representing a text file.
#
###############################################################################


use strict;
use Spreadsheet::WriteExcel;


if (! $ARGV[0])
{
	printUsageAndExit();
}


#delimiter - change if different from tab
my $delimiter = '\t';

#file with a list of text file names
my $fofFile = $ARGV[0];

#output file name
my $outFile = $ARGV[1];

open (IN, $fofFile) || die "$fofFile does not exist\n";

#Create workbook
my $workbook = Spreadsheet::WriteExcel->new("$outFile");

#Go thru each file in the $fofFile and process the data
while (<IN>){
	chomp $_;
	open (IN2, $_) || die "$_ does not exist\n";
	
	my $worksheet = $workbook->add_worksheet("$_");
	my @xlsData;
	#Add data from individual file into a sheet
	while (<IN2>){
		my @line = split ($delimiter, $_);
		push (@xlsData, [@line]);
	}
	#Print out to sheet
	$worksheet->write_col('A1',\@xlsData);
}
	
	


# If no argument is given, print usage and exit;

sub printUsageAndExit
{
	print STDERR "\nUsage: txt2xls.pl <file with list of file names> <output file name>\n\n";	
	exit;
}

