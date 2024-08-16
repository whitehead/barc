#! /usr/bin/env perl


#######################################################################
# Copyright(c) 2005 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Bingbing Yuan
#         Bioinformatics and Research Computing
#         barc@wi.mit.edu
#
###########################################################################


use warnings;
use Tie::IxHash;

# purpose:
# find the shared symbol (1st field) among multiple tab-delimited files,
# IDs such as gene symbols should be in the first column of all the files to be merged.
#          
# input file format: 
# multiple tab-delimited files
# Output:
# symbol order is the same as gene symbol (IDs) in the first file
# if a symbol in the first file does not exist in the 2nd file,
# 'UNMATCHED' will be added for this symbol
#
#



use strict;

if (! $ARGV[1]) {
  print "\nJoin two or more files using matching IDs in the first column of each file.
Usage: $0 file_1.txt file_2.txt > file.1+2.txt\n
Notes: 	
	Input files must be tab-delimited.
	The first column in each file should contain the IDs to be matched.
	IDs are case-sensitive and should be unique within each file.
	The order of rows in the merged file will match the order of IDs from the first input file.
	The order of columns in the merged file will follow the sequence specified in the command line.
	If an ID from the first file is not found in a file, the corresponding record will have 'UNMATCHED' in place of missing data.\n\n";
exit;
}

my $UNKNOWN = "UNMATCHED";

my $FILE_COUNT = $#ARGV;

my %JOINED = ();
tie %JOINED, 'Tie::IxHash';

my ($db2_href, $Size_2) = fileToHashOfArrayBySep($ARGV[1]);

open(FIRSTFILE, $ARGV[0]) || die "can not read $ARGV[0] $!\n";
while (<FIRSTFILE>)
  {
    my $found = 0;
    #    print $_;
    chomp();
    # ignore comments and empty spaces
    next if ($_ =~ /^\#/);
    next if ($_ =~ /^\s*$/);
    $_ =~ s/\015//g;

    # split record by tab
    my @arr = split(/\t/,  $_, -1);
    my $id = shift(@arr);
    
    if ( $db2_href ->{ $id } )
      {
	$found =1;
	for my $i (0..$#{ $db2_href->{$id} } )
	  {
	    #print "id=$id", "arr=", @arr, "other=", $db2_href->{$id}->[$i], "\n";
	    $JOINED{$id} = [@arr, $db2_href->{$id}->[$i] ];
	  }
      }
    if (! $found)
      {
	my @unknown = ();
	push @unknown, $UNKNOWN;
	# print the missing columns with blank
	for my $i (1..$Size_2)
	  {
	    push @unknown, "";
	  }
	$JOINED{$id} = [@arr, @unknown];
	
      }
  }
close(FIRSTFILE);


# To join the rest of files
if($#ARGV > 1) 
  {
    for my $f (2..$#ARGV) 
      {
	my ($db2_href, $Size_2) = fileToHashOfArrayBySep($ARGV[$f]);
	
	foreach my $id ( keys %JOINED )
	  {
	    my $found = 0;
	    my @arr = @{ $JOINED{$id} };
	    # print join("\t", @arr), "\n";
	    if ( $db2_href ->{ $id } )
	      {
		# print "matched\n";
		for my $i (0..$#{ $db2_href->{$id} } )
		  {
		    $JOINED{$id} = [@arr, $db2_href->{$id}->[$i] ];
		  }
		$found =1;
	      }
	    # print "\n\n";
	    if (! $found)
	      {
		# print "NOT matched\n";
		my @unknown = ();
		push @unknown, $UNKNOWN;
		# print the missing columns with blank
		for my $i (1..$Size_2)
		  {
		    push @unknown, "";
		  }
		$JOINED{$id} = [@arr, @unknown];
		
	      }
	  }
      }
  }

# print results
foreach my $gene (keys %JOINED)
  {
    print $gene, "\t", join("\t", @{$JOINED{$gene}} ), "\n";
  }












sub fileToHashOfArrayBySep {
  my $file = shift;
  my %hash = ();
  my $size;
  open (FL, $file) || die "Can not open $file\n";
  while (<FL>) {
#    print $_;
    chomp();
    # ignore comments and empty lines
    next if ($_ =~ /^\#/);
    next if ($_ =~ /^\s*$/);
    $_ =~ s/\015//g;

    #0610006A03Rik   2.1000

    # split each record by tab
    my @arr = split(/\t/, $_, -1);
    my $id  = shift(@arr);
    
    if (! defined $hash{ $id } )
    {
      @{ $hash{$id} } = ();
    }
    push @{ $hash{ $id } }, join("\t", @arr);
    
    if (! $size)
    {
      $size = $#arr;
    }
  }
  close (FL) || die "Can not close $file\n";
  return (\%hash, $size);
}



