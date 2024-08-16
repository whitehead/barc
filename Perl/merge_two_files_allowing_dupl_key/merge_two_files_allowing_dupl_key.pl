#!/usr/bin/perl

use warnings;

#######################################################################
#
# Copyright(c) 2023 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Description: Join two files by the first column
#
# Author: Bingbing Yuan, Bioinformatics and Research Computing, wibr-bioinformatics@wi.mit.edu 
#
# USAGE: join2filesByFirstColumn.pl file1.txt file2.txt
#
# Version: 1.0 -- Created Jan 2023
#
###########################################################################


use strict;

if (! $ARGV[4]) {
  print "Usage:
$0 file.txt anno_file 1st_file_id_col out_col field_seperator\n
Note: 	
Assume both input text files are tab-delimited.
The first column in the 2nd file must be unique and will be used as the key (ID) to match with the ID in the 1st file, specified as the 3rd argument
The corresponding value from the user-defined output column of the 2nd file, the 4th argument,  will be appended to the 1st file.
If a record in the 1st file contains multiple IDs, separated by a user-defined delimiter (specified as the last argument), each ID will be used for matching against the 2nd file.
The resulting merged file will maintain the same order of IDs as in the 1st file.
.\n"; 
 
exit;
}

my $SEP=$ARGV[4];
#print "SEPARATOR=$SEP=\n";
my $OUT_COL = $ARGV[3] -1;



my ($db2_href, $Size_2) = fileToHashOfArrayBySep($ARGV[1]);

open(FIRSTFILE, $ARGV[0]) || die "can not read $ARGV[0] $!\n";
while (<FIRSTFILE>) {
  my $found = 0;
  #print $_;
  chomp();
  # ignore comments and empty spaces
  next if ($_ =~ /^\#/);
  next if ($_ =~ /^\s*$/);
  # split record by tab
  my @arr = split(/\t/,  $_, -1);
  my $id_pos = $ARGV[2] -1;
  my @ids = split(/$SEP/, $arr[$id_pos]);
  my @out_items = ();
  #print "AAAAA ", join("==", @ids), "\n";
  
  foreach my $id (@ids) 
  {
  	if ( $db2_href ->{ $id } )
  	{
	  $found = 1;
	  push @out_items, $db2_href ->{ $id } -> [$OUT_COL];
    	}
  }
  # print "\n\n";
  if (! $found)
  {
    print join("\t", @arr), "\t", "\n";
  }
  else
  {
   print join("\t", @arr), "\t", join($SEP, @out_items), "\n";
  }
}
close(FIRSTFILE);


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

    #0610006A03Rik   2.1000

    # split each record by tab
    my @arr = split(/\t/, $_, -1);
    my $id  = $arr[0];
    
    if (! defined $hash{ $id } )
    {
      @{ $hash{$id} } = ();
    }
    push @{ $hash{ $id } }, @arr;
    
    if (! $size)
    {
      $size = $#arr;
    }
  }
  close (FL) || die "Can not close $file\n";
  return (\%hash, $size);
}
