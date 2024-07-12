#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2014-2021 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Given a tab-delimited file matrix, select desired rows and columns.
#
# George Bell, BaRC
# Version 1.1: 2/6/14 	# Removed $dir = "."; [GB]
# Version 1.2: 8/17/21 	# Let program run even if "columns to get" includes items not in matrix [GB]
#
#######################################################################

$exp_data_file = $ARGV[0];
$row_file = $ARGV[1];
$col_file = $ARGV[2];

if (! $exp_data_file && (! ($row_file || $col_file)))
{
	print STDERR "\nSubmatrix selector (select rows and/or columns from a matrix)\n";
	print STDERR "USAGE: $0 matrixFile rowIDfile columnIDfile > submatrixFile[output]\n";
	print STDERR "To select only by row IDs:\n";
	print STDERR "       $0 matrixFile rowIDfile > submatrixFile[output]\n\n";
	print STDERR "To select only by columns IDs:\n";
	print STDERR "       $0 matrixFile 0 columnIDfile > submatrixFile[output]\n\n";
	exit;
}

#####################  analyze  #####################

if ($row_file)	# Only select certain rows 
{
	open  (ROW_LIST, "$row_file") || die "Cannot open ROW_LIST file $row_file: $!";
	while (<ROW_LIST>)
	{
		chomp($_);

		s/^\s+//;
		s/\s+$//;
		s/\r\n/\n/g; # Convert Windows files
		s/\r/\n/g; # Convert Mac files
		
		if ($_)
		{
			push (@rows, $_);
			$rowsToGet{$_} = 1;
		}
	}
	close (ROW_LIST);
}

if ($col_file)	# Only select certain columns 
{
	open  (COL_LIST, "$col_file") || die "Cannot open COL_LIST file $col_file: $!";
	while (<COL_LIST>)
	{
		chomp($_);

		s/^\s+//;
		s/\s+$//;
		s/\r\n/\n/g; # Convert Windows files
		s/\r/\n/g; # Convert Mac files

		if ($_)
		{
			push (@cols, $_);
		}
	}
	close (COL_LIST);
	
	# Get column IDs from matrix file
	open  (EXP_DATA, "$exp_data_file") || die "Cannot open EXP_DATA file $exp_data_file: $!";
	while (<EXP_DATA>)
	{
		if ($. == 1)
		{
			chomp($_);

			s/^\s+//;
			s/\s+$//;
			s/\r\n/\n/g; # Convert Windows files
			s/\r/\n/g; # Convert Mac files

			@columnIDsInMatrix = split (/\t/, $_);
		}
	}
	close(EXP_DATA);
	
	# Added 8/17/2021 to remove columns in "to get" list but not present in matrix
	@cols = get_A_in_B(\@columnIDsInMatrix, \@cols);
}

# Read through big matrix file

open  (EXP_DATA, "$exp_data_file") || die "Cannot open EXP_DATA file $exp_data_file: $!";
$lineNum = 0;
while (<EXP_DATA>)
{
	chomp($_);
	
	s/^\s+//;
	s/\s+$//;
	s/\r\n/\n/g; # Convert Windows files
	s/\r/\n/g; # Convert Mac files

	@fields = split (/\t/, $_);
	
	# To get rid of double quotes added by Excel
	if (/\".+\"/)
	{
		for ($z = 0; $z <= $#fields; $z++)
		{
			if ($fields[$z] =~ /^\".*\"$/)
			{
				$fields[$z] =~ s/^\"(.*)\"$/$1/;
			}
		}
		$_ = join ("\t", @fields);
	}
	
	if (! $lineNum)	# first line; get header info
	{
		if (! $col_file)	# Print all columns
		{
			print "$_\n";
		}
		else	# Print headers for selected columns
		{			
			$colNamesTabbed = join ("\t", @cols);
			
			print $colNamesTabbed;
			print "\n";
		}
	
		# Link header for each column to column number
	
		@headerFields = split (/\t/, $_);
		
		for ($i = 0; $i <= $#headerFields; $i++)
		{
			$colNameToNum{$headerFields[$i]} = $i;
		}
	}
	
	if ($col_file)
	{
		checkChosenColsExist();
	}
	
	# Get this row (either this row is selected or all rows are requested)
	
	if (($fields[0] && $row_file && $rowsToGet{$fields[0]}) || ! $row_file)
	{
		if (! $col_file)	# Get all columns
		{
			if (! $row_file && $lineNum)
			{
				$id2info{$lineNum} = $_;
			}
			else
			{
				$id2info{$fields[0]} = $_;
			}
		}
		else
		{
			# Get the selected columns

			if ($row_file)	# Get only selected rows
			{
				$idKey = $fields[0];
			}
			else	# Get all rows -- and let redundant entries in first column
			{
				$idKey = $lineNum;
			}

			if ($row_file)	# start hash with key in first column
			{
				# Removed 7.6.07
				# $id2info{$idKey} = $fields[0];
			}
			else	# skip first column unless it's requested
			{
				$id2info{$idKey} = "";
			}

			foreach $selectedCol (@cols)
			{
				# Add a trailing tab if there's already something in the hash
				if ($id2info{$idKey})
				{
					$id2info{$idKey} .= "\t";
				}

				if ($colNameToNum{$selectedCol})
				{
					$id2info{$idKey} .= "$fields[$colNameToNum{$selectedCol}]";
				}
				elsif ($colNameToNum{$selectedCol} eq "0")
				{
					$id2info{$idKey} .= "$fields[$colNameToNum{$selectedCol}]";
				}
				else
				{
					$id2info{$idKey} .= "?";
				}
			}
		}
		
		if (! $row_file && $lineNum)	# Get all rows
		{
			push (@rows, $idKey);
			$rowsToGet{$fields[0]} = 1;
		}
	}
	
	$lineNum++;
}
close (EXP_DATA);

if ($row_file)
{
	foreach $inputGene (@rows)
	{
		if ($id2info{$inputGene})
		{
			print "$id2info{$inputGene}\n";
		}
		else
		{
			print "$inputGene\t?\n";
		}
	}
	close (SELECTED_OUT);
}
else	# Print all rows
{
	for ($i = 1; $i < $lineNum; $i++)
	{
		if ($id2info{$i})
		{
			print "$id2info{$i}\n";
		}
		else
		{
			print "$i\tNo data\n";
		}
	}
	close (SELECTED_OUT);
}

#####################  HTML out  #####################

print STDERR "All done.\n\n";

sub checkChosenColsExist
{
	foreach $colPresent (@headerFields)
	{
		$present{$colPresent} = 1;
	}
	foreach $colChosen (@cols)
	{
		if (! $present{$colChosen})
		{
			print STDERR "You've selected a column named \"$colChosen\"\n";
			print STDERR "that doesn't appear to be a column name in your input matrix file\n";
			exit;		
		}
	}
}

sub get_A_in_B
{
	# Added 8/17/2021 to remove columns in "to get" list but not present in matrix
	my ($inMatrix, $toGet) = @_;	
	my @inMatrix = @{ $inMatrix };
	my @toGet = @{ $toGet };
	
	%inMatrix = ();
	@toGetInMatrix = ();
	
	foreach $item (@inMatrix)
	{
		$inMatrix{$item} = 1;
	}

	foreach $item (@toGet)
	{
		if ($inMatrix{$item})
		{
			push @toGetInMatrix, $item;
			# print STDERR "Found column \"$item\" (present in matrix).\n";
		}
		else
		{
			print STDERR "Skipping column \"$item\" (not present in matrix).\n";
		}
	}
	return @toGetInMatrix;
}
