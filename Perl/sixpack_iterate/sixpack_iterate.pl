#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2010 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Translate each gene sequence in a multiple sequence file (using EMBOSS's sixpack)
# and save only the longest ORF for each input sequence,
# printing output to STDOUT
#
# George Bell, BaRC -- 22 January 2010
#
########################################################################


use Bio::SeqIO;

if (! $ARGV[0])
{
	print STDERR "\nTranslate a multiple sequence file in six frames,\n";
	print STDERR "  keeping the longest ORF for each gene sequence\n";
	print STDERR "  and printing to STDOUT\n";
	print STDERR "\nUSAGE: $0 inputGenesFile > outputProteinsFile\n\n";
	exit;
}

$inFile = $ARGV[0];

$in  = Bio::SeqIO->new('-file' => "$inFile" , '-format' => 'Fasta');

while ($seqobj = $in->next_seq())
{

	$seqID = $seqobj->display_id();
	$seqID =~ s/gnl\|exprUnit\|//;

	print STDERR "Processing $seqID\n";

	# Print sequence in output format
	$out = Bio::SeqIO->new('-format' => 'Fasta', '-file' => ">$seqID.fa");
	$out->write_seq($seqobj);

	# Run sixpack to translate sequence	
	# $sixpack_command = "sixpack -sequence $seqID.fa -orfminsize 30 -outseq $seqID.sixpack.fa -outfile sixpack.temp -auto";	
	$sixpack_command = "sixpack -sequence $seqID.fa -outseq $seqID.sixpack.fa -outfile sixpack.temp -auto";	
	`$sixpack_command`;
	
	$longestORFlength = 0;
	$longestORFid = $longestORF = "";
	
	###  Did sixpack produce output?
	###  If so, only keep the longest ORF
	
	if (-e "$seqID.sixpack.fa")
	{
		$inProtein  = Bio::SeqIO->new('-file' => "$seqID.sixpack.fa" , '-format' => 'Fasta');
		while ($seqobjProtein = $inProtein->next_seq())
		{
			# Get the ID (first word)
			$seqID_protein = $seqobjProtein->display_id();
			
			# Get the description (the rest of the header line)
			$description_protein = $seqobjProtein->description();
			# Clean up the description [optional]
			$description_protein =~ s/threshold 1, //;
			$description_protein =~ s/of SMED_\d+_V2 //;
			
			# Get the protein length from the header
			$proteinLength = $description_protein;
			$proteinLength =~ s/.+, (\d+)aa/$1/;
			
			# Is this ORF the longest one so far?
			if ($proteinLength > $longestORFlength)
			{
				$longestORFid = $seqID_protein;
				$longestORFdescription = $description_protein;
				$longestORFlength = $proteinLength;
				$longestORF = $seqobjProtein->seq();
			}
		}
		
		# Print the longest ORF for this contig
		print ">$seqID $longestORFdescription\n$longestORF\n";
		
		# Delete output file for this contig
		unlink("$seqID.sixpack.fa");
	}
	# Delete inut file for this contig
	unlink("$seqID.fa");
}

# Delete the shared temp file
unlink("sixpack.temp");
