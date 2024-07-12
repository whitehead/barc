#! /usr/bin/env perl
use warnings;

##################################################################
#
# Copyright(c) 2012 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#
##################################################################

################  User-supplied variables  #############

if (! $ARGV[0])
{
	print STDERR "\nConvert a directory of fasta (.fa) files into nib (BLAT index) files\n";
	print STDERR "USAGE: $0 fastaDir\n\n";
	exit;
}

# Directory of sequences
$fastaDir = $ARGV[0];

# Put nib files into 'blat' directory
$blatDir = "blat";

if (! -e $blatDir)
{
	mkdir $blatDir;
	print STDERR "Created directory \"$blatDir\" for nib files\n";
}

#########################################################

# Go to sequence directory and open it (i.e, read contents)
# chdir($fastaDir) || die "Cannot change to $fastaDir: $!";      # Go to $myDir
opendir(DIR, "$fastaDir") || die "Cannot open $fastaDir: $!";      # Open $myDir

foreach $seqFile (sort readdir(DIR))
{
	if ($seqFile =~ /\.fa$/)      # if file ends in .fa
	{
		print "Processing $seqFile\n";
		
		$blatFile = $seqFile;
		$blatFile =~ s/\.fa$/.nib/;

		#################  Run faToNib  #################

		# print "/usr/bin/faToNib -softMask $fastaDir/$seqFile $blatDir/$blatFile";
		`/nfs/BaRC_Public/BaRC_code/Perl/make_nib_files/faToNib -softMask $fastaDir/$seqFile $blatDir/$blatFile`;

	}
}
