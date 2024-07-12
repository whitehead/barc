#!/usr/bin/env perl

##################################################################
#
# Copyright(c) 2021 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#
# Convert the "attributes" of a gtf or gtf.gz file into tab-delimited text
#
# Version 1.0: 21 December 2021
# Version 1.1: 2 December 2022 -- Add option to choose feature type to extract
#
##################################################################

use warnings;

if (! $ARGV[0])
{
	print STDERR "\nExtract all gene/transcript annotations from field 9 of a GTF/GTF.GZ file\n";
	print STDERR "USAGE: $0 features.gtf [featureType] > features.txt\n\n";
	print STDERR "Ex:    $0 features.gtf gene > features.genes.txt\n\n";
	exit;
}

# Input file can be gtf or gtf.gz
# May work with some GFF files of a GTF style

$gtf_gff = $ARGV[0];

# Only look at rows for this feature
if ($ARGV[1])
{
	$FEATURE_TO_EXAMINE = $ARGV[1];
	print STDERR "\nReading all lines describing a \"$FEATURE_TO_EXAMINE\" feature ....\n\n";
} 
else
{
	$FEATURE_TO_EXAMINE = ".";
}

if ($gtf_gff =~ /\.gz/i)
{
	open (GTF_GFF, "zmore $gtf_gff |") || die "Cannot open $gtf_gff for reading: $!";
}
else
{
	open (GTF_GFF, "$gtf_gff") || die "Cannot open $gtf_gff for reading: $!";
} 
while (<GTF_GFF>)
{
	chomp;
	@f = split /\t/, $_;
	
	# Only look at $FEATURE_TO_EXAMINE features
	if (! /^#/)
	{
		if ($f[2] =~ $FEATURE_TO_EXAMINE)
		{
			$features = $f[8];

			if ($features =~ /gene_id (\S+)/)
			{
				$geneID = $1;
				$geneID =~ s/\"//g;
				$geneID =~ s/;$//g; 
			}

			if ($features =~ /transcript_id (\S+)/)
			{
				$transcriptID = $1;
				$transcriptID =~ s/\"//g;
				$transcriptID =~ s/;$//g; 

				$transcript2info{$transcriptID}{"gene_id"} = $geneID;
			}
			else
			{
				# print STDERR "Cannot find gene ID in row $.\n";
				next;
			}


			@features = split (/; /, $features);

			foreach $feature (@features)
			{
				if ($feature =~ /(\S+) (.+)/)
				{
					$tag = $1;
					$info = $2;
					$info =~ s/\"//g;

					if ($tag ne "transcript_id")
					{
						$gotTag{$tag}++;
						# print "$transcriptID\t$geneID\t$tag\t$info\n";

						if (! $transcript2info{$transcriptID}{$tag})
						{
							$transcript2info{$transcriptID}{$tag} = $info;
						} 
						else	# If we've seen this tag for this transcript before, add to the list of tags
						{
							$transcript2info{$transcriptID}{$tag} .= "; $info";
						}
					}
				}
				else
				{
					print STDERR "Cannot parse $feature\n";
				}
			}
		}
	}
}

# Sort tags by number of times we saw them
@tags = sort { $gotTag{$b} <=> $gotTag{$a} } keys(%gotTag);

for $tag (@tags) 
{
	print STDERR "Saw $tag $gotTag{$tag} times.\n";
}
print STDERR "\n";

print "Transcript.ID";
for $tag (@tags) 
{
	print "\t$tag";
}
print "\n";

for $transcriptID (sort keys %transcript2info ) 
{ 
	print "$transcriptID"; 
	for $tag (@tags) 
	{
		if ($transcript2info{$transcriptID}{$tag})
		{
			if ($transcript2info{$transcriptID}{$tag} =~ /; /)
			{
				$infoList = $transcript2info{$transcriptID}{$tag};
				@info = split /; /, $infoList;
				# Make this list of items unique
				@info = sort (do { my %seen; grep { !$seen{$_}++ } @info });
				$infoList = join "; ", @info;
			}
			else
			{
				$infoList = $transcript2info{$transcriptID}{$tag}
			}
		
			print "\t$infoList";
		}
		else
		{
			print "\t";
		}
	} 
	print "\n"; 
} 
