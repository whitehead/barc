#!/bin/bash

###
###  Create over.chain file from two closely related genomes
###    using the tools and method from UCSC Bioinformatics
###
###  Sample command: ./create_over.chain_for_liftOver.sh genomes/Ler.fa genomes/Col.fa 
###
###  Pipeline from https://hgwdev.gi.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt
###                https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/hg/doc/liftOver.txt
###                http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver
###
###  George Bell, BaRC
###  Version 1.0:  6 May 2022
###

if [  $# -le 3 ] 
	then
		thisScript=$0	# Get the name and path of this script
		printf "\nCreate an over.chain file for liftOver using two closely related genomes,\n"
		printf "  given two whole-genome fasta files.\n"
		printf "USAGE: $thisScript oldGenome newGenome chunkLength portNum\n\n"
		printf "Ex:    $thisScript oldHs.fa newHs.fa 3000 1234\n\n"
		exit 1
	fi

# This is the original genome (with the "old" coordinates that we want to start with)
oldGenomeFasta=$1

# This is the new genome (with the "new" coordinates that we want to convert to)
newGenomeFasta=$2

# When we split each chromosome, how long should each piece be?  Recommended = 3000
chunkLength=$3

# gfServer portNum
portNum=$4

####################################################################################################

# Drop .fa or .fasta from old genome
oldGenome=`basename $oldGenomeFasta | perl -pe 's/\.fa//; s/\.fasta//'`

# Drop .fa or .fasta from old genome
newGenome=`basename $newGenomeFasta | perl -pe 's/\.fa//; s/\.fasta//'`

# Create directory for all files and go there
oldToNew="${oldGenome}_to_${newGenome}"
projectDir="${oldToNew}_liftOver"

>&2 echo Creating $projectDir ....
mkdir $projectDir

>&2 echo Creating links to genomes in $projectDir and entering $projectDir ....
ln -s ../$oldGenomeFasta $projectDir/.
ln -s ../$newGenomeFasta $projectDir/.
cd $projectDir 

>&2 echo Converting genome files to 2bit format and making chrom.sizes files ....
oldGenomeFasta=`basename $oldGenomeFasta`
faToTwoBit $oldGenomeFasta $oldGenome.2bit
twoBitInfo $oldGenome.2bit $oldGenome.chrom.sizes

newGenomeFasta=`basename $newGenomeFasta`
faToTwoBit $newGenomeFasta $newGenome.2bit
twoBitInfo $newGenome.2bit $newGenome.chrom.sizes

# Split new genome into chunks
>&2 echo Splitting new genome into chunks of $chunkLength nt  ....
faSplit -lift=$newGenome.lft size $newGenome.fa -oneFile $chunkLength $newGenome.split

# Use gfServer and gfClient for BLAT
>&2 echo Starting gfServer on $HOSTNAME to get $oldGenome in memory ....
gfServer start $HOSTNAME $portNum $oldGenome.2bit &
>&2 echo Pausing for 1 minute to allow gfServer to complete indexing
sleep 1m
>&2 echo Continuing ....
>&2 echo Running gfClient to align $newGenome pieces to $oldGenome ....
splitPsl="${oldToNew}.split.psl"
# This (gfClient) is by far the longest step
# It could be sped up by splitting ther query file and running multiple concurrent jobs
#   but that would require keeping track of when all of the jobs are complete.
gfClient $HOSTNAME $portNum . $newGenome.split.fa $splitPsl

# Make backup of split psl file (just in case)
cp $splitPsl ${splitPsl}.COPY

>&2 echo Turning off gfServer ....
>&2 echo gfServer stop $HOSTNAME $portNum
gfServer stop $HOSTNAME $portNum

# Modify coordinates in psl files to match the original chr coordinates
psl="${oldToNew}.psl"
liftUp -pslQ $psl $newGenome.lft warn $splitPsl

>&2 echo Chaining and netting alignments ....
# Chain together axt alignments
axtChain -linearGap=medium -psl $psl $oldGenome.2bit $newGenome.2bit ${oldToNew}.chain

# Sort chains
chainMergeSort *.chain | chainSplit chain stdin

# Make alignment nets out of chains
chainNet ${oldToNew}.chain ${oldGenome}.chrom.sizes ${newGenome}.chrom.sizes chr.net /dev/null

# Create chain file with subset of chains that appear in the net
netChainSubset chr.net ${oldToNew}.chain ${oldToNew}.over.chain

>&2 printf "\nAll done!  See $projectDir/${oldToNew}.over.chain\n\n"
