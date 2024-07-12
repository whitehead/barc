#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Brief Description: Parse a MAF (genome-genome) alignment into a fasta file
# This works for the MAF files associated with hg38
#
# Copyright(c) 2011 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 1.0 -- Created 8 December 2009
# Version: 1.1 -- 13 September 2011 Added argument for genome
# Version: 1.2 -- 16 June 2021 Updated for current 100-way multiz species
#
###########################################################################

# MAF format info: http://eyebrowse.cit.nih.gov/goldenPath/help/maf.html

if (! $ARGV[2] || ! -e $ARGV[0])
{
	print STDERR "\nParse a MAF alignment into a fasta file\n\n";
	print STDERR "Usage:   $0 mafFile genomeBuild[UCSC] geneName strand > out.fa\n";
	print STDERR "Example: $0 myGene.maf hg38 myGene + > myGene.fa\n\n";
	exit;
}

$maf = $ARGV[0];
$REF_GENOME = $ARGV[1];
$geneName = $ARGV[2];
$strand = $ARGV[3];

getSpeciesIDs();
$lineNum = 1;

# print STDERR "Reference genome = $REF_GENOME\n";

# Get all the species in the file
open (MAF, $maf) || die "Cannot open $maf for reading\n";
while (<MAF>)
{
	chomp;
	
	# As of Feb 2009 MAF, format, can be
	# s otoGar1.scaffold_103029.1-407137
	
	if (/^s ([^.]+)\.(\S+) /)	# sequence alignment
	{
		$species = $1;
		$isSpecies{$species} = 1;
	}
}
close (MAF);

open (MAF, $maf) || die "Cannot open $maf for reading\n";
while (<MAF>)
{
	chomp;
	
	if (/^s /)	# sequence alignment
	{
		# s hg18.chr8                39378038 134 + 146274826 CCACTTTTCAAAAACAACGGTTTCAACTCATTTTCTACATTTCTCTACCTGTGCTTCTTATTACTACTGCCA---TTTTAATAAAACAGAAAAAATTGAGAGAGCTATGCTATAGAGGGGAAAATGAAAGTGAGAGG
		# s panTro1.chr7             40832012 133 + 149542033 CTGTTCCTCAGAAA-AACGGCCTTTAATCAGTTGCTACATTTTTCTACCTTTCCTCATTTTAACTGCCATCA---TTGGTCTTAAATGGAACAAAATGAAGAGATTTTGGCACTGAGCAAAAAGTGTAAGTGGAAGG

		@f = split /\s+/, $_;
		
		$speciesPlusChr = $f[1];

		# if ($speciesPlusChr =~ /([^.]+)\.+(.+)/)
		if ($speciesPlusChr =~ /([^.]+)\.*(.*)/)
		{
			$species = $1;
			
			$isSpecies{$species} = 1;
		}
		else
		{
			print STDERR "Couldn't get species and chr from $speciesPlusChr ($_)\n";
		}
		
		# $startCoord = $f[2];
		# $length = $f[3];
		# $strand = $f[4];
		# $chrLength = $f[5];
		$seq = $f[6];
		
		# Added April 17, 2007
		$seq = uc($seq);

		# Get this block of sequence
		$species2seqBlock{$species} = $seq;
	}
	elsif (/^\s*$/)
	{
		# a score=83911.0
		# print STDERR "At line $. -- concatenating this block.\n";
		
		concatenateSeqs();
	}	
	
	$lineNum++;
}
close (MAF);

# Get the last one
if (%species2seqBlock)
{
	concatenateSeqs();
}

printLastSeqData();

###############################

sub getSpeciesIDs
{
	%assembly2species = 
	("hg38" => "Human",
	"oryAfe1" => "Aardvark",
	"vicPac2" => "Alpaca",
	"allMis1" => "American_alligator",
	"dasNov3" => "Armadillo",
	"gadMor1" => "Atlantic_cod",
	"papAnu2" => "Baboon",
	"camFer1" => "Bactrian_camel",
	"eptFus1" => "Big_brown_bat",
	"pteAle1" => "Black_flying-fox",
	"octDeg1" => "Brush-tailed_rat",
	"melUnd1" => "Budgerigar",
	"hapBur1" => "Burton's_mouthbreeder",
	"otoGar3" => "Bushbaby",
	"eleEdw1" => "Cape_elephant_shrew",
	"chrAsi1" => "Cape_golden_mole",
	"felCat8" => "Cat",
	"galGal4" => "Chicken",
	"panTro4" => "Chimp",
	"chiLan1" => "Chinchilla",
	"criGri1" => "Chinese_hamster",
	"pelSin1" => "Chinese_softshell_turtle",
	"tupChi1" => "Chinese_tree_shrew",
	"latCha1" => "Coelacanth",
	"ficAlb2" => "Collared_flycatcher",
	"bosTau8" => "Cow",
	"macFas5" => "Crab-eating_macaque",
	"myoDav1" => "David's_myotis_bat",
	"canFam3" => "Dog",
	"turTru2" => "Dolphin",
	"capHir1" => "Domestic_goat",
	"loxAfr3" => "Elephant",
	"musFur1" => "Ferret",
	"fr3" => "Fugu",
	"nomLeu3" => "Gibbon",
	"mesAur1" => "Golden_hamster",
	"gorGor3" => "Gorilla",
	"chlSab2" => "Green_monkey",
	"cheMyd1" => "Green_seaturtle",
	"cavPor3" => "Guinea_pig",
	"eriEur2" => "Hedgehog",
	"equCab3" => "Horse",
	"orcOrc1" => "Killer_whale",
	"petMar2" => "Lamprey",
	"jacJac1" => "Lesser_Egyptian_jerboa",
	"anoCar2" => "Lizard",
	"anaPla1" => "Mallard_duck",
	"triMan1" => "Manatee",
	"calJac3" => "Marmoset",
	"oryLat2" => "Medaka",
	"geoFor1" => "Medium_ground_finch",
	"pteVam1" => "Megabat",
	"astMex1" => "Mexican_tetra_(cavefish)",
	"myoLuc2" => "Microbat",
	"mm10" => "Mouse",
	"hetGla2" => "Naked_mole-rat",
	"oreNil2" => "Nile_tilapia",
	"monDom5" => "Opossum",
	"ponAbe2" => "Orangutan",
	"odoRosDiv1" => "Pacific_walrus",
	"chrPic2" => "Painted_turtle",
	"ailMel1" => "Panda",
	"amaVit1" => "Parrot",
	"falPer1" => "Peregrine_falcon",
	"susScr3" => "Pig",
	"ochPri3" => "Pika",
	"ornAna1" => "Platypus",
	"micOch1" => "Prairie_vole",
	"neoBri1" => "Princess_of_Burundi",
	"punNye1" => "Pundamilia_nyererei",
	"oryCun2" => "Rabbit",
	"rn6" => "Rat",
	"rheMac3" => "Rhesus",
	"colLiv1" => "Rock_pigeon",
	"falChe1" => "Saker_falcon",
	"araMac1" => "Scarlet_macaw",
	"oviAri3" => "Sheep",
	"sorAra2" => "Shrew",
	"xipMac1" => "Southern_platyfish",
	"apaSpi1" => "Spiny_softshell_turtle",
	"lepOcu1" => "Spotted_gar",
	"speTri2" => "Squirrel",
	"saiBol1" => "Squirrel_monkey",
	"conCri1" => "Star-nosed_mole",
	"gasAcu1" => "Stickleback",
	"sarHar1" => "Tasmanian_devil",
	"echTel2" => "Tenrec",
	"tetNig2" => "Tetraodon",
	"panHod1" => "Tibetan_antelope",
	"pseHum1" => "Tibetan_ground_jay",
	"melGal1" => "Turkey",
	"macEug2" => "Wallaby",
	"lepWed1" => "Weddell_seal",
	"cerSim1" => "White_rhinoceros",
	"zonAlb1" => "White-throated_sparrow",
	"xenTro7" => "X._tropicalis",
	"takFla1" => "Yellowbelly_pufferfish",
	"taeGut2" => "Zebra_finch",
	"mayZeb1" => "Zebra_mbuna",
	"danRer10" => "Zebrafish");	

	%species2assembly = 
	("Human" => "hg38",
	"Aardvark" => "oryAfe1",
	"Alpaca" => "vicPac2",
	"American_alligator" => "allMis1",
	"Armadillo" => "dasNov3",
	"Atlantic_cod" => "gadMor1",
	"Baboon" => "papAnu2",
	"Bactrian_camel" => "camFer1",
	"Big_brown_bat" => "eptFus1",
	"Black_flying-fox" => "pteAle1",
	"Brush-tailed_rat" => "octDeg1",
	"Budgerigar" => "melUnd1",
	"Burton's_mouthbreeder" => "hapBur1",
	"Bushbaby" => "otoGar3",
	"Cape_elephant_shrew" => "eleEdw1",
	"Cape_golden_mole" => "chrAsi1",
	"Cat" => "felCat8",
	"Chicken" => "galGal4",
	"Chimp" => "panTro4",
	"Chinchilla" => "chiLan1",
	"Chinese_hamster" => "criGri1",
	"Chinese_softshell_turtle" => "pelSin1",
	"Chinese_tree_shrew" => "tupChi1",
	"Coelacanth" => "latCha1",
	"Collared_flycatcher" => "ficAlb2",
	"Cow" => "bosTau8",
	"Crab-eating_macaque" => "macFas5",
	"David's_myotis_bat" => "myoDav1",
	"Dog" => "canFam3",
	"Dolphin" => "turTru2",
	"Domestic_goat" => "capHir1",
	"Elephant" => "loxAfr3",
	"Ferret" => "musFur1",
	"Fugu" => "fr3",
	"Gibbon" => "nomLeu3",
	"Golden_hamster" => "mesAur1",
	"Gorilla" => "gorGor3",
	"Green_monkey" => "chlSab2",
	"Green_seaturtle" => "cheMyd1",
	"Guinea_pig" => "cavPor3",
	"Hedgehog" => "eriEur2",
	"Horse" => "equCab3",
	"Killer_whale" => "orcOrc1",
	"Lamprey" => "petMar2",
	"Lesser_Egyptian_jerboa" => "jacJac1",
	"Lizard" => "anoCar2",
	"Mallard_duck" => "anaPla1",
	"Manatee" => "triMan1",
	"Marmoset" => "calJac3",
	"Medaka" => "oryLat2",
	"Medium_ground_finch" => "geoFor1",
	"Megabat" => "pteVam1",
	"Mexican_tetra_(cavefish)" => "astMex1",
	"Microbat" => "myoLuc2",
	"Mouse" => "mm10",
	"Naked_mole-rat" => "hetGla2",
	"Nile_tilapia" => "oreNil2",
	"Opossum" => "monDom5",
	"Orangutan" => "ponAbe2",
	"Pacific_walrus" => "odoRosDiv1",
	"Painted_turtle" => "chrPic2",
	"Panda" => "ailMel1",
	"Parrot" => "amaVit1",
	"Peregrine_falcon" => "falPer1",
	"Pig" => "susScr3",
	"Pika" => "ochPri3",
	"Platypus" => "ornAna1",
	"Prairie_vole" => "micOch1",
	"Princess_of_Burundi" => "neoBri1",
	"Pundamilia_nyererei" => "punNye1",
	"Rabbit" => "oryCun2",
	"Rat" => "rn6",
	"Rhesus" => "rheMac3",
	"Rock_pigeon" => "colLiv1",
	"Saker_falcon" => "falChe1",
	"Scarlet_macaw" => "araMac1",
	"Sheep" => "oviAri3",
	"Shrew" => "sorAra2",
	"Southern_platyfish" => "xipMac1",
	"Spiny_softshell_turtle" => "apaSpi1",
	"Spotted_gar" => "lepOcu1",
	"Squirrel" => "speTri2",
	"Squirrel_monkey" => "saiBol1",
	"Star-nosed_mole" => "conCri1",
	"Stickleback" => "gasAcu1",
	"Tasmanian_devil" => "sarHar1",
	"Tenrec" => "echTel2",
	"Tetraodon" => "tetNig2",
	"Tibetan_antelope" => "panHod1",
	"Tibetan_ground_jay" => "pseHum1",
	"Turkey" => "melGal1",
	"Wallaby" => "macEug2",
	"Weddell_seal" => "lepWed1",
	"White_rhinoceros" => "cerSim1",
	"White-throated_sparrow" => "zonAlb1",
	"X._tropicalis" => "xenTro7",
	"Yellowbelly_pufferfish" => "takFla1",
	"Zebra_finch" => "taeGut2",
	"Zebra_mbuna" => "mayZeb1",
	"Zebrafish" => "danRer10");
	
	@expectedSpecies = qw(Aardvark Alpaca American_alligator Armadillo Atlantic_cod Baboon Bactrian_camel Big_brown_bat Black_flying-fox Brush-tailed_rat Budgerigar Burton's_mouthbreeder Bushbaby Cape_elephant_shrew Cape_golden_mole Cat Chicken Chimp Chinchilla Chinese_hamster Chinese_softshell_turtle Chinese_tree_shrew Coelacanth Collared_flycatcher Cow Crab-eating_macaque David's_myotis_bat Dog Dolphin Domestic_goat Elephant Ferret Fugu Gibbon Golden_hamster Gorilla Green_monkey Green_seaturtle Guinea_pig Hedgehog Horse Human Killer_whale Lamprey Lesser_Egyptian_jerboa Lizard Mallard_duck Manatee Marmoset Medaka Medium_ground_finch Megabat Mexican_tetra_(cavefish) Microbat Mouse Naked_mole-rat Nile_tilapia Opossum Orangutan Pacific_walrus Painted_turtle Panda Parrot Peregrine_falcon Pig Pika Platypus Prairie_vole Princess_of_Burundi Pundamilia_nyererei Rabbit Rat Rhesus Rock_pigeon Saker_falcon Scarlet_macaw Sheep Shrew Southern_platyfish Spiny_softshell_turtle Spotted_gar Squirrel Squirrel_monkey Star-nosed_mole Stickleback Tasmanian_devil Tenrec Tetraodon Tibetan_antelope Tibetan_ground_jay Turkey Wallaby Weddell_seal White_rhinoceros White-throated_sparrow X._tropicalis Yellowbelly_pufferfish Zebra_finch Zebra_mbuna Zebrafish);
}

sub printLastSeqData
{
	if ($strand eq "-")		# Correct 3' UTR be removing first nt [as of 3.25.07]
	{
		# reverse complement each sequence
		foreach $species (sort keys %isSpecies)
		{
			$rev_comp = $species2seq{$species};
			$rev_comp = reverse $rev_comp;
			$rev_comp =~ tr/ACGT/TGCA/;
			# Added April 17, 2007
			$rev_comp =~ tr/acgt/tgca/;
			
			$species2seq{$species} = $rev_comp;
			
			if ($assembly2species{$species})
			{
				$speciesNew = $assembly2species{$species};
				$species2seq{$speciesNew} = $species2seq{$species};
			}
		}
	}
	else
	{
		foreach $species (sort keys %isSpecies)
		{
			if ($assembly2species{$species})
			{
				$speciesNew = $assembly2species{$species};
				$species2seq{$speciesNew} = $species2seq{$species};
			}
		}
	}

	# Go through expected species
	foreach $species (@expectedSpecies)
	{
		if ($species2seq{$species})
		{
			$thisSeq = $species2seq{$species};
	
			print ">${geneName}_${species}\n$thisSeq\n";
			
			$isSpecies{$species} = 0;	# Mark that we already printed this
			if ($species2assembly{$species})
			{
				$isSpecies{$species2assembly{$species}} = 0;	# Mark that we already printed this
			}
		}
	}
	
	# Go through any user-selected species not in the expected set
	foreach $species (sort keys %isSpecies)
	{
		if ($isSpecies{$species})	# Didn't already get this
		{
			$thisSeq = $species2seq{$species};
	
			print ">${geneName}_${species}\n$thisSeq\n";
		}
	}
	
	# @expectedSpecies
}

sub concatenateSeqs
{
	$spacerSeq = $species2seqBlock{$REF_GENOME};
	# print STDERR "spacerSeq => $spacerSeq\n";
	$spacerSeq =~ s/[A-Za-z]/-/g;
	
	foreach $species (sort keys %isSpecies)
	{
		if ($species2seqBlock{$species})
		{
			$species2seq{$species} .= $species2seqBlock{$species};
		}
		else
		{
			if ($species2seq{$species})
			{
				$species2seq{$species} .= $spacerSeq;
			}
			else
			{
				$species2seq{$species} = $spacerSeq;
			}
		}
	}
	
	# Reset variable
	%species2seqBlock = ();
}
