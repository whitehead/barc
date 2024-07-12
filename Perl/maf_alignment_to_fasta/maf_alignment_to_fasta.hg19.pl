#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Brief Description: Parse a MAF (genome-genome) alignment into a fasta file
# This works for the MAF files associated with hg19
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
	("hg19" => "Human",
	"papHam1" => "Baboon",
	"otoGar3" => "Bushbaby",
	"panTro4" => "Chimp",
	"macFas5" => "Crab-eating_macaque",
	"nomLeu3" => "Gibbon",
	"gorGor3" => "Gorilla",
	"chlSab1" => "Green_monkey",
	"calJac3" => "Marmoset",
	"ponAbe2" => "Orangutan",
	"rheMac3" => "Rhesus",
	"saiBol1" => "Squirrel_monkey",
	"octDeg1" => "Brush-tailed_rat",
	"chiLan1" => "Chinchilla",
	"criGri1" => "Chinese_hamster",
	"tupChi1" => "Chinese_tree_shrew",
	"mesAur1" => "Golden_hamster",
	"cavPor3" => "Guinea_pig",
	"jacJac1" => "Lesser_Egyptian_jerboa",
	"mm10" => "Mouse",
	"hetGla2" => "Naked_mole-rat",
	"ochPri3" => "Pika",
	"micOch1" => "Prairie_vole",
	"oryCun2" => "Rabbit",
	"rn5" => "Rat",
	"speTri2" => "Squirrel",
	"vicPac2" => "Alpaca",
	"camFer1" => "Bactrian_camel",
	"eptFus1" => "Big_brown_bat",
	"pteAle1" => "Black_flying-fox",
	"felCat5" => "Cat",
	"bosTau7" => "Cow",
	"myoDav1" => "David's_myotis_bat",
	"canFam3" => "Dog",
	"turTru2" => "Dolphin",
	"capHir1" => "Domestic_goat",
	"musFur1" => "Ferret",
	"eriEur2" => "Hedgehog",
	"equCab2" => "Horse",
	"orcOrc1" => "Killer_whale",
	"pteVam1" => "Megabat",
	"myoLuc2" => "Microbat",
	"odoRosDiv1" => "Pacific_walrus",
	"ailMel1" => "Panda",
	"susScr3" => "Pig",
	"oviAri3" => "Sheep",
	"sorAra2" => "Shrew",
	"conCri1" => "Star-nosed_mole",
	"panHod1" => "Tibetan_antelope",
	"lepWed1" => "Weddell_seal",
	"cerSim1" => "White_rhinoceros",
	"oryAfe1" => "Aardvark",
	"eleEdw1" => "Cape_elephant_shrew",
	"chrAsi1" => "Cape_golden_mole",
	"loxAfr3" => "Elephant",
	"triMan1" => "Manatee",
	"echTel2" => "Tenrec",
	"dasNov3" => "Armadillo",
	"monDom5" => "Opossum",
	"ornAna1" => "Platypus",
	"sarHar1" => "Tasmanian_devil",
	"macEug2" => "Wallaby",
	"melUnd1" => "Budgerigar",
	"galGal4" => "Chicken",
	"ficAlb2" => "Collared_flycatcher",
	"anaPla1" => "Mallard_duck",
	"geoFor1" => "Medium_ground_finch",
	"amaVit1" => "Parrot",
	"falPer1" => "Peregrine_falcon",
	"colLiv1" => "Rock_pigeon",
	"falChe1" => "Saker_falcon",
	"araMac1" => "Scarlet_macaw",
	"pseHum1" => "Tibetan_ground_jay",
	"melGal1" => "Turkey",
	"zonAlb1" => "White-throated_sparrow",
	"taeGut2" => "Zebra_finch",
	"allMis1" => "American_alligator",
	"pelSin1" => "Chinese_softshell_turtle",
	"latCha1" => "Coelacanth",
	"cheMyd1" => "Green_seaturtle",
	"anoCar2" => "Lizard",
	"chrPic1" => "Painted_turtle",
	"apaSpi1" => "Spiny_softshell_turtle",
	"xenTro7" => "X._tropicalis",
	"gadMor1" => "Atlantic_cod",
	"hapBur1" => "Burton's_mouthbreeder",
	"fr3" => "Fugu",
	"petMar2" => "Lamprey",
	"oryLat2" => "Medaka",
	"astMex1" => "Mexican_tetra_(cavefish)",
	"oreNil2" => "Nile_tilapia",
	"neoBri1" => "Princess_of_Burundi",
	"punNye1" => "Pundamilia_nyererei",
	"xipMac1" => "Southern_platyfish",
	"lepOcu1" => "Spotted_gar",
	"gasAcu1" => "Stickleback",
	"tetNig2" => "Tetraodon",
	"takFla1" => "Yellowbelly_pufferfish",
	"mayZeb1" => "Zebra_mbuna",
	"danRer7" => "Zebrafish");	

	%species2assembly = 
	("Human" => "hg19",
	"Baboon" => "papHam1",
	"Bushbaby" => "otoGar3",
	"Chimp" => "panTro4",
	"Crab-eating_macaque" => "macFas5",
	"Gibbon" => "nomLeu3",
	"Gorilla" => "gorGor3",
	"Green_monkey" => "chlSab1",
	"Marmoset" => "calJac3",
	"Orangutan" => "ponAbe2",
	"Rhesus" => "rheMac3",
	"Squirrel_monkey" => "saiBol1",
	"Brush-tailed_rat" => "octDeg1",
	"Chinchilla" => "chiLan1",
	"Chinese_hamster" => "criGri1",
	"Chinese_tree_shrew" => "tupChi1",
	"Golden_hamster" => "mesAur1",
	"Guinea_pig" => "cavPor3",
	"Lesser_Egyptian_jerboa" => "jacJac1",
	"Mouse" => "mm10",
	"Naked_mole-rat" => "hetGla2",
	"Pika" => "ochPri3",
	"Prairie_vole" => "micOch1",
	"Rabbit" => "oryCun2",
	"Rat" => "rn5",
	"Squirrel" => "speTri2",
	"Alpaca" => "vicPac2",
	"Bactrian_camel" => "camFer1",
	"Big_brown_bat" => "eptFus1",
	"Black_flying-fox" => "pteAle1",
	"Cat" => "felCat5",
	"Cow" => "bosTau7",
	"David's_myotis_bat" => "myoDav1",
	"Dog" => "canFam3",
	"Dolphin" => "turTru2",
	"Domestic_goat" => "capHir1",
	"Ferret" => "musFur1",
	"Hedgehog" => "eriEur2",
	"Horse" => "equCab2",
	"Killer_whale" => "orcOrc1",
	"Megabat" => "pteVam1",
	"Microbat" => "myoLuc2",
	"Pacific_walrus" => "odoRosDiv1",
	"Panda" => "ailMel1",
	"Pig" => "susScr3",
	"Sheep" => "oviAri3",
	"Shrew" => "sorAra2",
	"Star-nosed_mole" => "conCri1",
	"Tibetan_antelope" => "panHod1",
	"Weddell_seal" => "lepWed1",
	"White_rhinoceros" => "cerSim1",
	"Aardvark" => "oryAfe1",
	"Cape_elephant_shrew" => "eleEdw1",
	"Cape_golden_mole" => "chrAsi1",
	"Elephant" => "loxAfr3",
	"Manatee" => "triMan1",
	"Tenrec" => "echTel2",
	"Armadillo" => "dasNov3",
	"Opossum" => "monDom5",
	"Platypus" => "ornAna1",
	"Tasmanian_devil" => "sarHar1",
	"Wallaby" => "macEug2",
	"Budgerigar" => "melUnd1",
	"Chicken" => "galGal4",
	"Collared_flycatcher" => "ficAlb2",
	"Mallard_duck" => "anaPla1",
	"Medium_ground_finch" => "geoFor1",
	"Parrot" => "amaVit1",
	"Peregrine_falcon" => "falPer1",
	"Rock_pigeon" => "colLiv1",
	"Saker_falcon" => "falChe1",
	"Scarlet_macaw" => "araMac1",
	"Tibetan_ground_jay" => "pseHum1",
	"Turkey" => "melGal1",
	"White-throated_sparrow" => "zonAlb1",
	"Zebra_finch" => "taeGut2",
	"American_alligator" => "allMis1",
	"Chinese_softshell_turtle" => "pelSin1",
	"Coelacanth" => "latCha1",
	"Green_seaturtle" => "cheMyd1",
	"Lizard" => "anoCar2",
	"Painted_turtle" => "chrPic1",
	"Spiny_softshell_turtle" => "apaSpi1",
	"X._tropicalis" => "xenTro7",
	"Atlantic_cod" => "gadMor1",
	"Burton's_mouthbreeder" => "hapBur1",
	"Fugu" => "fr3",
	"Lamprey" => "petMar2",
	"Medaka" => "oryLat2",
	"Mexican_tetra_(cavefish)" => "astMex1",
	"Nile_tilapia" => "oreNil2",
	"Princess_of_Burundi" => "neoBri1",
	"Pundamilia_nyererei" => "punNye1",
	"Southern_platyfish" => "xipMac1",
	"Spotted_gar" => "lepOcu1",
	"Stickleback" => "gasAcu1",
	"Tetraodon" => "tetNig2",
	"Yellowbelly_pufferfish" => "takFla1",
	"Zebra_mbuna" => "mayZeb1",
	"Zebrafish" => "danRer7");
	
	@expectedSpecies = qw(Baboon Bushbaby Chimp Crab-eating_macaque Gibbon Gorilla Green_monkey Human Marmoset Orangutan Rhesus Squirrel_monkey Brush-tailed_rat Chinchilla Chinese_hamster Chinese_tree_shrew Golden_hamster Guinea_pig Lesser_Egyptian_jerboa Mouse Naked_mole-rat Pika Prairie_vole Rabbit Rat Squirrel Alpaca Bactrian_camel Big_brown_bat Black_flying-fox Cat Cow David's_myotis_bat Dog Dolphin Domestic_goat Ferret Hedgehog Horse Killer_whale Megabat Microbat Pacific_walrus Panda Pig Sheep Shrew Star-nosed_mole Tibetan_antelope Weddell_seal White_rhinoceros Aardvark Cape_elephant_shrew Cape_golden_mole Elephant Manatee Tenrec Armadillo Opossum Platypus Tasmanian_devil Wallaby Budgerigar Chicken Collared_flycatcher Mallard_duck Medium_ground_finch Parrot Peregrine_falcon Rock_pigeon Saker_falcon Scarlet_macaw Tibetan_ground_jay Turkey White-throated_sparrow Zebra_finch American_alligator Chinese_softshell_turtle Coelacanth Green_seaturtle Lizard Painted_turtle Spiny_softshell_turtle X._tropicalis Atlantic_cod Burton's_mouthbreeder Fugu Lamprey Medaka Mexican_tetra_(cavefish) Nile_tilapia Princess_of_Burundi Pundamilia_nyererei Southern_platyfish Spotted_gar Stickleback Tetraodon Yellowbelly_pufferfish Zebra_mbuna Zebrafish);
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
