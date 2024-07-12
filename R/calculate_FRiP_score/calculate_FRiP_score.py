#!/usr/bin/env python3

##################################################################
#
# Copyright(c) 2020 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
# Author: Austin Laramee
#         Bioinformatics and Research Computing
#         barc@wi.mit.edu
#
# Version: 1.0 -- June 18, 2020
#
# Using a BAM file and a matching peaks file (narrowPeak or bed format),
# calculate the FRiP (fraction of reads in peaks) score.
#
# For interpretation of this score, see
#   Landt et al. https://pubmed.ncbi.nlm.nih.gov/22955991/
#
##################################################################

import sys
import os
import os.path
from os import path
# Ensure that the arguments were entered
try:
    BAMFullFile = sys.argv[1]
    NPFullFile = sys.argv[2]
except IndexError:
    sys.stderr.write("\nCalculate the FRiP score from a BAM file and a narrowPeak or bed file.\n\n")
    sys.stderr.write("USAGE: calculate_FRiP_score.py Path/To/File.bam Path/To/File.narrowPeak")
    sys.stderr.write("\n\nTo fractionally count multimapped reads:\n")
    sys.stderr.write("       calculate_FRiP_score.py Path/To/File.bam Path/To/File.narrowPeak -M\n\n")
    exit()

sys.stderr.write("\nCalculating the FRiP score from a .bam and a .narrowPeak file ...\n")

BAMLocations = os.path.split(os.path.abspath(BAMFullFile))
BAMPathName = str(BAMLocations[0])
BAMFileName = str(BAMLocations[1])
# Verify that the arguments include 1 bam file and 1 narrowPeak file
if path.exists(BAMFullFile):
    if ".bam" == BAMFullFile[len(BAMFullFile)-4:len(BAMFullFile)]:
        pass
    else:
        sys.stderr.write("The file you have selected, '" + BAMFullFile + "', is not a '.bam' file. This script will now terminate without calculating the FRiP score")
        exit()
else:
    sys.stderr.write("THE FILE '" + str(BAMFullFile) + "' CANNOT BE LOCATED. This script will now terminate without calculating the FRiP score")
    exit()
NPLocations = os.path.split(os.path.abspath(NPFullFile))
NPPathName = str(NPLocations[0])
NPFileName = str(NPLocations[1])
if path.exists(NPFullFile):
    if ".narrowPeak" == NPFullFile[len(NPFullFile)-11:len(NPFullFile)]:
        pass
#    else:
#        sys.stderr.write("The file you have selected, '" + NPFullFile + "', is not a '.narrowPeak' file. This script will now terminate without calculating the FRiP score")
#        exit()
else:
    sys.stderr.write("THE FILE '" + str(NPFullFile) + "' CANNOT BE LOCATED. This script will now terminate without calculating the FRiP score")
    exit()
try:
    CountMMR = sys.argv[3]
    if CountMMR == "-M":
        CountMMR = True
    else:
        CountMMR = False
except IndexError:
    CountMMR = False
    pass
SAFFile = NPFullFile[0:(len(NPFullFile)-10)] + "saf"
# print(SAFFile)
KeepTempFiles = False
SigFigs = 4
salt = BAMFileName[0:len(BAMFileName)-4]
os.system("awk 'OFS=" + '"\\t"' + " {print $4, $1, $2+1, $3, " + '"+"' + "}' " + NPFullFile + " > " + SAFFile)
if CountMMR is True:
#   featureCounts_command = "featureCounts -M -f -p -a " + SAFFile + " -F SAF -o TempInfo_" + salt + ".txt " + BAMFullFile
    featureCounts_command = "featureCounts -M -f -p -a " + SAFFile + " -F SAF -o TempInfo_" + salt + ".txt " + BAMFullFile + " 2> /dev/null"
    sys.stderr.write("Running:\n" + featureCounts_command + "\n")
    os.system(featureCounts_command)
else:
#   featureCounts_command = "featureCounts -p -a " + SAFFile + " -F SAF -o TempInfo_" + salt + ".txt " + BAMFullFile
    featureCounts_command = "featureCounts -p -a " + SAFFile + " -F SAF -o TempInfo_" + salt + ".txt " + BAMFullFile + " 2> /dev/null"
    sys.stderr.write("Running:\n" + featureCounts_command + "\n")   
    os.system(featureCounts_command)
with open("TempInfo_" + salt + ".txt.summary") as outputfile:
    for i, line in enumerate(outputfile):
        if "Assigned" in line:
            assigned = float(line.split("	")[1])
        if "Unassigned_MultiMapping" in line:
            multimapped = float(line.split("	")[1])
        if "Unassigned_NoFeatures" in line:
            assignednofeatures = float(line.split("	")[1])
total = assigned + multimapped + assignednofeatures
FRiP = float(assigned)/float(total)
# os.system("clear")
print("FRiP score = " + str(round(FRiP*100, SigFigs)) + "%")

if KeepTempFiles is False:
    os.remove("TempInfo_" + salt + ".txt.summary")
    os.remove("TempInfo_" + salt + ".txt")
    os.remove(SAFFile)
