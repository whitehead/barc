#! /usr/bin/env Rscript

# This script allows easy command line access to computing b-scores
# (see C. Brideau et al. "Improved Statistical Methods for Hit Selection
# in High-Throughput Screening", Journal of Biomolecular Screening 8,
# 634-647 (2003)) from raw single plate matrices for high throughput
# screens.

# Version 1.0, 15 June, 2021
# Troy Whitfield - Bioinformatics and Research Computing, Whitehead Institute

argv <- commandArgs()
if('--args' %in% argv){
   nignore <- which (argv == "--args")
   argv <- argv[-(1:nignore)]
}
if((length(argv) != 3) | (argv[1] == "/usr/lib/R/bin/exec/R")) {
   message("Convert a raw input matrix (e.g. from high throughput")
   message("screening experiments) into a matrix of B-scores, thereby")
   message("correcting for plate effects.\n")
   message("USAGE: ./getBscores.R input scale output\n")
   message("Example: ./getBscores.R input.csv additive output.csv\n")
   message("The input and output files are each a matrix in comma")
   message("separated variable format.\n")
   message("The 'scale' keyword can be either 'additive' (default)")
   message("or 'multiplicative'.  Multiplicative data (e.g. ratios)")
   message("will be (base 2) logaritmically transformed before b-scores")
   message("are computed.")
   q()
}

message("Computing B-scores.\n")

# Assign variables from input.
infile<-argv[1]
scale<-argv[2]
outfile<-argv[3]

optScale<-c("additive","multiplicative","default")
if (scale == "default"){
   scale<-"additive"
} else if (!(scale %in% optScale)){
   scale<-"additive"
}

# Read input matrix.
raw<-read.table("input.csv",header=FALSE,sep=",")

if (scale == "multiplicative"){
   if (any(raw) <= 0){
      q("Exiting.  Input matrices on the multiplicative scale must have all positive entries.")
   } else {
     raw<-log(raw,base=2)
   }
}

# Estimate residuals using Tukey median polish.
sink("/dev/null") # Suppress default reporting.
MP<-medpolish(raw,eps=0.001,maxiter=10000,trace.iter=TRUE,na.rm=FALSE)
sink()

# Scale the residuals by MAD.
bscores<-MP$residuals/mad(MP$residuals)

# Write the results to file.
write.table(bscores,outfile,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")
message("Done. B-scores written to file.")
