#! /usr/bin/env Rscript

message("\nFunction to convert a redundant MAList into a MAList with uniuque probes.")
message("USAGE: From within R, source(\"makeNonRedundantMAList.R\"); makeNonRedundantMAList(...)\n")

# This method can be used to summarize a redundant MAList
# by getting the median of each pair of M and A values 
# across multiple spots of the same probe (identified 
# as having identical $genes$ProbeName fields)
# 
# A summarized MA list is returned.
#
# This was created for Agilent yeast arrays 
# with 7 copies of each spot
#
# Similar code should work to summarize a RGList
# but this has not yet been implemented.
# 
# George Bell - 14 Sept 2009
#
# To use the function,
# 1 - Copy and paste this function into a R session
# 2 - Run a command like
#     my.MAList.nonRedundant = makeNonRedundantMAList(my.MAList.redundant)


makeNonRedundantMAList = function(my.MAList)
{
	if (!is(my.MAList, "MAList"))
	{
		stop("Your input must be a MAList.")
	}

	if (is.null(my.MAList$genes$ProbeName))
	{
		stop("\nYour input MAList must contain a $genes$ProbeName field, as this us used as the key for summarization.\n\n")
	}

	cat("Summarizing (by medians) your MAList by $genes$ProbeName key...\n") 

	num.samples = nrow(my.MAList$targets)

	ProbeNames.nr = unique(my.MAList$genes$ProbeName)
	# How many are there?
	num.unique.probes = length(ProbeNames.nr)

	# Make a new MAList data structure to hold the non-redundant MA data
	MAList.nr = new("MAList")

	# What fields does the original MAList contain??
	# names(my.MAList)	# [1] "targets" "genes"   "source"  "M"       "A"

	# Set $targets (array names) to those of original structure
	MAList.nr$targets = my.MAList$targets
	# Get number of columns of data
	num.cols = dim(my.MAList$targets)[1]

	# Set $source to that of original structure (ex: "agilent") 
	MAList.nr$source = my.MAList$source

	# Set $genes to non-redundant list of ProbeName, GeneName, SystematicName, ControlType
	# MA.norm.nr$genes = unique(MA.norm$genes[,c(5:7,4)])
	# Make sure we have these fields in the original data structure
	MAList.nr$genes = unique(my.MAList$genes[,c("ProbeName", "GeneName", "SystematicName", "ControlType")])

	# Sort structure by $ProbeName
	MAList.nr$genes = MAList.nr$genes[order(MAList.nr$genes$ProbeName),]

	# Make matrices to hold summarized M and A values
	A.medians = matrix(data=NA, ncol=ncol(my.MAList), nrow = num.unique.probes)
	colnames(A.medians) = colnames(my.MAList)
	rownames(A.medians) = MAList.nr$genes$ProbeName
	M.medians = matrix(data=NA, ncol=ncol(my.MAList), nrow = num.unique.probes)
	colnames(M.medians) = colnames(my.MAList)
	rownames(M.medians) = MAList.nr$genes$ProbeName
	
	# Iterate through probes one at a time
	for (i in 1:length(MAList.nr$genes$ProbeName))
	{
		my.probe = MAList.nr$genes$ProbeName[i]

		# Get all A or M values for this probe
		A.this.probe = my.MAList$A[my.MAList$genes$ProbeName == my.probe,]
		M.this.probe = my.MAList$M[my.MAList$genes$ProbeName == my.probe,]

		# If there's more than one probe for this gene
		if (length(A.this.probe) > num.samples)
		{
			# Get A value as median by probe
			A.medians[i,] = apply(A.this.probe, 2, median)

			# Get M value as median by gene
			M.medians[i,] = apply(M.this.probe, 2, median)
		} else {		# Just one probe for this gene; nothing to do 
			A.medians[i,] = A.this.probe
			M.medians[i,] = M.this.probe
		}

		# Put medians into data structure
		rownames(A.medians)[i] = my.probe
		rownames(M.medians)[i] = my.probe
	}

	# Change name to fill final "MAList" data structure
	MAList.nr$A = A.medians
	MAList.nr$M = M.medians
	
	cat("Completed summaries of", length(MAList.nr$genes$ProbeName), "probes.\n") 

	MAList.nr
}
