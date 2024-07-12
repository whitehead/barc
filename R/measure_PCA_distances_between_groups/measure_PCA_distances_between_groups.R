#! /usr/bin/env Rscript

###
###  Compare distances between samples in PC space
###  Input is an output file from PCA (and a design file indicating which samples are replicates of which groups)
###
###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  Sample command: ./measure_PCA_distances_between_groups.R My_experiment.PCA_distances.txt My_experiment_design.txt My_experiment
###  
###  Version 1.0: 17 June 2021
###


offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./measure_PCA_distances_between_groups.R" }

PCA.file = commandArgs()[6 - offset]
design.file = commandArgs()[7 - offset]
experiment.name = commandArgs()[8 - offset]

if (is.na(experiment.name))
{
	message("\nGiven the coordinates of a group of samples in PC space (from a PCA),")
	message("  calculate distances between all pairs of points, and then determine")
	message("  distance between each sample type (summarizing across replicates).\n")
	message("USAGE: ", this.script, " PCA.coordinates.file design.file experiment.name")
	message("Ex:    ", this.script, " My_experiment.PCA_distances.txt My_experiment_design.txt My_experiment\n")
	message("Notes:")
	message("  The rows of the PCA coordinates and the design files should match.\n")
	quit()
}

first.PCA.to.use = 1
second.PCA.to.use = 2

# Try measuring distances
pca.values = read.delim(PCA.file)
design = read.delim(design.file, row.names=1)
pca.values.PC1.PC2 = pca.values[,c(first.PCA.to.use, second.PCA.to.use)]

# One key command: calculate distance between each sample and each other sample
PC1.PC2.distances = as.matrix(dist(pca.values.PC1.PC2, method = "euclidean", upper=F))

# Get a vector of all of our groups
group.IDs = design[,1]
groups = as.vector(unique(group.IDs))

# Make a matrix to hold our distances
biggest.group.size = max(table(group.IDs))
all.distances.between.groups = matrix(data=NA, ncol=sum(1:length(groups)), nrow=biggest.group.size*biggest.group.size)
colnames(all.distances.between.groups)=1:ncol(all.distances.between.groups)
comparison.num = 1

for (i in 1:length(groups))
{
	for (j in 1:length(groups))
	{
		if (i >= j)
		{
			group.1 = groups[i]
			group.2 = groups[j]
			this.comparison = paste(group.1, group.2, sep=".vs.")
			# print(this.comparison)
			colnames(all.distances.between.groups)[comparison.num] = this.comparison

			# Get subset of these distances
			distances.this.comparison = as.vector(PC1.PC2.distances[group.IDs == group.1 , group.IDs == group.2])
			all.distances.between.groups[1:length(distances.this.comparison), comparison.num] = distances.this.comparison

			comparison.num = comparison.num + 1
		}
	}
}

medians.by.comparison = apply(all.distances.between.groups, 2, median, na.rm=T)
comparison.by.increasing.median = names(sort(medians.by.comparison))

pdf(paste(experiment.name, "comparisons.distances.boxplot.pdf", sep="."), w=11, h=8.5, useDingbats=FALSE)
par(las=2, mai=c(2,1,1,1))
boxplot(all.distances.between.groups[, comparison.by.increasing.median], col="wheat", main=paste(experiment.name, ": distances in PCA space between groups", sep=""),
	ylab="PC distance")
foo = dev.off()

write.table(all.distances.between.groups, file=paste(experiment.name, "group_comparisons.distances.txt", sep="."), sep="\t", quote=F, row.names=F)

message("\nAll done! -- see table and boxplot output files.  Table values can now be compared using ANOVA.\n")
