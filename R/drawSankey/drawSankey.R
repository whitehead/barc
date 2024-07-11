#! /usr/bin/env Rscript

# Sankey diagrams can illustrate flux in dynamic processes.  In this application
# a Sankey diagram is generated to depict the relationships between two sets
# of class labels.

# A good way to make Sankey diagrams with R is to use the riverplot
# package, which requires nodes and edges as input.  The function below
# generates nodes and edges from an input matrix.  This input matrix is
# expected to be a 3-column data.frame.  The makeRiver function is used
# below to convert nodes and edges into a riverplot object before plotting.

# This is the expected format for matrix input:

# Column 1: object identifier.
# Column 2: class labels 1.
# Column 3: class labels 2.

# Note that missing class labels will be converted to NA and will not
# contribute to edge widths.  To capture the birth and death of objects,
# define and label the appropriate null states in the input file (i.e. fill in
# the missing values).

# Finally, the colours for the sources and sinks in the Sankey diagram with
# two sets of class labels are given by the input arguments col1 and col2,
# which can be valid R colour names or RGB codes.

# Version 1.0, 24 Sept., 2020
# Version 1.1, 11 July, 2024, accommodate v.0.6 -> v.0.10 changes in riverplot. 
# Troy Whitfield - Bioinformatics and Research Computing, Whitehead Institute

argv <- commandArgs()
if('--args' %in% argv){
   nignore <- which (argv == "--args")
   argv <- argv[-(1:nignore)]
}
if((length (argv) != 4) | (argv[1] == "/usr/lib/R/bin/exec/R")) {
   message("Draw a 2-level Sankey diagram.")
   message("USAGE: ./drawSankey.R input colour1 colour2 output.pdf")
   message("Example: ./drawSankey.R riverData.dat orange green sankeyPlot.pdf\n")
   message("Note: The input should be a 3 column matrix with unique row names/IDs in column 1")
   message("and two sets of class labels in columns 2 and 3.")
   message("Colours can be valid R colour names or RGB codes and output should be a user-specified pdf file name.\n")

   # stop("Follow the expected input format above.")
   q()
}

message("Drawing a Sankey diagram for N objects with two sets of class labels.")
#message("For example, N cells labelled by treatment (class 1) and response (class 2).")
message("Colouring for the diagram has been set by the user.\n")

suppressMessages(require("riverplot"))
suppressMessages(require("Cairo")) # Prettier fonts for plotting.

# Assign variables from input.
infile<-argv[1]
col1<-argv[2]
col2<-argv[3]
outfile<-argv[4]

# Define useful functions.
riverPrepTwoLevel<-function(dat,col1,col2){
	levs<-length(dat[1,])
	initStates<-na.omit(unique(dat[,2])) # Assumes that column 1 is obj. ID.
	finalStates<-na.omit(unique(dat[,levs])) # NA class assignments excluded
	ninitStates<-length(initStates)
	nfinalStates<-length(finalStates)
# Create nodes for initial and final states.	
	initNodes<-c(paste0("i",1),1,initStates[1])
	if (ninitStates > 1){
	   for (i in 2:ninitStates){
	       initNodes<-rbind(initNodes,c(paste0("i",i),1,initStates[i]))
	   }
        }
	finalNodes<-c(paste0("f",1),2,finalStates[1])
	if (nfinalStates > 1){
	   for (i in 2:nfinalStates){
	       finalNodes<-rbind(finalNodes,c(paste0("f",i),2,finalStates[i]))
	   }
        }
	nodes<-rbind(initNodes,finalNodes)
# Compute weights for edges.
	edges<-array(NA,c(length(initNodes[,3])*length(finalNodes[,3]),3))
	edgeCount<-0
	for (i in 1:length(initNodes[,3])){
    	    for (j in 1:length(finalNodes[,3])){
    	    	edgeCount<-edgeCount+1
		flux<-length(which((dat[,2] == initStates[i]) & dat[,3] == finalStates[j]))
		edges[edgeCount,1]<-initNodes[i,1]
		edges[edgeCount,2]<-finalNodes[j,1]
    		edges[edgeCount,3]<-flux
    	    }
	 }

	 colnames(nodes)<-c("ID","x","labels")
	 nodes<-as.data.frame(nodes)
	 nodes$ID<-as.character(nodes$ID)
	 nodes$x<-as.numeric(nodes$x)
	 nodes$labels<-as.character(nodes$labels)

	 colnames(edges)<-c("N1","N2","Value")
	 edges<-as.data.frame(edges)
	 edges$N1<-as.character(edges$N1)
	 edges$N2<-as.character(edges$N2)
	 edges$Value<-as.numeric(as.character(edges$Value))

# Use two scalar input colours to assign styles.
         nstyles<-list()
	 for (i in 1:ninitStates){
	     nstyles[[i]]<-list(col=col1,lty=0,textpos=2)
	 }
	 for (j in (i+1):(ninitStates+nfinalStates)){
	     nstyles[[j]]<-list(col=col2,lty=0,textpos=4)
	 }
	 attr(nstyles,"names")<-nodes$ID # Convert to named list.
	 madeRiver<-makeRiver(nodes,edges,node_styles=nstyles)
	 return(madeRiver)
}

make.true.NA<-function(x){
	if(is.character(x)||is.factor(x)){
		is.na(x) <- x=="NA"; x
	} else {x}
}

# Read in data file.
data<-read.table(infile,header=TRUE,fill=TRUE,na.strings="",stringsAsFactors=FALSE,sep="\t")

# Make sure that all "NA"s are true NA.
data[] <- lapply(data, make.true.NA)

# Create a riverplot object from the input data.
myRiver<-riverPrepTwoLevel(data,col1,col2)

sink("/dev/null") # Suppress default reporting from riverplot.
cairo_pdf(outfile,width = 8.5, height = 8.0, family = "Helvetica", bg="white", pointsize=14)
riverplot(myRiver,gravity="center",srt=0,nodewidth=0,plot_area=0.75,fix.pdf=FALSE)
dev.off()
sink()
