#! /usr/bin/Rscript

require("ggplot2")

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
if (length(args)==1){
  sName = "./"
} else {
  sName = args[2]
}

myfile <- paste(sName,".tsv", sep="")

mydata <- read.csv(myfile,header=TRUE, sep="\t")

pdf(gsub(" ","",paste(sName,".graph.pdf")), width = 25, height = 25)
myplot <- ggplot(mydata, aes(x = SampleId, y = Reads))+geom_col(aes(fill = Legend), width = 0.75)
myplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()