library(argparse, quietly=TRUE)

parser <- ArgumentParser(description="What cell types available for epigenetic heatmap")

source("epigeneticFeatures.R")

histoneorder <- epigenetic_features()
cat(histoneorder, sep='\n')

