source("epigeneticHeatMapMaker.R")

library(argparse, quietly=TRUE)
library(DBI, quietly=TRUE)
library(RMySQL, quietly=TRUE)
library(dplyr)

library(BSgenome)
library(intSiteRetriever)

parser <- ArgumentParser(description="Make epigenetic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-o", "--output_dir", type="character", default="heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg18", 
    help="reference genome used for all samples(only hg18 supported at present)")
parser$add_argument("-s", "--sites_group", type="character", default="intsites_miseq", 
    help="which group to use for connection")
parser$add_argument("-a", "--annotation_path", type="character", default="/media/THING1/dryga/Epigenetic/hg18", 
    help="epigenetic annotation path with RData files, e.g. H2BK120ac, H3K4ac, NRSF")

args <- parser$parse_args()

referenceGenome <- args$ref_genome
heat_map_result_dir <- args$output_dir 
annotation <- args$annotation_path

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, "not found")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

connection <- dbConnect(MySQL(), group=args$sites_group)
info <- dbGetInfo(connection)
connection <- src_sql("mysql", connection, info = info)

make_epi_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, connection, annotation)

