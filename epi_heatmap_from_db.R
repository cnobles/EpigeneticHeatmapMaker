library(argparse, quietly=TRUE)

parser <- ArgumentParser(description="Make epigenetic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
parser$add_argument("-o", "--output_dir", type="character", default="epi_heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg18", 
    help="reference genome used for all samples(only hg18 and mm8 are supported at present)")
parser$add_argument("-s", "--sites_group", type="character", default="intsites_miseq.read", 
    help="which group to use for connection")
parser$add_argument("-t", "--cell_types", type="character", 
    help="file with cell types to use: each cell type on separate line(epi_cell_types.R prints all avalable types)")

args <- parser$parse_args()

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

library('yaml')

# Load configuration file
if (!file.exists(args$c)) stop("the configuration file can not be found.")
config <<- yaml.load_file(args$c)

### source("epigeneticFeatures.R")
### source("epigeneticHeatMapMaker.R")

source(file.path(codeDir, "epigeneticFeatures.R"))
source(file.path(codeDir, "epigeneticHeatMapMaker.R"))

library(DBI, quietly=TRUE)
library(RMySQL, quietly=TRUE)
library(dplyr)

library(BSgenome)
library(intSiteRetriever)
library(devtools)
source_url('https://raw.githubusercontent.com/BushmanLab/genomicHeatmapMaker/master/utils.R')

available_epigenetic_ref_genomes <- c ("hg18", "mm8")

referenceGenome <- args$ref_genome
if ( ! referenceGenome %in% available_epigenetic_ref_genomes) {
    message("Epigenetic annotation only available for:")
    message(paste(available_epigenetic_ref_genomes, collapse=" "))
    stop(0)
}
heat_map_result_dir <- args$output_dir 
### annotation <- args$annotation_path
annotation <- config$epigeneticDataDirectory
annotation <- file.path(annotation, referenceGenome)

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, "not found")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

# Connect to my database
if (config$dataBase == 'mysql'){
   stopifnot(file.exists("~/.my.cnf"))
   stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))
   dbConn <- dbConnect(MySQL(), group=config$mysqlConnectionGroup)
   info <- dbGetInfo(dbConn)
   connection <- src_sql("mysql", dbConn, info = info)
}else if (config$dataBase == 'sqlite') {
   dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteIntSitesDB)
   info <- dbGetInfo(dbConn)
   connection <- src_sql("sqlite", dbConn, info = info)
   dbConn2 <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteSampleManagement)
   info2 <- dbGetInfo(dbConn2)
   connection2 <- src_sql("sqlite", dbConn2, info = info2)
} else { stop('Can not establish a connection to the database') }


histoneorder_for_heatmap <- epigenetic_features()
if ( ! is.null(args$cell_types)) {
    # use only these features
    historder_subset <- read.csv(args$cell_types, header=FALSE, stringsAsFactors=FALSE)$V1
    stopifnot(all(historder_subset %in% histoneorder_for_heatmap)) 
    histoneorder_for_heatmap <- historder_subset
}
message("epigenetic heatmap for the following cell types:")
cat(histoneorder_for_heatmap, sep='\n')

make_epi_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, 
    connection, annotation, histoneorder_for_heatmap)

