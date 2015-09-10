make_epi_heatmap <- function(sampleName_GTSP, referenceGenome, output_dir, connection, annotation) {
    freeze <- referenceGenome
    insfile.name <- output_dir

    windows <- c('10000')
    if ( ! "label" %in% colnames(sampleName_GTSP)) {
        sampleName_GTSP$label <- sampleName_GTSP$GTSP
    }
    sampleName_GTSP <- select(sampleName_GTSP, sampleName, GTSP, label)
    sampleName_GTSP$refGenome <- rep(referenceGenome, nrow(sampleName_GTSP))
    setName <- sampleName_GTSP$sampleName
    Alias <- sampleName_GTSP$label

    # datasets
    meth.sets <- c('H2BK5me1','H3R2me1','H3R2me2','H3K4me1','H3K4me2','H3K4me3','H3K9me1','H3K9me2','H3K9me3','H3K27me1','H3K27me2','H3K27me3','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H4R3me2','H4K20me1','H4K20me3','PolII','H2AZ','CTCF','H3K4me1-HeLaS3_Untr','H3K4me1-HeLaS3_INFS','H3K4me3-HeLaS3_Untr','H3K4me3-HeLaS3_INFS','H3K27me3-CD133','H3K4me1-CD133','H3K9me1-CD133','H4K20me1-CD133','H3K27me1-CD133','H3K36me3-CD133','H3K4me3-CD133','H3K9me3-CD133','PolII-CD133','H2AZ-CD133','H3K9me2-HeLa','H3K9me3-HeLa','H3K27me3-HeLa','H3K36me3-HeLa','H4K20me3-HeLa','CTCF-HeLa','H4K20me1-HeLa','H2BK5me1-HeLa','H3K4me1-HeLa','H2AZ_hsalt-HeLa','H2AZ_lsalt-HeLa','H3K4me3-HeLa')
    act.sets <- c('H2AK5ac','H2AK9ac','H2BK120ac','H2BK12ac','H2BK20ac','H2BK5ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K36ac','H3K4ac','H3K9ac','H4K12ac','H4K16ac','H4K5ac','H4K8ac','H4K91ac','H3K27ac-HeLa') 
    tf.sets <- c('NRSF','STAT1-HeLaS3_Untr','STAT1-HeLaS3_INFS','STAT1-HeLaS3_INFS-2','ActivatedNucleosomes','RestingNucleosomes','Rest-CD4-HDAC1','Rest-CD4-HDAC2','Rest-CD4-HDAC3','Rest-CD4-HDAC6','Act-CD4-HDAC6','Rest-CD4-CBP','Rest-CD4-PCAF','Rest-CD4-p300','Rest-CD4-MOF','Rest-CD4-Tip60','Act-CD4-Tip60','H3_3-HeLa','p300-HeLa')
    inoutNuc <- 'yes'

    #### Load required libraries & functions ####
    libs <- c("survival", "pipeUtils", "colorspace", "hiAnnotator", 
              "plyr", "reshape2")
    loaded <- sapply(libs, library, character.only=TRUE, quietly=TRUE)

    #### set parameters for annotations & analysis ####
    dbConn <- connection
    annotPath <- annotation
    if( ! file.exists(annotPath)) {
        stop("Required annotation directory doesn't exist: ", annotPath)
    }

    add_label <- function(sites, sampleName_GTSP) {
        sites_GTSP <- merge(sites, sampleName_GTSP)
        sites_GTSP$sampleName <- sites_GTSP$label
        sites_GTSP$refGenome <- NULL # not needed downstream
        sites_GTSP$GTSP <- NULL # not needed downstream
        sites_GTSP$label <- NULL
        sites_GTSP
    }

    sites <- getUniqueSites(sampleName_GTSP, connection)
    sites$type <- "insertion"
    sites <- add_label(sites, sampleName_GTSP)
    
    mrcs <- getMRCs(sampleName_GTSP, connection)
    mrcs$type <- "match"
    mrcs <- add_label(mrcs, sampleName_GTSP)

    sites_mrcs <- rbind(sites, mrcs)

    sites <- makeGRanges(sites_mrcs, soloStart=TRUE, 
        chromCol="chr", strandCol="strand", startCol='position', freeze=freeze)

    #species <- switch(freeze, hg18="Homo sapiens", mm8="Mus musculus")
    #sites <- keepStandardChromosomes(sites, style="UCSC", species=species)

    #### Set the order of histone modifications to be displayed on the heatmap ####
    display.order <- c("H3K9me2", "H3K9me2-HeLa", "H3K9me3", "H3K9me3-MEF", 
                       "H3K9me3-CD133", "H3K9me3-HeLa", "H3K27me2", "H3K27me3", 
                       "H3K27me3-MEF", "H3K27me3-CD133", "H3K27me3-HeLa", "H3K14ac",
                       "H2AK9ac", "H3K23ac", "H3K36me1", "H3K36me3", "H3K36me3-MEF",
                       "H3K36me3-CD133", "H3K36me3-HeLa", "H3K27me1", "H3K27me1-CD133",
                       "H3R2me1", "H3R2me2", "H4R3me2", "H4K20me3", "H4K20me3-HeLa",
                       "H2AK5ac", "H4K16ac", "H4K12ac", "CTCF", "CTCF-HeLa", 
                       "CTCF-Jurkat", "PolII", "PolII-CD133", "H3K79me1", "H3K79me2",
                       "H3K79me3", "H4K20me1", "H4K20me1-CD133", "H4K20me1-HeLa", 
                       "H2BK5me1", "H2BK5me1-HeLa", "H3K4me1", "H3K4me1-CD133", 
                       "H3K4me1-HeLa", "H3K4me1-HeLaS3_Untr", "H3K4me1-HeLaS3_INFS",
                       "H3K4me2", "H3K9me1", "H3K9me1-CD133", "H4K8ac", "H4K5ac", 
                       "H2AZ", "H2AZ-CD133", "H2AZ_hsalt-HeLa", "H2AZ_lsalt-HeLa", 
                       "H3K4me3", "H3K4me3-MEF", "H3K4me3-CD133", "H3K4me3-HeLa",
                       "H3K4me3-HeLaS3_Untr", "H3K4me3-HeLaS3_INFS", "H2BK12ac", 
                       "H3K9ac", "H3K18ac", "H3K27ac", "H3K27ac-HeLa", "H2BK5ac", 
                       "H3K36ac", "H3K4ac", "H2BK20ac", "H4K91ac", "H2BK120ac", 
                       "NRSF", "STAT1-HeLaS3_Untr", "STAT1-HeLaS3_INFS", 
                       "STAT1-HeLaS3_INFS-2", "ActivatedNucleosomes", 
                       "ActivatedNucleosomes.inout", "RestingNucleosomes", 
                       "RestingNucleosomes.inout", "Rest-CD4-HDAC1", 
                       "Rest-CD4-HDAC2", "Rest-CD4-HDAC3", "Rest-CD4-HDAC6", 
                       "Act-CD4-HDAC6", "Rest-CD4-CBP", "Rest-CD4-PCAF", 
                       "Rest-CD4-p300", "Rest-CD4-MOF", "Rest-CD4-Tip60", 
                       "Act-CD4-Tip60", "H3_3-HeLa", "p300-HeLa", 
                       "Brd4_2_1H_minus_100w_20s-Uwe", "Brd4_6_0H_minus_100w_20s-Uwe")
    histoneorder <- intersect(display.order, c(meth.sets, act.sets, tf.sets))
    windows <- as.integer(windows)
    windows_names <- structure(windows, names=getWindowLabel(windows))

    ## for mapping hiAnnotator::cleanColname colnames to original names later ##
    histone_clean <- structure(histoneorder, names=cleanColname(histoneorder))

    #### See if combinations of setName-epigeneticFactor-window exist ####
    ## get counts of all sites in samples ##
    todocombos <- expand.grid(setName=unique(sites$sampleName), histones=histoneorder, 
                              windows=windows, stringsAsFactors=FALSE)
    todocombos$histone_window <- with(todocombos, 
                                      paste(histones, getWindowLabel(windows),sep="."))

    todocombos$todo <- TRUE

    rows <- todocombos$todo
    todoHistones <- unique(todocombos$histones[rows])
    todoWindows <- unique(todocombos$windows[rows])

    if(length(sites)>0) {
        for(f in todoHistones) {
            message(f)
            load(file.path(annotPath, paste0(f, ".RData")))
            
            if(length(sites)>1e6) {
                sites <- getFeatureCounts(sites, epigenData, f, width=todoWindows,
                                          doInChunks=TRUE)
            } else {
                sites <- getFeatureCounts(sites, epigenData, f, width=todoWindows)
            }
            
            if(tolower(inoutNuc)=="yes" & freeze=="hg18" & 
                   f %in% c("ActivatedNucleosomes","RestingNucleosomes")) {        
                methinout <- paste(f,"inout",sep=".")
                epigenData$name <- "bore"
                sites <- getSitesInFeature(sites, epigenData, methinout, asBool=TRUE)
            }
            rm(epigenData)
        }
    }

    sites <- as.data.frame(sites)

    # not DRY: the same code is used in genomic heatmap
    get_annotation_columns <- function(sites) {
      granges_column_names <- c("seqnames", "start", "end", "width", "strand")
      int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
      required_columns <- unique(c(
        granges_column_names, int_site_column_names, "type"))
      stopifnot(all(required_columns %in% names(sites)))
      setdiff(names(sites), required_columns)
    }

    message("Make Heatmap")
    
    annotation_columns <- get_annotation_columns(sites)
    
    rset <- with(sites, ROC.setup(
      rep(TRUE, nrow(sites)), type, siteID, sampleName))
    roc.res <- ROC.strata(annotation_columns, rset, add.var=TRUE, sites)

    dcol <- c(rep(rgb(185,185,0,max=185), 70), 
              rgb(185:0,185:0,0,max=185),
              rgb(0,0,0:220,max=220), 
              rep(rgb(0,0,220,max=220), 35))
    ROCSVG(roc.res, output_dir, colScale=dcol, overlayCols="white")
}
