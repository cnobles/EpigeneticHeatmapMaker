#### Test heatmap with following params ####
if(FALSE) {
    freeze <- 'hg18'
    insfile.name <- 'CD4_HIV-XMRV_293t_HERV-MLV-ASLV'
    windows <- c('10000')
    setName <- c('Beitzel-ASLV-293T-tva','Brady-HERV-293T-MseI-MultiPCR4-July08','Brady-HERV-293T-MseI-MultiPCR3-July08','Brady-HERV-293T-MseI-MultiPCR2-July08','Brady-HERV-293T-MseI-SinglePCR-July08','Brady-HERV-293T-ApoI-July08','Brady-HERV-293T-ApoI-June08','Brady-HERV-293T-MseI-MultiPCR1-July08','Brady-HERV-293T-MseI-MultiPCR4-June08','Brady-HERV-293T-MseI-MultiPCR3-June08','Brady-HERV-293T-MseI-MultiPCR2-June08','Brady-HERV-293T-MseI-MultiPCR1-June08','Brady-HERV-293T-MseI-SinglePCR-June08','Brady-HERV-293T-MseI-Dec2007','Wang-VA2-June-201-ExVivo-ApoI','Wang-VA2-June-201-ExVivo-Avr','Wang-VA2-June-202-ExVivo-Avr','Wang-VA2-June-202-ExVivo-ApoI','Wang-VA2-June-203-ExVivo-ApoI','Wang-VA2-June-203-ExVivo-Avr','Roth-XMRV-CD4-20100129Well1Tsp','Roth-XMRV-CD4-20100129Well1Mse','Roth-XMRV-CD4-20100129Well2Tsp','Roth-XMRV-CD4-20100129Well2Mse','Barr-HIV-VSV-Mac','Barr-HIVR5-Mac')
    Alias <- c('d)ASLV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','e)HERV-293T','b)HIV-CD4','b)HIV-CD4','b)HIV-CD4','b)HIV-CD4','b)HIV-CD4','b)HIV-CD4','c)XRMV-CD4','c)XRMV-CD4','c)XRMV-CD4','c)XRMV-CD4','a)HIV-MAC','a)HIV-MAC')
    meth.sets <- c('H2BK5me1','H3R2me1','H3R2me2','H3K4me1','H3K4me2','H3K4me3','H3K9me1','H3K9me2','H3K9me3','H3K27me1','H3K27me2','H3K27me3','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H4R3me2','H4K20me1','H4K20me3','PolII','H2AZ','CTCF','H3K4me1-HeLaS3_Untr','H3K4me1-HeLaS3_INFS','H3K4me3-HeLaS3_Untr','H3K4me3-HeLaS3_INFS','H3K27me3-CD133','H3K4me1-CD133','H3K9me1-CD133','H4K20me1-CD133','H3K27me1-CD133','H3K36me3-CD133','H3K4me3-CD133','H3K9me3-CD133','PolII-CD133','H2AZ-CD133','H3K9me2-HeLa','H3K9me3-HeLa','H3K27me3-HeLa','H3K36me3-HeLa','H4K20me3-HeLa','CTCF-HeLa','H4K20me1-HeLa','H2BK5me1-HeLa','H3K4me1-HeLa','H2AZ_hsalt-HeLa','H2AZ_lsalt-HeLa','H3K4me3-HeLa')
    act.sets <- c('H2AK5ac','H2AK9ac','H2BK120ac','H2BK12ac','H2BK20ac','H2BK5ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K36ac','H3K4ac','H3K9ac','H4K12ac','H4K16ac','H4K5ac','H4K8ac','H4K91ac','H3K27ac-HeLa') 
    tf.sets <- c('NRSF','STAT1-HeLaS3_Untr','STAT1-HeLaS3_INFS','STAT1-HeLaS3_INFS-2','ActivatedNucleosomes','RestingNucleosomes','Rest-CD4-HDAC1','Rest-CD4-HDAC2','Rest-CD4-HDAC3','Rest-CD4-HDAC6','Act-CD4-HDAC6','Rest-CD4-CBP','Rest-CD4-PCAF','Rest-CD4-p300','Rest-CD4-MOF','Rest-CD4-Tip60','Act-CD4-Tip60','H3_3-HeLa','p300-HeLa')
    inoutNuc <- 'yes'
}

#### Load required libraries & functions ####
libs <- c("survival", "pipeUtils", "colorspace", "hiAnnotator", 
          "plyr", "reshape2")
loaded <- sapply(libs, library, character.only=TRUE, quietly=TRUE)

#### set parameters for annotations & analysis ####
dbConn <- connectToIntSitesDB(user="3y996DKX7i", password="7j32ue9j3l")
annotPath <- file.path("/usr/local/Annotations/Epigenetic",freeze)
if(!file.exists(annotPath)) {
    stop("Required annotation directory doesn't exist: ", annotPath)
}

if(length(setName)!=length(Alias)) { 
    stop("Number of setName do not match number of Aliases.")
}
setsConAlias <- data.frame(setName, Alias, stringsAsFactors=FALSE)
print(setsConAlias)

#### load up the integration sites ####
sites <- getSitesFromDB(dbConn, setName=setsConAlias$setName, freeze, mrcs=TRUE)
sites <- merge(sites, setsConAlias, all.x = TRUE)
stopifnot(!any(is.na(sites$Alias)))

## Dereplicate sites within same alias
dereplication <- FALSE
if(dereplication) {
    sites.ins <- droplevels(subset(sites, type=='insertion'))
    sites.mrc <- droplevels(subset(sites, type=='match'))
    rows <- duplicated(with(sites.ins, paste0(Alias,Chr,Ort,Position)))
    sites.ins <- sites.ins[!rows,]
    sites <- rbind(sites.ins, 
                   subset(sites.mrc, Sequence %in% sites.ins$Sequence))
    sites <- droplevels(sites)
    rm("sites.ins","sites.mrc","rows")
}

if(nrow(sites)==0) { stop("No results found for given datasets.") }

## Limit the number of MRCs to the given limit if defined
limitMRCs <- FALSE
if(limitMRCs) {
    mrcPerSitesLimit <- 3
    sites.mrc <- droplevels(subset(sites, type=='match'))    
    rows <- ave(rep(TRUE, nrow(sites.mrc)), sites.mrc$Sequence, 
                FUN=function(x) cumsum(x) <= mrcPerSitesLimit)
    sites.mrc <- sites.mrc[rows,]
    sites <- rbind(droplevels(subset(sites, type=='insertion')), sites.mrc)
    rm("sites.mrc","rows")
}

sites <- makeGRanges(sites, soloStart=TRUE, chromCol="Chr", strandCol="Ort",
                     freeze=freeze)

species <- switch(freeze, hg18="Homo sapiens", mm8="Mus musculus")
sites <- keepStandardChromosomes(sites, style="UCSC", species=species)

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
todocombos <- expand.grid(setName=unique(sites$setName), histones=histoneorder, 
                          windows=windows, stringsAsFactors=FALSE)
todocombos$histone_window <- with(todocombos, 
                                  paste(histones, getWindowLabel(windows),sep="."))
counts <- count(as.data.frame(sites)[,c("setName","BID")], "setName")
todocombos <- merge(todocombos, counts)

sql <- sprintf("SELECT name, size, histone FROM psl_histones WHERE freeze='%s' AND histone IN ('%s') AND name IN ('%s')", freeze,
               paste(unique(todocombos$histone_window), collapse="','"), 
               paste(unique(todocombos$setName), collapse="','"))
res <- dbGetQuery(dbConn, sql)

## get counts of all sites in samples already done and compare with new sites ##
if(nrow(res)>0) {    
    todocombos$todo <- !with(todocombos, paste0(setName,freq,histone_window)) %in%
        with(res, paste0(name,size,histone))
} else {
    todocombos$todo <- TRUE
}

rows <- todocombos$todo
todoHistones <- unique(todocombos$histones[rows])
todoWindows <- unique(todocombos$windows[rows])
todoSets <- unique(todocombos$setName[rows])

#### Do combinations of setName-epigeneticFactor-window for newbies ####
sites <- sites[sites$setName %in% todoSets]

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
    
    sites <- as.data.frame(sites)
    
    ## Load the newbies to the DB ##
    ## load by histone & sets to conserve bandwidth and memory ##
    cols <- c(grep(paste(getWindowLabel(windows), collapse="|"), 
                   names(sites), value=TRUE), 
              grep("inout$", names(sites), value=TRUE))   
    for(f in cols) {
        message("Loading to DB ", f)
        for(x in unique(sites$setName)) {
            toload <- subset(sites, setName==x)[,c("Sequence","type",f)]
            size <- nrow(toload)
            
            ## undo addition of '_' within the histone name due to cleanColname ##
            if(!grepl("inout$",f)) {
                newcol <- histone_clean[[sub("(.+)\\.\\d+.+", "\\1", f)]]
                newcol <- paste0(newcol,sub(".+(\\.\\d+.+)", "\\1", f))
                if(newcol!=f) {
                    ## no need to add more column if colname didnt change! ##
                    toload[,newcol] <- toload[,f]
                    toload[,f] <- NULL
                }                
            } else {
                newcol <- f
            }
            toload <- rawToChar(serialize(toload,NULL,ascii = TRUE))
            sql <- sprintf("REPLACE INTO psl_histones VALUES ('%s',%d,'%s','%s','%s')",
                           x, size, freeze, newcol, toload)
            dbSendQuery(dbConn,sql)
        }
    }
    rm(sites)
}

#### get the data and make the ROC heatmap ####
message("Acquiring data for Heatmap")

if(tolower(inoutNuc)=="yes" & freeze=="hg18") {
    sql <- sprintf("SELECT name, histone, counts FROM psl_histones WHERE freeze='%s' AND histone IN ('%s') AND name IN ('%s')", freeze,
                   paste(c(unique(todocombos$histone_window),
                           "ActivatedNucleosomes_inout", 
                           "RestingNucleosomes_inout"), collapse="','"), 
                   paste(unique(todocombos$setName), collapse="','"))
} else {
    sql <- sprintf("SELECT name, histone, counts FROM psl_histones WHERE freeze='%s' AND histone IN ('%s') AND name IN ('%s')", freeze,
                   paste(unique(todocombos$histone_window), collapse="','"), 
                   paste(unique(todocombos$setName), collapse="','"))    
}
sites <- dbGetQuery(dbConn, sql)

## do unserialization ##
sites <- split(sites, sites$name)
sites <- lapply(names(sites), function(x) 
    cbind(setName=x,
          do.call(cbind, lapply(sites[[x]]$counts, 
                          function(y) unserialize(charToRaw(y))))))
sites <- rbind.fill(sites)

## make sure total sites obtained from DB is same as onces retreived earlier ##
counts <- arrange(counts, setName)
counts2 <- arrange(count(sites,"setName"),setName)
stopifnot(identical(counts$freq, counts2$freq))

sites <- merge(sites, setsConAlias, all.x = TRUE)
stopifnot(!any(is.na(sites$Alias)))
sites <- arrange(sites, Alias, Sequence, type)

message("Make Heatmap")
## order columns by histone order ##
cols <- grep(paste(getWindowLabel(windows), collapse="|"), names(sites), 
             value=TRUE)
bore <- sapply(strsplit(cols,"\\."),"[[",1)
rows <- match(histoneorder, bore)
cols <- cols[rows]
if(tolower(inoutNuc)=="yes" & freeze=="hg18") {
    cols <- c(cols, grep("inout$", names(sites), value=TRUE))
}

sites$good.row <- TRUE
rset <- with(sites, ROC.setup(good.row, type, Sequence, Alias))
roc.res <- ROC.strata(cols, rset, add.var=TRUE, sites)
#heatmap.dir <- file.path("../../Other/MethHeatMaps",insfile.name)
heatmap.dir = file.path("../Other/MethHeatMaps", insfile.name)
system(paste("mkdir", heatmap.dir, sep=" "))

dcol <- c(rep(rgb(185,185,0,max=185), 70), 
          rgb(185:0,185:0,0,max=185),
          rgb(0,0,0:220,max=220), 
          rep(rgb(0,0,220,max=220), 35))
ROCSVG(roc.res, heatmap.dir, colScale=dcol, overlayCols="white")

#### save the sites frame ####
## create a summary file of sites + counts to be displayed in the php script, 
## and write the command ran to /Library/WebServer/Documents/Insipid/RhmpFiles/heatMapCMDhistory.txt
filename <- paste(heatmap.dir,"/",insfile.name,".sites.frame.RData",sep="")
save(sites, file=filename, compress=TRUE)

summaryFrame <- ddply(sites, .(setName, Alias), summarize,
                      insertion=length(unique(Sequence)), 
                      match=length(Sequence))

cols <- c("setName","insertion","match","Alias")
summaryFrame <- summaryFrame[,cols]
colnames(summaryFrame) <- c("Setname","Insertions","Matches","Alias")
filename <- file.path(heatmap.dir,paste0(insfile.name,"-summary.txt"))
summaryFrame <- arrange(summaryFrame, Alias, Setname)
write.table(summaryFrame, filename, sep="\t", na="", row.names=FALSE, quote=FALSE)

cmd <- paste0("echo \"insfile.name='",insfile.name,
              "'; windows <- c('",paste(windows,sep="",collapse="','"),
              "'); setName <- c('",paste(setName,sep="",collapse="','"),
              "'); Alias <- c('",paste(Alias,sep="",collapse="','"),
              "'); meth.sets <- c('",paste(meth.sets,sep="",collapse="','"),
              "'); act.sets <- c('",paste(act.sets,sep="",collapse="','"),
              "'); tf.sets <- c('",paste(tf.sets,sep="",collapse="','"),
              "'); inoutNuc <- '",inoutNuc,"'; freeze <- '",freeze,
              "';\" | cat - methHeatMapMaker.R | R --no-save --no-restore")

write("********************************************************",append=TRUE,
      file="/Library/WebServer/Documents/Insipid/RhmpFiles/heatMapCMDhistory.txt")

write(cmd, append=TRUE,
      file="/Library/WebServer/Documents/Insipid/RhmpFiles/heatMapCMDhistory.txt")
