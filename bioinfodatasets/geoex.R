# R example for working with GEO data using GEOquery and limma for
# Seattle SIGKDD data club.
# (riffing on auto generated R script for GEO2R service on GEO site)

## libraries for everything
library(Biobase)
library(GEOquery)
library(limma)

###
# Series info:
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6275
#
# Note: we are getting the already normalized data. Different procedure
# to start with raw data
###

## get data
# 1. series data
gset <- getGEO("GSE6275", GSEMatrix =TRUE)
gset <- gset[[1]]

# 2. NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)

# reset column names
fvarLabels(gset) <- make.names(fvarLabels(gset))

##
# Data set has 36 samples corresponding to different brain regions of three
# strains of mice. Samples from each strain are pooled into three groups.

# group names for all samples
sml <- c("G0","G0","G1","G2","G2","G0","G0","G1","G2","G2","G0","G0","G1","G2","G2","G0","G0","G1","G2","G2","G0","G1","G1","G2","G0","G1","G1","G2","G0","G1","G1","G2","G0","G1","G1","G2");

# strain names corresponding to group names
labels <- c("stargazer","tottering","lethargic")

###
# plot the samples as is

oex <- exprs(gset)[,order(sml)]
osml <- sml[order(sml)]
ofl <- as.factor(sml)

palette(c("#dfeaf4","#f4dfdf","#f2cb98", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE6275", '/', annotation(gset), " selected samples", sep ='')
boxplot(oex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=ofl)
legend("topleft", labels, fill=palette(), bty="n")

###
# limma setup

# 1. log2 transform
ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex)

# 2. setup model for individual expression
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)

# 3. model for comparisons across three groups
cont.matrix <- makeContrasts(G2-G0, G1-G0, G2-G1, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# 4. grab table of best 250 results
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

###
# use NCBI annotation (which genes)
# note: uses platform data from beginning

ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))

##dump the file 
# but change the directory to yours first!
#setwd("~/Documents/projects/dataclub/demo")
#write.csv(tT,"toptable.csv",row.names=FALSE)

## inspect normalized values
#barplot(oex["1427328_a_at",])