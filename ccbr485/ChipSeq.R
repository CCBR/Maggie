source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("csaw")
biocLite("GenomicAlignments")
biocLite("Rsamtools")
library(Rsubread)
library(limma)
library(edgeR)
library(Rsamtools)

setwd("/Volumes/data/ChipSeq/fastq")
load("/Volumes/data/ChipSeq/fastq/ChipSeq_092613.RData")

# Run Alignments and generate bam files
ref="/data/maggiec/RNASeq/Genomes/mm9/chr_all3.fa"
buildindex(basename="reference_index",reference=ref)
read1 = read.table("samples2.txt")
align(index="reference_index",readfile1=read1)

# Bring in bam files and run featurecounts: 
targets <- readTargets()

#> targets
#Bam CellType Group
#1      Jmjd3WT_H3K27.sort.bam       WT     1
#2 Sample_DKO_4SP_K27.sort.bam       KO     1
#3    Sample_SM_DKO_2.sort.bam       KO     2
#4     Sample_SM_WT_2.sort.bam       WT     2
#5   Sample_OTII_Cont.sort.bam       WT     3
#6    Sample_OTII_DKO.sort.bam       KO     3

#saf <- "mouse_prom_3000.txt"
#saf <- "mouse_prom_2000.txt"
saf <- "mouse_prom_2000_rev.txt"
fc <- featureCounts(files=targets$Bam,nthreads=32,annot.ext=saf,reportReads=FALSE)
fc$counts[rownames(fc$counts)=="S1pr1",]

##Create DGEList from counts:
#Group <- factor(targets$Group, levels = c(1,2,3))
Group <- factor(targets$Group, levels = c(1,2))
CellType <- factor(targets$CellType, levels = c("WT","KO"))

y <- DGEList(counts=fc$counts[,1:6], group=CellType, genes=fc$annotation[,1:6])
y <- DGEList(counts=fc$counts[,1:4], group=CellType, genes=fc$annotation[,1:4])
y <- calcNormFactors(y)
#isexpr <- rowSums(cpm(y) > 1) >= 2
#y <- y[isexpr,]
#y_rpkm <- rpkm(y,y$genes)

design <- model.matrix(~Group+CellType, data=y$samples)
y <- estimateDisp(y, trend="none", robust=TRUE)
summary(y$prior.df)
sqrt(y$common.disp)

pdfPath=paste('bcv_final.pdf',sep = "")
pdf(pdfPath)
plotBCV(y)
dev.off()

###Calculate offsets for each pair
#OS1=calcNormOffsetsforChIP(y$counts[,1],y$counts[,2])
#OS2=calcNormOffsetsforChIP(y$counts[,4],y$counts[,3])
#OS3=calcNormOffsetsforChIP(y$counts[,5],y$counts[,6])
# OS1=normalizeChIPtoInput(y$counts[,1],y$counts[,2])
# OS2=normalizeChIPtoInput(y$counts[,4],y$counts[,3])
# OS3=normalizeChIPtoInput(y$counts[,5],y$counts[,6])
# 
# a = rep(0,22668)
# OS=cbind(a,OS1,OS2,a,a,OS3)
# OS=cbind(a, log(OS1$scaling.factor), log(OS2$scaling.factor),a,a,
#          log(OS3$scaling.factor))
# y$offset=OS

###Do GLM statistical testing #######
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

fit <- glmFit(y, design, prior.count=2)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=4)
#lrt <- glmLRT(fit,coef=3)
topTags(lrt)
finalres = topTags(lrt,sort="none",n=Inf)

o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags, lowess=TRUE)
abline(h=c(-1, 1), col="blue")

y_logcpm=cpm(y, log=TRUE, prior.count=2)
A=rowMeans(y_logcpm)
M=rowMeans(y_logcpm[,4:6])-rowMeans(y_logcpm[,1:3])

y_norm=normalize.loess(y_logcpm, subset = sample(1:(dim(y_logcpm)[1])))
nM=rowMeans(y_norm[,4:6])-rowMeans(y_norm[,1:3])
nA=rowMeans(y_norm)
plot(nA,nM)
reg1 <- lm(nM~nA)
abline(reg1, col = "red")

#Read in gene names:
#genetable=read.table("mouse_prom_3000.txt", sep = "\t")
genetable=read.table("mouse_prom_2000.txt", sep = "\t")
colnames(genetable)=c("TranscriptID","Chr","Start","End","Strand","Gene")

finaltab=(cbind(genetable$Gene,finalres$table,lrt$fitted.values))
write.table(finaltab,file="H3K27_ChipSeq_2000.txt",row.names=FALSE,sep="\t",quote=FALSE)

##################################################################
##If without replicates, run the following:
##Run Normalization and calculate pvals based on individual response-input pairs:
y3=y[,1:2]
isexpr3 = rowSums(cpm(y3) > 0) >=1
x3 = y3[isexpr3,]
x3 <- calcNormFactors(x3)

Norm3 = normalizeChIPtoInput(x3$counts[,1],x3$counts[,2])
names(Norm3$pmid.value) <- x3$genes[,1]

y1=y[,3:4]
isexpr1 = rowSums(cpm(y1) > 0) >=1
x1 = y1[isexpr1,]
Norm = normalizeChIPtoInput(x1$counts[,2],x1$counts[,1])
names(Norm$pmid.value) <- x$genes[,1]

y2=y[,5:6]
isexpr2 = rowSums(cpm(y2) > 0) >=1
x2 = y2[isexpr2,]
Norm2 = normalizeChIPtoInput(x2$counts[,2],x2$counts[,1])
names(Norm2$pmid.value) <- x2$genes[,1]

##Calculate the offsets:
Normval = calcNormOffsetsforChIP(x1$counts[,2],x1$counts[,1])

##Calculate fold changes:
predlfc<-predFC(x1,design,dispersion=dispersion,offset=Normval,prior.count=0.5)


# generate counts for a two group experiment with n=2 in each group and 100 genes
dispersion <- 0.01
#y <- matrix(rnbinom(400,size=1/dispersion,mu=4),nrow=100)
y <- DGEList(y,group=c(1,1,2,2))
design <- model.matrix(~group, data=y$samples)

#estimate the predictive log fold changes
predlfc<-predFC(y,design,dispersion=dispersion,prior.count=0.5)
logfc <- predFC(y,design,dispersion=dispersion,prior.count=0)
logfc.truncated <- pmax(pmin(logfc,100),-100)

#plot predFC's vs logFC's
plot(predlfc[,2],logfc.truncated[,2],xlab="Predictive log fold changes",ylab="Raw log fold changes")
abline(a=0,b=1)

#####Using new bioconductor package: csaw 

setwd("/Users/maggiec/Documents/Bosselut")
targets <- readTargets()
#bam.files=c(targets$Bam[2],targets$Bam[3],targets$Bam[5],targets$Bam[6])
bam.files=targets$Bam

#Sample:
require(csaw) 
data<-windowCounts(bam.files,ext=110) 
binned<-windowCounts(bam.files,bin=TRUE,width=10000) 
normfacs<-normalizeChIP(binned$counts,lib.size=binned$totals) 
require(edgeR) 
y<-DGEList(data$counts,lib.size=data$totals,norm.factors=normfacs) 
y<-estimateDisp(y,design) 
results<-glmQLFTest(y,design,robust=TRUE) 
merged<-mergeWindows(data$region,tol=1000L) 
tabcom<-combineTests(merged$id,results$table)

##########Cross-correlation plots for fragment length est ########

n<-10000 
SM_WT<-correlateReads("Sample_SM_WT_2.sort.bam",n,dedup=TRUE,cross=TRUE) 
OT_WT<-correlateReads("Sample_OTII_Cont.sort.bam",n,dedup=TRUE,cross=TRUE) 
SM_KO<-correlateReads("Sample_SM_DKO_2.sort.bam",n,dedup=TRUE,cross=TRUE) 
OT_KO<-correlateReads("Sample_OTII_DKO.sort.bam",n,dedup=TRUE,cross=TRUE) 

plot(0:n,SM_WT,col="blue",ylim=c(0,0.1),xlim=c(0,1000), 
     xlab="Delay (bp)", ylab="CCF", pch=16, cex=0.5) 
points(0:n,OT_WT,col="red",pch=16,cex=0.5) 
points(0:n,SM_KO,col="forestgreen",pch=16,cex=0.5)
points(0:n,OT_KO,col="magenta",pch=16,cex=0.5)

legend("topright",col=c("blue","red","forestgreen","magenta"), 
       c("SM_WT", "OT_WT", "SM_KO","OT_KO"), pch=16)

param<-list(bam.files=bam.files,dedup=TRUE,minq=100,pet="none",ext=frag.len) 
demo<-countWindows(param,filter=50) 
demo<-countWindows(param,bin=TRUE,width=1000)


################ Get Regions from +/- 2KB ###############################

regions=read.delim("mouse_prom_2000_rev.txt",header=FALSE,sep="\t")
#regions=read.delim("mouse_prom_chr_2000.txt",header=FALSE,sep="\t")
prom_regions=GRanges(regions$V2,IRanges(regions$V3,regions$V4))
reg.counts<-regionCounts(bam.files,pet="none",prom_regions)

#############  Get SICER Regions ####################

SM=read.delim("Sample_SM_DKO_2-vs-Sample_SM_WT_2-W200-G200-E100-union.island",header=FALSE,sep="\t")
SM_peaks=GRanges(SM$V1,IRanges(SM$V2,SM$V3))

OTII=read.delim("Sample_OTII_DKO-vs-Sample_OTII_Cont-W200-G200-E100-union.island",header=FALSE,sep="\t")
OTII_peaks=GRanges(OTII$V1,IRanges(OTII$V2,OTII$V3))

peak_insct=findOverlaps(SM_peaks,OTII_peaks)
peak_SMOT=subsetByOverlaps(SM_peaks,OTII_peaks , type = "any",ignore.strand=TRUE)

promoter_peaks=subsetByOverlaps(prom_regions,peak_SMOT , type = "any",ignore.strand=TRUE)
reg.counts.peaks<-regionCounts(bam.files,pet="none",promoter_peaks)

#############Try to get the most abundant count regions first ##########

frag.len <- 150
#bin.size <- 2000L
binned <- windowCounts(bam.files, ext=frag.len, minq=100, width=200)
#binned <- windowCounts(bam.files, bin=TRUE, width=bin.size)
normfacs <- normalizeChIP(binned$counts, lib.size=binned$totals)
bin.ab <- aveLogCPM(binned$counts, lib.size=binned$total*normfacs)
threshold <- median(bin.ab)

width <- median(width(binned$region))
eff.win.size <- width + 150
adjustment <- log2(bin.size/eff.win.size)
threshold <- threshold - adjustment

log.min.fc <- log2(10)
threshold <- threshold+log.min.fc

keep <- binned$count >= threshold
hist(bin.ab-adjustment, xlab="Adjusted bin log-CPM", breaks=100, main="")
abline(v=threshold, col="red")

suppressWarnings(keep2 <- overlapsAny(binned$region, prom_regions))
suppressWarnings(keep3 <- overlapsAny(binned$region, merged_peaks))

suppressWarnings(keep4 <- overlapsAny(binned$region, promoter_peaks))


############### Draw plots ######################

ab<-aveLogCPM(reg.counts$counts,lib.size=reg.counts$total)
keep<-ab<=quantile(ab,p=0.9)
normfacs=normalizeChIP(reg.counts$counts[keep,],lib.size=reg.counts$total)

adj.counts<-cpm(reg.counts$counts,log=TRUE)
par(mfrow=c(1,3),mar=c(5,4,2,1.5))

pdfPath="raw_ma.pdf"
pdf(pdfPath)
for(i in 1:(length(bam.files)-1)){
  cur.x <- adj.counts[,1] 
  cur.y <- adj.counts[,1+i] 
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y, 
      xlab="A", ylab="M", main=paste("1 vs", i+1)) 
      all.dist <- diff(log2(normfacs[c(i+1, 1)])) 
      abline(h=all.dist, col="red") 
  }

####  Getting Normalization factors (TMM)
win.counts<-windowCounts(bam.files,bin=TRUE,width=4000L) 
ab<-aveLogCPM(win.counts$counts,lib.sizes=win.counts$total) 
keep<-rank(ab)>0.99*length(ab) 
win.counts.norm<-normalizeChIP(win.counts$count[keep,],lib.sizes=win.counts$totals) 
win.counts.norm

####### Applying Loess to SM data ###################

#> bam.files
#[1] "Sample_SM_WT_2.sort.bam"   
#[2] "Sample_OTII_Cont.sort.bam" 
#[3] "Sample_SM_DKO_2.sort.bam" 
#[4] "Sample_OTII_DKO.sort.bam" 

#####use window counts, then overlap with promoter regions - unbiased way of finding peaks
#win.counts<-windowCounts(win.samples,bin=TRUE,width=10000L)
#prom_regions=GRanges(regions$V2,IRanges(regions$V3,regions$V4))
#suppressWarnings(keep<-overlapsAny(win.counts$region,prom_regions))
#sum(keep)  # like length(keep[keep==TRUE])

###Alternate: count within regions themselves, same as before
win.counts<-regionCounts(bam.files, pet="none", prom_regions)  # Do all files together
#win.counts<-regionCounts(win.samples,pet="none",prom_regions)

keep<-rank(ab)>0.95*length(ab) 
#win.off<-normalizeChIP(win.counts$counts[keep,],lib.size=win.counts$totals,type="loess") 
win.off<-normalizeChIP(reg.counts$counts[keep,],lib.size=reg.counts$totals,type="loess") 
win.off<-normalizeChIP(win.counts$counts,lib.size=win.counts$totals,type="loess") 
head(win.off)

par(mfrow=c(1,2)) 
#aval<-ab[keep] 
ab<-aveLogCPM(cbind(win.counts$counts[,3],win.counts$counts[,1]),
              lib.size=c(win.counts$total[3],win.counts$total[1]))
ab<-aveLogCPM(cbind(win.counts$counts[,4],win.counts$counts[,2]),
              lib.size=c(win.counts$total[4],win.counts$total[2]))
aval<-ab 
o<-order(aval) 
#adjc<-cpm(win.counts$counts[keep,],log=TRUE,lib.size=win.counts$total) 
adjc<-cpm(cbind(win.counts$counts[,3],win.counts$counts[,1]),
          log=TRUE,lib.size=c(win.counts$total[3],win.counts$total[1])) 
adjc<-cpm(cbind(win.counts$counts[,4],win.counts$counts[,2]),
          log=TRUE,lib.size=c(win.counts$total[4],win.counts$total[2])) 
mval<-adjc[,1]-adjc[,2] 
fit<-loessFit(x=aval,y=mval) 
smoothScatter(aval,mval,ylab="M",xlab="AveragelogCPM", 
    main="Raw", ylim=c(-2,2), xlim=c(0, 7)) 
lines(aval[o],fit$fitted[o],col="red") 
#re.adjc<-log2(win.counts$counts[keep,]+0.5)-win.off/log(2)
re.adjc<-log2(win.counts$counts+0.5)-win.off/log(2)
mval<-re.adjc[,1]-re.adjc[,2] 
fit<-loessFit(x=aval,y=mval) 
#smoothScatter(ab[keep],re.adjc[,1]-re.adjc[,2],ylab="M", 
#    xlab="Average logCPM", main="Normalized", ylim=c(-2,2), xlim=c(0, 7)) 
smoothScatter(ab,re.adjc[,1]-re.adjc[,2],ylab="M", 
    xlab="Average logCPM", main="Normalized", ylim=c(-2,2), xlim=c(0, 7)) 
lines(aval[o],fit$fitted[o],col="red")

#######  Run Differential Expression #############
#> targets
#Bam CellType Group
#1   Sample_SM_WT_2.sort.bam       WT     1
#2 Sample_OTII_Cont.sort.bam       WT     2
#3  Sample_SM_DKO_2.sort.bam       KO     1
#4  Sample_OTII_DKO.sort.bam       KO     2

Group <- factor(targets$Group, levels = c(1,2))
CellType <- factor(targets$CellType, levels = c("WT","KO"))
design <- model.matrix(~Group+CellType, data=y$samples)
design

y <- DGEList(win.counts$counts,lib.size=win.counts$totals)
y$offset <- win.off
y <- calcNormFactors(y)

y<-estimateDisp(y,design) 
results<-glmQLFTest(y,design,robust=TRUE) 
#merged<-mergeWindows(data$region,tol=1000L) 
#tabcom<-combineTests(merged$id,results$table)

names(prom_regions)<- regions$V1
strand(prom_regions)<-regions$V5

colnames(results$fitted.values) = paste(targets$Bam,"fitted_values",sep="_")
colnames(y$counts) = paste(targets$Bam,"rawcounts",sep="_")
colnames(results$offset) = paste(targets$Bam,"offsets",sep="_")

finalresults<-cbind(names(prom_regions),regions$V6,results$fitted.values,results$table,y$counts,results$offset)
write.table(finalresults,"ChipSeq_results_2000.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")



####Final Analysis for all H3K27me3 files#####

#First run promoter regions:

perl promoter_parse.pl 
sort -k2n -k3n mouse_prom_3000.txt -o mouse_prom_3000.txt 

##Run in biowulf:
rm(list=ls())

library(csaw)
library(edgeR)
library(GenomicRanges)

targets <- readTargets()
bam.files=targets$Bam
regions=read.delim("mouse_prom_2000_rev.txt",header=FALSE,sep="\t")
prom_regions=GRanges(regions$V2,IRanges(regions$V3,regions$V4))
reg.counts<-regionCounts(bam.files,param=readParam(pet="none"),prom_regions)

####Run on laptop##################
rm(list=ls())
setwd("/data/maggiec/ChipSeq/fastq")
#load("chip.RData")
load("newChip.RData")
#>ls()
#[1] "reg.counts" "regions"    "targets"   

#>targets:
Bam        Genotype Subset
1  Sample_WT_H3K27_DP_2x.sort.bam       wild-type     DP     
2     Sample_WT_H3K27_DP.sort.bam*       wild-type     DP
3      Sample_WT_DP2_K27.sort.bam*       wild-type     DP
4     Sample_KO_H3K27_DP.sort.bam        Jmjd3 ko     DP
5      Sample_UTX_DP_K27.sort.bam*          Utx ko     DP
6   Sample_UTX_DP_K27_2x.sort.bam*          Utx ko     DP
7      Sample_DKO_DP_K27.sort.bam             dko     DP
8          Jmjd3KO_H3K27.sort.bam       wild-type CD4_SP
9          Jmjd3WT_H3K27.sort.bam        Jmjd3 ko CD4_SP
10    Sample_UTX_4SP_K27.sort.bam          Utx ko CD4_SP
11    Sample_DKO_4SP_K27.sort.bam             dko CD4_SP
12            Sample_DKO.sort.bam             dko CD4_SP
13        Sample_SM_WT_2.sort.bam OT-II wild-type CD4_SP
14      Sample_OTII_Cont.sort.bam OT-II wild-type CD4_SP
15       Sample_SM_DKO_2.sort.bam       OT-II dKO CD4_SP
16       Sample_OTII_DKO.sort.bam       OT-II dKO CD4_SP

$totals
[1]  4910484 32461694 24547758 26111731  7970835 16715574 19240687 22288059 20581212 24041029 22356572
[12] 20123858 29302201 20938289 27538979 22020639

#########################################################################

Draw Raw MA plots to look for outliers

########################################################################
rm(list=ls())
load("newChip.RData")


#ab<-aveLogCPM(reg.counts$counts,lib.size=reg.counts$total)
#keep<-ab<=quantile(ab,p=0.9)
#normfacs=normalizeChIP(reg.counts$counts[keep,],lib.size=reg.counts$total)

ab <- aveLogCPM(asDGEList(reg.counts))
keep <- ab <= quantile(ab, p=0.9)
normalize(reg.counts[keep,])
normfacs=normalize(reg.counts[keep,])

#adj.counts<-cpm(reg.counts,log=TRUE)
adj.counts<-cpm(asDGEList(reg.counts),log=TRUE)
par(mfrow=c(1,3),mar=c(5,4,2,1.5))

pdfPath="raw_ma_rev.pdf"
pdfPath="raw_ma2_rev.pdf"
pdfPath="raw_ma3_rev.pdf"
pdf(pdfPath)
j=1
for(i in j:(length(bam.files)-1)){
  cur.x <- adj.counts[,j] 
  cur.y <- adj.counts[,1+i] 
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y, 
                xlab="A", ylab="M", main=paste(j, "vs", i+1)) 
  all.dist <- diff(log2(normfacs[c(i+1, j)])) 
  abline(h=all.dist, col="red") 
}
dev.off()

#########################################

Remove outlier samples - 2,3,5,6 (for simplification, removed 1-7)

Get loess corrected values and plot 

##########################################

k=1
#win.off<-normalizeChIP(reg.counts$counts,lib.size=reg.counts$totals,type="loess") 
#win.off<-normalizeChIP(reg.counts$counts[,k:16],lib.size=reg.counts$totals[k:16],type="loess",iterations = 3) 
win.off<-normalize(reg.counts[,k:16],type="loess",iterations = 3) 
head(win.off)

#cpmnorm <- log2(reg.counts$counts[,k:16]+0.5) - win.off/log(2)
cpmnorm <- log2(assay(reg.counts[,k:16])+0.5) - win.off/log(2)
#cpmnorm = cpm(results$fitted.values,log=TRUE)
#cpmnorm=cpmnorm-results$offset/log(2)

par(mfrow=c(1,3),mar=c(5,4,2,1.5))
pdfPath="norm_ma_all_new_rev.pdf"
pdf(pdfPath)
j=8
num=8
for(i in j:(length(targets[,1]))){
  cur.x <- cpmnorm[,num] 
  cur.y <- cpmnorm[,i] 
  smoothScatter(x=(cur.x+cur.y)/2, y=cur.x-cur.y, 
                xlab="A", ylab="M", main=paste(num, "vs", i)) 
  #smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y, 
  #              xlab="A", ylab="M", main=paste(num, "vs", i+num)) 
  #smoothScatter((normed[,8]+normed[,13])/2, normed[,8]-normed[,13])
  #abline(h=0, col="red")  
  all.dist <- 0
  abline(h=all.dist, col="red") 
}
dev.off()

#############  Get SICER Regions ####################

#SM=read.delim("Sample_SM_DKO_2-vs-Sample_SM_WT_2-W200-G200-E100-union.island",header=FALSE,sep="\t")
SM=read.delim("SM/SM.island",header=FALSE,sep="\t") #import file after "chr" removed
SM_peaks=GRanges(SM$V1,IRanges(SM$V2,SM$V3))

#OTII=read.delim("Sample_OTII_DKO-vs-Sample_OTII_Cont-W200-G200-E100-union.island",header=FALSE,sep="\t")
OTII=read.delim("OTII/OTII.island",header=FALSE,sep="\t") #import file after "chr" removed
OTII_peaks=GRanges(OTII$V1,IRanges(OTII$V2,OTII$V3))

peak_insct=findOverlaps(SM_peaks,OTII_peaks)
combined_peaks = cbind(SM[queryHits(peak_insct),], OTII[subjectHits(peak_insct),])
#peak_SMOT=overlapsAny(SM_peaks,OTII_peaks , type = "any",ignore.strand=TRUE)

comb_min=pmin(combined_peaks[,2],combined_peaks[,5])
comb_max=pmax(combined_peaks[,3],combined_peaks[,6])
merged_peaks=cbind(as.character(combined_peaks[,1]),comb_min,comb_max)

final_merged_peaks=GRanges(merged_peaks[,1],
    IRanges(as.numeric(merged_peaks[,2]),as.numeric(merged_peaks[,3])))

##Check:
final_merged_peaks[seqnames(final_merged_peaks) == "X" 
                   & start(final_merged_peaks) > 17720996 
                   & end(final_merged_peaks) < 18725818,]


###Add the DP peaks as well:
DP=read.delim("DP/DP.island",header=FALSE,sep="\t") 
DP_peaks=GRanges(DP$V1,IRanges(DP$V2,DP$V3))

peak_union=union(final_merged_peaks,DP_peaks)

#dedup_final_merged_peaks=final_merged_peaks[duplicated(final_merged_peaks)==FALSE]
#dedup_prom_regions=prom_regions[duplicated(prom_regions)==FALSE]

dedup_final_merged_peaks=peak_union[duplicated(peak_union)==FALSE]
dedup_prom_regions=prom_regions[duplicated(prom_regions)==FALSE]

#promoter_peaks_SM=subsetByOverlaps(prom_regions,SM_peaks , type = "any",ignore.strand=TRUE)
promoter_peaks_OT=subsetByOverlaps(prom_regions,OTII_peaks , type = "any",ignore.strand=TRUE)

peak_insct2=findOverlaps(dedup_prom_regions,dedup_final_merged_peaks)
promoter_peaks=subsetByOverlaps(dedup_prom_regions,dedup_final_merged_peaks , type = "any",ignore.strand=TRUE)

#reg.counts.peaks<-regionCounts(bam.files,pet="none",promoter_peaks)
####################################################

#####Isolate promoter regions with SICER peaks

#####################################################

#> dim(reg.counts$counts)
#[1] 32800    16
#> length(prom_regions)
#[1] 32800
#> length(dedup_prom_regions)
#[1] 25654
#> length(promoter_peaks)
#[1] 2544

#prom_regions[subsetByOverlaps(prom_regions,promoter_peaks , type="equal")==TRUE]
promoter=findOverlaps(prom_regions,promoter_peaks, type="equal")
promoter_OT=findOverlaps(prom_regions,promoter_peaks_OT, type="equal")

#reg.counts$counts2=reg.counts$counts[queryHits(promoter),]
reg.counts2=reg.counts[queryHits(promoter),]

#############################################

#####Analysis using edgeR

##################################################################
#k=1
k=13
#y <- DGEList(reg.counts$counts,lib.size=reg.counts$totals)
#y <- DGEList(reg.counts$counts2[,k:16],lib.size=reg.counts$totals[k:16])
#y <- DGEList(cpmnorm)

#Genotype <- factor(targets$Genotype, levels = c("wild-type","Jmjd3 ko", "Utx ko", "dko", "OT-II wild-type", "OT-II dKO"))
#Genotype <- factor(targets$Genotype[k:16], levels = c("wild-type","Jmjd3 ko", "Utx ko", "dko", "OT-II wild-type", "OT-II dKO"))
Genotype <- factor(targets$Genotype[k:16], levels = c("OT-II dKO","OT-II wild-type"))
#Group <- factor(c(1,2,1,2), levels=c("1","2"))
#Subset <- factor(targets$Subset, levels = c("DP","CD4_SP"))
design <- model.matrix(~Genotype, data=y$samples)
#design <- model.matrix(~Group+Genotype, data=y$samples)
#design <- model.matrix(~Subset+Genotype, data=y$samples)

design

#y$offset <- win.off[queryHits(promoter),]
y$offset <- win.off[queryHits(promoter),k:16]
y <- calcNormFactors(y)

y<-estimateDisp(y,design) 
results<-glmQLFTest(y,design,robust=TRUE) 

#y <- asDGEList(reg.counts2[,k:16])
y <- asDGEList(reg.counts2[,k:16], norm.factors=normfacs[k:16])
y <- estimateDisp(y, design)
y$offset <- win.off[queryHits(promoter),k:16]
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results)

############################################

#Plot density plots

#############################################

pdfPath="density.pdf"
pdf(pdfPath)
plot (density(cpmnorm[,1]), col="red",ylim=c(0,0.3),xlim=c(-4,12)) 
lines (density(cpmnorm[,2]), col="green")
lines (density(cpmnorm[,3]), col="blue")
lines (density(cpmnorm[,4]), col="orange")
lines (density(cpmnorm[,5]), col="purple")
lines (density(cpmnorm[,6]), col="yellow")
lines (density(cpmnorm[,7]), col="black")
lines (density(cpmnorm[,8]), col="grey")
lines (density(cpmnorm[,9]), col="pink")
dev.off()



############################################

names(prom_regions)<- regions$V1
strand(prom_regions)<-regions$V5

#colnames(cpmnorm) = paste(targets$Bam,"norm_cpm",sep="_")
colnames(cpmnorm) = paste(targets$Bam[1:16],"norm_cpm",sep="_")
#colnames(y$counts) = paste(targets$Bam,"rawcounts",sep="_")
#colnames(results$offset) = paste(targets$Bam,"offsets",sep="_")

peaks_OTonly = as.character(regions$V6)[queryHits(promoter_OT)]

as.character(regions$V6)[queryHits(promoter)] %in% peaks_OTonly

compresults=as.data.frame(topTags(results,n=Inf,sort.by="none"))
finalresults<-cbind(names(prom_regions)[queryHits(promoter)],
          as.character(regions$V6)[queryHits(promoter)],
          as.character(regions$V6)[queryHits(promoter)] %in% peaks_OTonly,
          as.character(regions$V5)[queryHits(promoter)],
          cpmnorm[queryHits(promoter),],compresults)
colnames(finalresults)[1] = "RefSeq"
colnames(finalresults)[2] = "Genes"
colnames(finalresults)[3] = "OTII peak"
colnames(finalresults)[4] = "Strand"
finalresults[finalresults$Genes=="S1pr1",]
#NM_007901  S1pr1	8.159113405	8.630411362	8.715935467	
#    8.662441045	8.666934306	8.552539018	9.300885669	
#    8.855297979	6.876174122	7.427865283	10.36109683	
#    8.889909327	6.392206344	6.223978045	9.347964645	
#    8.860431979	
#    2.813657865	3.543396393	86.83832915	2.99E-06
finalresults[finalresults$Genes=="Zbtb7b",]
finalresults[finalresults$Genes=="Prnd",]

#write.table(finalresults,"ChipSeq_DP2_results_2000_loessvals_pvals.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(finalresults,"ChipSeq_DP2_results_2000_loessvals_pvals_rev.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")


cpmdata=cbind(names(prom_regions),as.character(regions$V6),cpmnorm)
colnames(cpmdata)[1]="RefSeq"
colnames(cpmdata)[2]="Gene"
#write.table(cpmdata,"ChipSeq_allcpm_loess.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(cpmdata,"ChipSeq_allcpm_loess_rev.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

###############################################

#write out all offsets:
k=1
y$offset <- win.off[queryHits(promoter),k:16]
colnames(y$offset) = bam.files
finaloffsets<-cbind(names(prom_regions)[queryHits(promoter)],
                    as.character(regions$V6)[queryHits(promoter)],
                    y$offset)

dimnames(finaloffsets)[[2]][1] = "RefSeq"
dimnames(finaloffsets)[[2]][2] = "Genes"
head(finaloffsets)

#write.table(finaloffsets,"ChipSeq_alloffsets_loess.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(finaloffsets,"ChipSeq_alloffsets_loess_rev2.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

###############################################

#write out bedfile formats:

chr=seqnames(prom_regions[queryHits(promoter),])
startend=ranges(prom_regions[queryHits(promoter),])
start=start(startend)-1
end=end(startend)-1
bedfile=cbind(as.data.frame(chr),start,end,as.data.frame(startend)[,4])
colnames(bedfile)[1]<- "chr"

write.table(bedfile,"ChipSeq_bed_promoters_rev.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#############################################################################

#### Find Distribution of Peaks in Promoter vs Intergenic vs Intragenic #####

# Get Sicer Regions

SMWT=read.delim("SM/SM_WT_peaks.bed",header=FALSE,sep="\t") #import file after "chr" removed
SMKO=read.delim("SM/SM_DKO_peaks.bed",header=FALSE,sep="\t") #import file after "chr" removed
  
SM_WT_peaks=GRanges(SMWT$V1,IRanges(SMWT$V2,SMWT$V3))
SM_KO_peaks=GRanges(SMKO$V1,IRanges(SMKO$V2,SMKO$V3))

OTWT=read.delim("OTII/OTII_WT_peaks.bed",header=FALSE,sep="\t") #import file after "chr" removed
OTKO=read.delim("OTII/OTII_DKO_peaks.bed",header=FALSE,sep="\t") #import file after "chr" removed

OT_WT_peaks=GRanges(OTWT$V1,IRanges(OTWT$V2,OTWT$V3))
OT_KO_peaks=GRanges(OTKO$V1,IRanges(OTKO$V2,OTKO$V3))

WT_peak_insct=findOverlaps(SM_WT_peaks,OT_WT_peaks)
KO_peak_insct=findOverlaps(SM_KO_peaks,OT_KO_peaks)

WT_combined_peaks = cbind(SMWT[queryHits(WT_peak_insct),], OTWT[subjectHits(WT_peak_insct),])
KO_combined_peaks = cbind(SMKO[queryHits(KO_peak_insct),], OTKO[subjectHits(KO_peak_insct),])

WT_comb_min=pmin(WT_combined_peaks[,2],WT_combined_peaks[,5])
WT_comb_max=pmax(WT_combined_peaks[,3],WT_combined_peaks[,6])
WT_merged_peaks=cbind(as.character(WT_combined_peaks[,1]),WT_comb_min,WT_comb_max)
WT_final_merged_peaks=GRanges(WT_merged_peaks[,1],
    IRanges(as.numeric(WT_merged_peaks[,2]),as.numeric(WT_merged_peaks[,3])))

KO_comb_min=pmin(KO_combined_peaks[,2],KO_combined_peaks[,5])
KO_comb_max=pmax(KO_combined_peaks[,3],KO_combined_peaks[,6])
KO_merged_peaks=cbind(as.character(KO_combined_peaks[,1]),KO_comb_min,KO_comb_max)
KO_final_merged_peaks=GRanges(KO_merged_peaks[,1],
    IRanges(as.numeric(KO_merged_peaks[,2]),as.numeric(KO_merged_peaks[,3])))

length(WT_final_merged_peaks) # 49394  
length(KO_final_merged_peaks) # 43203

WT_dedup_final_merged_peaks=WT_final_merged_peaks[duplicated(WT_final_merged_peaks)==FALSE]
KO_dedup_final_merged_peaks=KO_final_merged_peaks[duplicated(KO_final_merged_peaks)==FALSE]

length(WT_dedup_final_merged_peaks) # 44430  intersected, more stringent
length(KO_dedup_final_merged_peaks) # 33004

WT_peak_union=union(SM_WT_peaks,OT_WT_peaks)
KO_peak_union=union(SM_KO_peaks,OT_KO_peaks)
length(WT_peak_union)  # 55798 - more but not just overlap but merge (no dups)
length(KO_peak_union)  # 64637

WT_final_merged_peaks[seqnames(WT_final_merged_peaks) == "1" 
                & start(WT_final_merged_peaks) > 3659200 
                & end(WT_final_merged_peaks) < 5000000,]

WT_dedup_final_merged_peaks[seqnames(WT_dedup_final_merged_peaks) == "1" 
                      & start(WT_dedup_final_merged_peaks) > 3659200 
                      & end(WT_dedup_final_merged_peaks) < 5000000,]

WT_peak_union[seqnames(WT_peak_union) == "1" 
                      & start(WT_peak_union) > 3659200 
                      & end(WT_peak_union) < 5000000,]

WT_dedup_peak_union=WT_peak_union[duplicated(WT_peak_union)==FALSE] #Doesn't remove anything

##### Count the overlaps with Promoters, Intragenic and Intergenic Regions: ###

#Get Promoter Peaks:
regions=read.delim("mouse_prom_2000_rev.txt",header=FALSE,sep="\t")
prom_regions=GRanges(regions$V2,IRanges(regions$V3,regions$V4))
dedup_prom_regions=prom_regions[duplicated(prom_regions)==FALSE]

#> length(prom_regions)
#[1] 24014
#> length(dedup_prom_regions)
#[1] 23855

#Count number of promoter peaks that have overlap with either WT or KO peaks:
table(overlapsAny(dedup_prom_regions,WT_dedup_final_merged_peaks))
FALSE  TRUE 
15834  8021 
table(overlapsAny(dedup_prom_regions,KO_dedup_final_merged_peaks))
FALSE  TRUE 
16542  7313

#Count number of WT or KO SICER peaks that have overlap with promoter peaks:
table(overlapsAny(WT_dedup_final_merged_peaks,dedup_prom_regions))
FALSE  TRUE 
34592  9838 
table(overlapsAny(KO_dedup_final_merged_peaks,dedup_prom_regions))
FALSE  TRUE 
23562  9442 

subsetByOverlaps(KO_dedup_final_merged_peaks, dedup_prom_regions)

#### Get mm9 

library(GenomicFeatures)
library(ChIPpeakAnno)

txdb <- makeTranscriptDbFromUCSC(genome="mm9", tablename="refGene") 
exonRanges <- exonsBy(txdb, "tx")
