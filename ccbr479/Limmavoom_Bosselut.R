?
#First make directory with hard links to bam files
#ln {source} {link}
#cd into directory, launch R, paste in the following code:

setwd("/Users/maggiec/GitHub/Maggie/ccbr479/data")
source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")
library(Rsubread)
library(limma)

#For libraries under development:
library(BiocInstaller) 
useDevel()
biocLite("ggplot")
biocLite("Rsubread")
biocLite("EDASeq")
biocLite("RUVSeq")

#bamfiles = list.files(path=".",pattern="bam$", full=TRUE) 
#fc <- featureCounts(bamfiles,annot.inbuilt='mm10')
options(digits=20)
targets <- readTargets()
load("Gencode_fc.RData")

#Ensembl:
gtf <- "/data/maggiec/RNASeq/GTF/Mus_musculus/Ensembl/NCBIM37/Annotation/Archives/archive-2012-03-09-05-33-07/Genes/genes.gtf"
gtf <- "/data/maggiec/RNASeq/GTF/Mus_musculus.GRCm38.77.gtf"
gtf <- "/fdb/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf" 

#RefSeq:
gtf <- "/data/maggiec/RNASeq/GTF/RefSeq_mm10.gtf"
gtf <- "/fdb/igenomes/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf"

#UCSC genes:
gtf <- "/data/maggiec/RNASeq/GTF/UCSC_mm10.gtf"
gtf <- "/fdb/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"

#Gencode genes:
gtf <- "/data/maggiec/RNASeq/Genomes/mm10/gencode.vM4.all.gtf"

#Paired end, stranded data (first ran without paired/ss):
fc <- featureCounts(files=targets$Bam,annot.inbuilt="mm10",nthreads=32)

fc_single <- featureCounts(files=targets$Bam,isGTFAnnotationFile=TRUE,
                           annot.ext=gtf,nthreads=32)

fc_paired_r <- featureCounts(files=targets$Bam,isGTFAnnotationFile=TRUE,nthreads=32,
          annot.ext=gtf,GTF.attrType="gene_name",strandSpecific=2,isPairedEnd=TRUE)

#Unstranded: Nextera is unstranded!
fc <- featureCounts(files=targets$bam,isGTFAnnotationFile=TRUE,nthreads=32,
      annot.ext=gtf,GTF.attrType="gene_name",strandSpecific=0,isPairedEnd=TRUE)
#Forward-stranded:
# fc <- featureCounts(files=targets$bam,isGTFAnnotationFile=TRUE,nthreads=32,
#                     annot.ext=gtf,strandSpecific=1,isPairedEnd=TRUE)
# 
# #Reverse-stranded:
# fc <- featureCounts(files=targets$bam,isGTFAnnotationFile=TRUE,nthreads=32,
#                     annot.ext=gtf,strandSpecific=2,isPairedEnd=TRUE)

## To save:
save(targets,fc, gtf, file="UCSC_fc.RData")
write.table(fc$count,file="UCSC_featurecounts.txt",row.names=FALSE,sep="\t",quote=FALSE)

save(targets,fc, gtf, file="RefSeq_fc.RData")
write.table(fc$count,file="RefSeq_featurecounts.txt",row.names=FALSE,sep="\t",quote=FALSE)

save(targets,fc, gtf, file="Ensembl_fc.RData")
write.table(fc$count,file="Ensembl_featurecounts.txt",row.names=FALSE,sep="\t",quote=FALSE)

save(targets,fc, gtf, file="Gencode_fc.RData")
write.table(fc$count,file="Gencode_featurecounts.txt",row.names=FALSE,sep="\t",quote=FALSE)


###########Load Data from Biowulf##########################
#setwd("/Users/maggiec/Documents/Bosselut2")
setwd("/Volumes/CCRIFXCCR/CCBR-479/Unaligned/bam2/")
load("UCSC_fc.RData")
load("Gencode_fc.RData")
targets <- readTargets()

targets
# bam    Sample  Phenotype   Genotype  Mouse
# 1  Sample_CRE_CTRL1-STAR.fin2.bam CRE_CTRL1    wt YFP+  Foxp3-cre   8619
# 2  Sample_CRE_CTRL2-STAR.fin2.bam CRE_CTRL2    wt YFP+  Foxp3-cre    124
# 3  Sample_CRE_CTRL3-STAR.fin2.bam CRE_CTRL3    wt YFP+  Foxp3-cre    124
# 4  Sample_CRE_CTRL4-STAR.fin2.bam CRE_CTRL4    wt YFP+  Foxp3-cre   8800
# 5  Sample_CRE_CTRL5-STAR.fin2.bam CRE_CTRL5    wt YFP+  Foxp3-cre   8800
# 6  Sample_CRE_CTRL7-STAR.fin2.bam CRE_CTRL7    wt YFP+  Foxp3-cre   8619
## 7       Sample_EXP2-STAR.fin2.bam      EXP2 dbl mutant dbl mutant 511714
## 8       Sample_EXP3-STAR.fin2.bam      EXP3 dbl mutant dbl mutant 132614
# 9       Sample_EXP4-STAR.fin2.bam      EXP4 dbl mutant dbl mutant 232614
# 10      Sample_EXP5-STAR.fin2.bam      EXP5 dbl mutant dbl mutant 232614
# 11      Sample_EXP6-STAR.fin2.bam      EXP6 dbl mutant dbl mutant 432614
### 12      Sample_EXP7-STAR.fin2.bam      EXP7 dbl mutant dbl mutant 332614
# 13      Sample_EXP8-STAR.fin2.bam      EXP8 dbl mutant dbl mutant 142114
# 14      Sample_EXP9-STAR.fin2.bam      EXP9 dbl mutant dbl mutant 142114
# 15 Sample_GEN_CTRL1-STAR.fin2.bam GEN_CTRL1    wt YFP- dbl mutant 142114
# 16 Sample_GEN_CTRL2-STAR.fin2.bam GEN_CTRL2    wt YFP- dbl mutant 142114
# 17 Sample_GEN_CTRL3-STAR.fin2.bam GEN_CTRL3    wt YFP- dbl mutant 232614
# 18 Sample_GEN_CTRL4-STAR.fin2.bam GEN_CTRL4    wt YFP- dbl mutant 232614
# 19 Sample_GEN_CTRL5-STAR.fin2.bam GEN_CTRL5    wt YFP- dbl mutant 432614
# 20 Sample_GEN_CTRL6-STAR.fin2.bam GEN_CTRL6    wt YFP- dbl mutant 432614
## 21 Sample_GEN_CTRL7-STAR.fin2.bam GEN_CTRL7    wt YFP- dbl mutant 332614
#^^ 22 Sample_GEN_CTRL8-STAR.fin2.bam GEN_CTRL8    wt YFP- dbl mutant 332614 # replace with below
# 22 Sample_CTRL8-STAR.fin2.bam GEN_CTRL8 dbl mutant  dbl mutant  332614
## 23 Sample_THP_CTRL1-STAR.fin2.bam THP_CTRL1    Thpok-d    Thpok-d   8799
## 24 Sample_THP_CTRL2-STAR.fin2.bam THP_CTRL2    Thpok-d    Thpok-d   8799

library(edgeR)
options(digits=2)
x <- DGEList(counts=fc$counts, genes=fc$annotation)

##########Remove Outlier #7#####################
#x <- DGEList(counts=fc$counts[,-12], genes=fc$annotation[,1:6])

targets$Group[targets$Phenotype=="dbl mutant"]="Thpok_LRF_min"
targets$Group[targets$Phenotype=="wt YFP+"]="Ctrl"
targets$Group[targets$Phenotype=="wt YFP-"]="Thpok_LRF_plus"
targets$Group[targets$Phenotype=="Thpok-d"]="Thpok_d"
targets$Group[targets$Sample=="GEN_CTRL8"]="Thpok_LRF_min"
targets$Genotype[targets$Sample=="GEN_CTRL8"]="dbl mutant"

###First outlier, but don't run:
Pheno <- factor(targets$Phenotype[-12])
Geno <- factor(targets$Genotype[-12])
Mouse <- factor(targets$Mouse[-12])
Group <- factor(targets$Group[-12])

table(Mouse,Group)  ## Look at technical reps:
# Group
# Mouse    Ctrl Thpok_d Thpok_LRF_min Thpok_LRF_plus
# 124       2       0             0              0
# 8619      2       0             0              0
# 8799      0       2             0              0
# 8800      2       0             0              0
# 132614    0       0             1              0
# 142114    0       0             2              2
# 232614    0       0             2              2
# 332614    0       0             1              1
# 432614    0       0             1              2
# 511714    0       0             1              0
# Make 332614 a technical replicate of 432614 for Thpok_LRF_min

targets1=targets[c(-1,-2,-3,-4,-5,-6,-7,-8,-12,-21,-23,-24),]
targets1$Mouse[targets1$Sample=="GEN_CTRL8"]="432614"
targets1$Sample[targets1$Sample=="GEN_CTRL8"]="EXP_new"

####Use this one instead to remove outliers, leaves 6 per group, and relabel:
targets2=targets[c(-7,-8,-12,-21,-23,-24),]
targets2$Mouse[targets2$Sample=="GEN_CTRL8"]="432614"
targets2$Sample[targets2$Sample=="GEN_CTRL8"]="EXP_new"
targets2

Mouse2 <- factor(targets2$Mouse)
Group2 <- factor(targets2$Group)
#Mouse2[18] = "432614"
table(Mouse2,Group2) 

design=model.matrix(~0+Group2)
colnames(design) <- levels(Group2)
cont.matrix <- makeContrasts(Thpok_LRF_min_Ctrl=Thpok_LRF_min-Ctrl,
                             Thpok_LRF_plus_Ctrl=Thpok_LRF_plus-Ctrl,
                             Thpok_LRF_min_plus=Thpok_LRF_min-Thpok_LRF_plus, 
                             Diff=(Thpok_LRF_min-Thpok_LRF_plus)-(Thpok_LRF_min-Ctrl),
                             levels=design)
# Group
# Mouse    Ctrl Thpok_LRF_min Thpok_LRF_plus
# 124       2             0              0
# 142114    0             2              2
# 232614    0             2              2
# 432614    0             2              2
# 8619      2             0              0
# 8800      2             0              0

# table(Mouse,Group)
# Group
# Mouse    Thpok_LRF_min Thpok_LRF_plus
# 142114             2              2
# 232614             2              2
# 432614             2              2

x <- DGEList(counts=fc$counts[,c(-7,-8,-12,-21,-23,-24)], genes=fc$annotation)
isexpr <- rowSums(cpm(x)>0.5) >= 2
hasannot <- rowSums(is.na(x$genes))==0
x <- x[isexpr & hasannot,keep.lib.sizes=FALSE]
x <- calcNormFactors(x)

#fc1=fc$counts[,c(-7,-8,-12,-21,-23,-24)]
#tfc1=t(fc1)
#filter <- apply(fc1, 1, function(x) length(x[x>5])>=2)
#fc1filt <- fc1[filter,]
#genes <- rownames(fc1filt)

#v1 <- voom(as.matrix(fc1filt),design,plot=TRUE,normalize="quantile")

#dim(v1)
#[1] 11947    18
#[1] 15111    18  (including new lncRNA)

v <- voom(x,design,plot=TRUE,normalize="quantile")
dim(v)

library(ggplot2)
library(reshape)
library(plotly)

fc1filt=as.data.frame(fc1filt)
fc1filtnames=sapply(strsplit(colnames(fc1filt),split="\\."),function(x) x[1])
colnames(fc1filt)=paste(targets2$Sample,targets2$Phenotype)

df.m <- melt(as.data.frame(fc1filt))
ggplot(df.m) + 
  geom_density(aes(x = value, colour = variable)) + labs(x = NULL) +
  theme(legend.position='right') + scale_x_log10() + ggtitle("Raw Counts")

par(mar=c(10,7,1,1))
boxplot(log(value)~variable,las=2,data=df.m,main="Raw Signal", 
        ylab="Counts",col=as.numeric(as.factor(targets2$Sample)))


edf=as.matrix(v$E)
tedf= t(edf)
pca=prcomp(tedf,scale.=T)
tedf1 = data.frame(tedf)
Phenotype=targets2$Group
#Phenotype=sapply(strsplit(targets2$Sample,split="_"),function(x) x[1])
cell_rep=paste(Phenotype,targets2$Mouse)
tedf1$group = as.factor(Phenotype)

#plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
#pca$x[,1:3]  #look at first 3 PC's

plot3d(pca$x[,1:3],col = as.integer(tedf1$group),type="s",size=2)
group.v<-as.vector(cell_rep)
text3d(pca$x, pca$y, pca$z, group.v, cex=1.0, adj = 1.2) 


##EDASeq 
x <- as.factor(sumfc[,2])
set <- newSeqExpressionSet(as.matrix(filtered),
        phenoData = data.frame(x, row.names=colnames(filtered)))


colnames(v1Efilt)=colnames(v1$E)
jpeg(filename="MA_plot.jpeg")
par(mfrow=c(2,3))
l = dim(v1Efilt)[2]-1
for (i in 1:l) {
  j = i+1
  MA = cbind(v1Efilt[,i],v1Efilt[,j])
  colnames(MA) = c(colnames(v1Efilt[,i]),colnames(v1Efilt[,j]))
  plotMA(MA,main=colnames(MA))
  abline(h=0,col="red")
}
dev.off()

######Aggregate the data (mean or sum)
##Decided against using average, use sum instead, below)

fc1=fc$counts[,c(-1,-2,-3,-4,-5,-6,-7,-8,-12,-21,-23,-24)]
tfc1=t(fc1)
meanfc=aggregate(tfc1,list(Mouse,Group),FUN="mean")

meanfc[,1:2]
# Group.1        Group.2
# 1      124           Ctrl
# 2     8619           Ctrl
# 3     8800           Ctrl
## 4     8799        Thpok_d
# 5   132614  Thpok_LRF_min
# 6   142114  Thpok_LRF_min
# 7   232614  Thpok_LRF_min
# 8   432614  Thpok_LRF_min
# 9   511714  Thpok_LRF_min
# 10  142114 Thpok_LRF_plus
# 11  232614 Thpok_LRF_plus
## 12  332614 Thpok_LRF_plus
# 13  432614 Thpok_LRF_plus

Group=as.character(meanfc[c(-4,-8,-12),2])
Mouse=as.character(meanfc[c(-4,-8,-12),1])
Group=as.factor(Group)
Mouse=as.factor(Mouse)

##Try Quantile
tmeanfc=t(meanfc[c(-4,-8,-12),3:23420])
colnames(tmeanfc)=Group
v1 <- voom(as.matrix(tmeanfc),design,plot=TRUE,normalize="quantile")


#########Get Sum instead #######

tfc1=t(x$counts)
sumfc=aggregate(tfc1,list(Mouse2,Group2),FUN="sum")
sumfc[,1:2]

# Group.1        Group.2
# 1     124           Ctrl
# 2    8619           Ctrl
# 3    8800           Ctrl
# 4  142114  Thpok_LRF_min
# 5  232614  Thpok_LRF_min
# 6  432614  Thpok_LRF_min
# 7  142114 Thpok_LRF_plus
# 8  232614 Thpok_LRF_plus
# 9  432614 Thpok_LRF_plus
Group=as.factor(sumfc[4:9,2])
Mouse=as.factor(sumfc[4:9,1])

Group=as.factor(sumfc[1:9,2])  # Redo, include Ctrl
Mouse=as.factor(sumfc[1:9,1])

sumfc[,1:2]
#Group.1        Group.2
# 1  142114  Thpok_LRF_min
# 2  232614  Thpok_LRF_min
# 3  432614  Thpok_LRF_min
# 4  142114 Thpok_LRF_plus
# 5  232614 Thpok_LRF_plus
# 6  432614 Thpok_LRF_plus

Group=as.factor(sumfc[,2])
Mouse=as.factor(sumfc[,1])

dim(sumfc)
#[1]     9 43282
#tsumfc=t(sumfc[,3:23422])
#tsumfc=t(sumfc[,3:43282])
tsumfc=t(sumfc[,3:13754])
dim(tsumfc)
#tsumfc=tsumfc[3:23422,]
head(tsumfc)
colnames(tsumfc)=paste(sumfc[,1],sumfc[,2],sep=".")
head(tsumfc)
filter <- apply(tsumfc, 1, function(x) length(x[x>5])>=2)
sumfcfilt <- tsumfc[filter,]
head(sumfcfilt)
dim(sumfcfilt)
genes <- rownames(sumfcfilt)

####Histogram: Look at distribution of data before quantile:
library(reshape)
df.m <- melt(as.data.frame(sumfcfilt))
ggplot(df.m) + geom_density(aes(x = value,
  colour = variable)) + labs(x = NULL) +
  theme(legend.position='top') + scale_x_log10()

###Do Normalization: cpm is further normalized by quantile
##v1 = scaled by lib.size then quantile norm

design=model.matrix(~0+Group)
colnames(design) <- levels(Group)
cont.matrix <- makeContrasts(Thpok_LRF_min_Ctrl=Thpok_LRF_min-Ctrl,
                             Thpok_LRF_plus_Ctrl=Thpok_LRF_plus-Ctrl,             
                             Thpok_LRF_min_plus=Thpok_LRF_min-Thpok_LRF_plus, 
                             Diff=(Thpok_LRF_min-Thpok_LRF_plus)-(Thpok_LRF_min-Ctrl),
                             levels=design)
v1 <- voom(sumfcfilt,design,plot=TRUE,normalize="quantile")

#####First set of Scaling factors - changed to actscale (down below):
head(log2((sumfcfilt/colSums(sumfcfilt))*1000000))

sf = v1$E/log2((sumfcfilt/colSums(sumfcfilt))*1000000)
head(sf)  #Check scaling factors - correct version (first scaled by sum of assigned counts)
write.table(sf,file="scaling_factors.txt",
            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

#####Tried scaling first by total library counts - offset doesn't work as well)
lmfitcoeff=vector()
scale=matrix(nrow=dim(v1$E)[1],ncol=dim(v1$E)[2])
scaleoff=matrix(nrow=dim(v1$E)[1],ncol=dim(v1$E)[2])
for (i in 1:9){
x=log2(sumfcfilt[,i]/v1$targets[i,1])
y=(v1$E[,i])
ylim=y[is.infinite(x)==FALSE]
xlim=x[is.infinite(x)==FALSE]
#plot(xlim,ylim)

lmfit=lm(ylim~xlim)
abline(lmfit,col="red")
lmfitcoeff[i]=lmfit$coefficients[1]
scaleoff[,i]=2^(y-lmfitcoeff[i])/2^x
scale[,i]=2^y/2^x
}

scaled_dat=(sumfcfilt)*(scaleoff)
scaleval=(2^v1$E)*100  # Real data (aiming for)
actscale=(sumfcfilt/scaleval)

scale_fit=(sumfcfilt/actscale) ###Need to use this formula to calculate normalized counts
scale_fit2=(sumfcfilt/sf)

scaled_dat[rownames(scaled_dat)=="Zbtb7b"]
scaleval[rownames(scaleval)=="Zbtb7b"]
scale_fit[rownames(scale_fit)=="Zbtb7b"]
scale_fit2[rownames(scale_fit2)=="Zbtb7b"]
sumfcfilt[rownames(sumfcfilt)=="Zbtb7b"]

plot(scaleval[,1],scale_fit2[,1])  # Closer to scaleval
plot(scaleval[,1],scale_fit[,1])

df.m <- melt(as.data.frame(sf))  # look at original sf calculated using sumfcfilt counts 
df.m <- melt(as.data.frame(actscale))  # look at actual scaling factor
ggplot(df.m) +  xlim(0,2) + geom_density(aes(x = value,
  colour = variable)) 
boxplot(value~variable,data=df.m)

#### Use scale factor (actscale to check if Quantile Normalized)

df.m <- melt(as.data.frame(scale_fit)) 
ggplot(df.m) + geom_density(aes(x = value,
  colour = variable)) + labs(x = NULL) +
  theme(legend.position='top') + scale_x_log10()
boxplot(value~variable,data=df.m)

write.table(actscale,file="scaling_factors2.txt",
            row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)

#####Look after quantile normalization:

head(v1$weights)
df.n <- melt(as.data.frame(v1$E))
ggplot(df.n) + geom_density(aes(x = value,
    colour = variable)) + labs(x = NULL) +
  theme(legend.position='top')

edf=as.matrix(v1$E)
tedf= t(edf)
pca=prcomp(tedf,scale.=T)
tedf1 = data.frame(tedf)
Phenotype=sumfc[,2]
#Phenotype=sapply(strsplit(targets2$Sample,split="_"),function(x) x[1])
cell_rep=paste(Phenotype,sumfc[,1],sep=".")
tedf1$group = as.factor(Phenotype)

#plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
#pca$x[,1:3]  #look at first 3 PC's

plot3d(pca$x[,1:3],col = as.integer(tedf1$group),type="s",size=2)
group.v<-as.vector(cell_rep)
text3d(pca$x, pca$y, pca$z, group.v, cex=1.0, adj = 1.2) 

#####Run GLM Model##########

design=model.matrix(~0+Group)
colnames(design) <- levels(Group)
cont.matrix <- makeContrasts(Thpok_LRF_min_Ctrl=Thpok_LRF_min-Ctrl,
                             Thpok_LRF_plus_Ctrl=Thpok_LRF_plus-Ctrl,             
                             Thpok_LRF_min_plus=Thpok_LRF_min-Thpok_LRF_plus, 
                             Diff=(Thpok_LRF_min-Thpok_LRF_plus)-(Thpok_LRF_min-Ctrl),
                             levels=design)

fit <- lmFit(v1,design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2,coef=ncol(design))
results <- decideTests(fit2)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

finalres=topTable(fit2,sort="none",n=Inf)

FC = 2^finalres[,1:4]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC) = paste(colnames(FC),"FC",sep=".")
Thpok_LRF_min_Ctrl_pval=p.adjust(fit2$p.value[,1], method='BH')
Thpok_LRF_plus_Ctrl_pval=p.adjust(fit2$p.value[,2], method='BH')
Thpok_LRF_min_plus_pval=p.adjust(fit2$p.value[,3], method='BH')
Diff_pval=p.adjust(fit2$p.value[,4], method='BH')
pvalall=cbind(Thpok_LRF_min_Ctrl_pval,Thpok_LRF_plus_Ctrl_pval,
              Thpok_LRF_min_plus_pval,Diff_pval)
finalres = cbind(v1$E,FC,pvalall)

finalres$logFC = -1*finalres$logFC
finalres[rownames(finalres)=="Zbtb7b",]  #Thpok
finalres[rownames(finalres)=="Foxp3",]  #Foxp3

#write.table(v1$E,file="UCSC_voom_Q_Ctrl.txt",
#            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
#write.table(finalres,file="UCSC_FCpval_Ctrl.txt",
#            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

write.table(v1$E,file="Gencode_voom_Q_Ctrl.txt",
            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(finalres,file="Gencode_FCpval_Ctrl.txt",
            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

####  Try TMM Normalization ##################
x <- calcNormFactors(x)
x$sample
v <- voom(x,design,plot=TRUE)
write.table(v$E,file="UCSC_voom.txt",
    row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

#####Try Normalization method########

library("RUVSeq")
dim(tsumfc)
filter <- apply(tsumfc, 1, function(x) length(x[x>5])>=2)
filtered <- tsumfc[filter,]
genes <- rownames(filtered)

x <- as.factor(sumfc[,2])
set <- newSeqExpressionSet(as.matrix(filtered),
    phenoData = data.frame(x, row.names=colnames(filtered)))

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x],  las = 2, cex.lab = 0.8)
plotPCA(set, col=colors[x], cex=1.2)

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], las = 2, cex.lab = 0.8)
plotPCA(set, col=colors[x], cex=1.2)

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)

# x         W_1
# 124.Ctrl                        Ctrl  0.07430946
# 8619.Ctrl                       Ctrl  0.39389371
# 8800.Ctrl                       Ctrl  0.06827974
# 142114.Thpok_LRF_min   Thpok_LRF_min  0.24603239
# 232614.Thpok_LRF_min   Thpok_LRF_min  0.31701027
# 432614.Thpok_LRF_min   Thpok_LRF_min -0.01729187
# 142114.Thpok_LRF_plus Thpok_LRF_plus -0.16831538
# 232614.Thpok_LRF_plus Thpok_LRF_plus -0.79414149
# 432614.Thpok_LRF_plus Thpok_LRF_plus -0.11977682

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x],las = 2, cex.lab = 0.8)
plotPCA(set2, col=colors[x], cex=0.8, xlim=c(-0.5,0.5))

differences <- matrix(data=c(1:3, 4:6), byrow=TRUE, nrow=3)
differences
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)

plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x],las = 2, cex.lab = 0.8)
plotPCA(set3, col=colors[x], cex=0.8, xlim=c(-0.5,0.5))

#Try RUVr
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

set4 <- RUVr(set, genes, k=1, res)
pData(set4)

plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[x],las = 2, cex.lab = 0.8)
plotPCA(set4, col=colors[x], cex=0.8, xlim=c(-0.5,0.5))


######Choose normalization method RUVr######################

design=model.matrix(~Mouse+Group)
colnames(design)
normExp=normCounts(set4)

v <- voom(normExp,design,plot=TRUE)
fit <- lmFit(v,design) 
fit <- eBayes(fit) 
topTable(fit,coef=ncol(design))
head(topTable(fit,coef=4,sort="none",n=Inf))
finalres=topTable(fit,coef=4,sort="none",n=Inf)


FC = 2^(-fit$coefficients[,4])
FC = ifelse(FC<1,-1/FC,FC)

finalres = cbind(v$E,FC,finalres)
finalres$logFC = -1*finalres$logFC
finalres[rownames(finalres)=="Zbtb7b",]  #Thpok
# 142114.Thpok_LRF_min 232614.Thpok_LRF_min 432614.Thpok_LRF_min 142114.Thpok_LRF_plus
# Zbtb7b                   58                   18                    7                 10763
# 232614.Thpok_LRF_plus 432614.Thpok_LRF_plus        FC    logFC AveExpr        t
# Zbtb7b                  7298                 12750 -473.0061 8.885715 2.11337 6.206964
# P.Value  adj.P.Val         B
# Zbtb7b 0.0003457501 0.04666803 -2.796885

finalres[rownames(finalres)=="Zbtb7a",]  #LRF
# 142114.Thpok_LRF_min 232614.Thpok_LRF_min 432614.Thpok_LRF_min 142114.Thpok_LRF_plus
# Zbtb7a                 8788                 8491                10651                 11984
# 232614.Thpok_LRF_plus 432614.Thpok_LRF_plus        FC     logFC  AveExpr       t
# Zbtb7a                 11779                 11631 -1.308252 0.3876408 6.656433 1.75685
# P.Value adj.P.Val         B
# Zbtb7a 0.1198189 0.4426547 -4.990207

write.table(finalres,file="YFPplus_minus_final_results.txt",row.names=TRUE,sep="\t",quote=FALSE)

######Choose normalization method Quantile######################

design=model.matrix(~Mouse+Group)
colnames(design)
normExp=filtered

v <- voom(normExp,design,plot=TRUE,normaliz="quantile")
fit <- lmFit(v,design) 
fit <- eBayes(fit) 
topTable(fit,coef=ncol(design))
head(topTable(fit,coef=4,sort="none",n=Inf))
finalres=topTable(fit,coef=4,sort="none",n=Inf)
hist(finalres$P.Value)

FC = 2^(-fit$coefficients[,4])
FC = ifelse(FC<1,-1/FC,FC)

finalres = cbind(v$E,FC,finalres)
finalres$logFC = -1*finalres$logFC
finalres[rownames(finalres)=="Zbtb7b",]  #Thpok
# 142114.Thpok_LRF_min 232614.Thpok_LRF_min 432614.Thpok_LRF_min 142114.Thpok_LRF_plus
# Zbtb7b            -1.674561            -1.874066            -3.000973              6.335185
# 232614.Thpok_LRF_plus 432614.Thpok_LRF_plus       FC 142114.Thpok_LRF_min 232614.Thpok_LRF_min
# Zbtb7b              6.782066              6.670319 -432.389            -1.674561            -1.874066
# 432614.Thpok_LRF_min 142114.Thpok_LRF_plus 232614.Thpok_LRF_plus 432614.Thpok_LRF_plus
# Zbtb7b            -3.000973              6.335185              6.782066              6.670319
# FC    logFC  AveExpr        t      P.Value adj.P.Val         B     ID
# Zbtb7b -432.389 8.756186 2.206328 5.172056 0.0007344314 0.1003251 -3.296954 Zbtb7b

finalres[rownames(finalres)=="Zbtb7a",]  #LRF
# 142114.Thpok_LRF_min 232614.Thpok_LRF_min 432614.Thpok_LRF_min 142114.Thpok_LRF_plus
# Zbtb7a             6.430526             6.308101             6.699206              7.044582
# 232614.Thpok_LRF_plus 432614.Thpok_LRF_plus        FC 142114.Thpok_LRF_min
# Zbtb7a              6.503431              6.900989 -1.263558             6.430526
# 232614.Thpok_LRF_min 432614.Thpok_LRF_min 142114.Thpok_LRF_plus 232614.Thpok_LRF_plus
# Zbtb7a             6.308101             6.699206              7.044582              6.503431
# 432614.Thpok_LRF_plus        FC     logFC  AveExpr        t   P.Value adj.P.Val         B
# Zbtb7a              6.900989 -1.263558 0.3374921 6.647806 1.341463 0.2149816 0.6761335 -5.392384
# ID
# Zbtb7a Zbtb7a

genename="Zbtb7b"
probeid=finalres[grepl(genename,perl=TRUE,finalres$ID)==TRUE,]$ID[1]
genedata=finalres[finalres$ID==probeid,1:6]
genedata = as.numeric(genedata)
stripchart((genedata)~Group,
           #stripchart((genedata)~pData(esetSub)$Dx,
           method="stack",
           col=c("red","blue"),
           offset=1.0, las=2,
           cex.axis=0.8,
           vertical = TRUE,
           ylab=genename)
dev.off()

head(finalres)
write.table(finalres,file="YFPplus_minus_final_results.txt",row.names=TRUE,sep="\t",quote=FALSE)

##############################For GO###################### 
# Write files for GO analysis:
##TR4:
finalres$ID = rownames(finalres)
FCresults = cbind(finalres$ID,finalres$FC,finalres$P.Value)
colnames(FCresults)=c("ID","FC","pval")

write.table(FCresults,file="ThpokLRF_results.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(finalres$ID,file="RNASeq.chip",quote=FALSE,row.names=FALSE,col.names=FALSE)

########## Do only dbl mutants ####################

subtargets=targets1[7:21,]

[,1]         [,2]         [,3]    
[1,] "wt YFP+"    "Foxp3-cre"  "8619"  
[2,] "wt YFP+"    "Foxp3-cre"  "124"   
[3,] "wt YFP+"    "Foxp3-cre"  "124"   
[4,] "wt YFP+"    "Foxp3-cre"  "8800"  
[5,] "wt YFP+"    "Foxp3-cre"  "8800"  
[6,] "wt YFP+"    "Foxp3-cre"  "8619"  
[7,] "dbl mutant" "dbl mutant" "132614"
[8,] "dbl mutant" "dbl mutant" "232614"
[9,] "dbl mutant" "dbl mutant" "232614"
[10,] "dbl mutant" "dbl mutant" "432614"
[11,] "dbl mutant" "dbl mutant" "332614"
[12,] "dbl mutant" "dbl mutant" "142114"
[13,] "dbl mutant" "dbl mutant" "142114"
[14,] "wt YFP-"    "dbl mutant" "142114"
[15,] "wt YFP-"    "dbl mutant" "142114"
[16,] "wt YFP-"    "dbl mutant" "232614"
[17,] "wt YFP-"    "dbl mutant" "232614"
[18,] "wt YFP-"    "dbl mutant" "432614"
[19,] "wt YFP-"    "dbl mutant" "432614"
[20,] "wt YFP-"    "dbl mutant" "332614"
[21,] "wt YFP-"    "dbl mutant" "332614"
[22,] "Thpok-d"    "Thpok-d"    "8799"  
[23,] "Thpok-d"    "Thpok-d"    "8799"  


#Voom: Run only on mice with dbl mutant
xSub=x[,7:21]
PhenoSub <- as.factor(subtargets$Phenotype)
GenoSub <- as.factor(subtargets$Genotype)
MouseSub <- as.factor(subtargets$Mouse)

design <- model.matrix(~MouseSub*PhenoSub)

##########Do Only the dbl mutants#########

design <- model.matrix(~0+PhenoSub)  
colnames(design) <- levels(PhenoSub) 
dupcor <- duplicateCorrelation(v1,design,block=MouseSub) 
dupcor$consensus.correlation

fit <- lmFit(v1,design,block=Mouse,
             correlation=dupcor$consensus.correlation) 

fit2 <- contrasts.fit(fit, contrasts) 
fit2 <- eBayes(fit2, trend=TRUE) 
summary(decideTests(fit2, method="global"))


### Do TMM normalization and check variance###
xSub <- calcNormFactors(xSub)
xSub$sample
v <- voom(xSub,design,plot=TRUE)
#v <- voom(counts,design,plot=TRUE)


###Do Quantile Normalization on Counts and get differential exp####
fc1=fc$counts[,-12]
fc.Counts.Sub=fc1[,7:21]
table(MouseSub,PhenoSub)

#MouseSub dbl mutant wt YFP- (technical replicates)
#132614          1       0
#142114          2       2
#232614          2       2
#332614          0       2
#432614          1       2
#511714          1       0

v <- voom(fc.Counts.Sub,design,plot=TRUE,normalize="quantile")
fit <- lmFit(v,design) 
fit <- eBayes(fit) 
topTable(fit,coef=ncol(design))

## Get dispersion stats
d <- estimateDisp(xSub, trend="none", robust=TRUE)
summary(d$prior.df)

## Draw the dispersion plot
jpeg(filename = "disp_plot_2.jpeg", width = 480, height = 480, units = "px", pointsize = 12,
     quality = 100)
plotBCV(d)
dev.off()

###### Filter data to where expression values are mostly more than 10 cpm

#isexpr <- rowSums(cpm(x) > 10) >= 2
log(20,2)
isexpr=(rowSums(v1$E > 4.3)) >=2
length(isexpr[isexpr==TRUE])
head(v1$E[isexpr==TRUE,])
v.expr <- v1$E[isexpr,]
#x_rpkm <- rpkm(x,x$genes$Length)

xsumm = cbind(xSub$genes,xSub$counts,x_rpkm)
write.table(xsumm,file="UCSC_notnormdata.txt",row.names=FALSE,sep="\t",quote=FALSE)
write.table(x$samples,file="total2.txt",row.names=FALSE,sep="\t",quote=FALSE)

#############PCA Plot########
#pca=prcomp(xSub$counts,scale.=T)

vdata=t(v1$E)
pca=prcomp(vdata)
pca1=princomp(v1$E)

### Get loadings #######

PC1_genes=order(abs(pca$rotation[,1]),decreasing=TRUE)[1:500]
PC1_genenames=rownames(pca$rotation)[PC1_genes]
PC1_genenames
head(2^fit2$coefficients[PC1_genes,])
head(fit2$p.value[PC1_genes,])

PC1_genes1=order(abs(pca1$scores[,1]),decreasing=TRUE)[1:500]
PC1_genenames1=rownames(pca1$scores)[PC1_genes1]
PC1_genenames1
head(2^fit2$coefficients[PC1_genes1,])
head(fit2$p.value[PC1_genes1,])

#PM=paste(PhenoSub,MouseSub,sep="_")
#PM=Group
PM=paste(Mouse,Group,sep="_")
library(rgl)
open3d() 
plot3d(pca$x[,1:3],col=as.numeric(as.factor(Group2)), 
       type="s",size=2)
group.v<-as.vector(Group2)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Qindiv.pdf","pdf")

###Plot PCA in 2D

tedf.class = as.character(Group2)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=tedf.class,
         ellipse=TRUE,circle=TRUE,var.axes=FALSE) 


######  Volcano Plot
for (i in 1:3) {
  cont=colnames(cont.matrix)[i]
  filename = paste("volcano",cont,"jpeg",sep=".")
  jpeg(file=filename)
  
  lod = fit2[["p.value"]][,i]
  lod = -log10(lod)
  mtstat = fit2[["t"]][,i]
  m = fit2[["coefficients"]][,i] 
  o1 = order(abs(m),decreasing=TRUE)[1:25]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:25]
  o=union(o1,o2)
  smoothScatter(m,lod,main=paste("Moderated t",i,sep="_"),
                xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")
  points(m[o2],lod[o2],pch=1,col="red")
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

##############  Run Heatmap ##############################
head(v1$E[PC1_genes,])
heatmap(v1$E[PC1_genes,], col=topo.colors(75), 
        ColSideColors=as.character(unclass(Group)),
         cexRow=0.5)

### Get Correlation Matrix ###########

library('corrplot')
PM=paste(Mouse,Group,sep="_")
PM
colnames(v1$E)=PM
#corv=cor(v1$E)
corv=cor(v.expr)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                           "cyan", "#007FFF", "blue","#00007F"))
corrplot(corv, order="hclust", col=col1(1000), 
         hclust.method="centroid", 
         is.corr=FALSE,addrect = 2,cl.lim=c(0.9,1.0))

#Use newly normalized Data
normExp=normCounts(set4)

#corv=cor(v1$E)
corv=cor(normExp)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                           "cyan", "#007FFF", "blue","#00007F"))
corrplot(corv, order="hclust", col=col1(1000), 
         hclust.method="centroid", 
         is.corr=FALSE,addrect = 2,cl.lim=c(0.9,1.0))

normExp.df=as.data.frame(normExp)
pairs(normExp.df)

#v1.df=as.data.frame(v1$E)
v1.expr.df=as.data.frame(v.expr)
pairs(v1.expr.df)
      
###Get MA plot###
logNE=log(normExp,2)
colnames(logNE)=colnames(normExp)

dim(logNE)
#[1] 11338     6

jpeg(filename="MA_plot.jpeg")
par(mfrow=c(2,3))
l = dim(logNE)[2]-1
for (i in 1:l) {
   j = i+1
   MA = cbind(logNE[,i],logNE[,j])
   colnames(MA) = c(colnames(logNE[,i]),colnames(logNE[,j]))
   plotMA(MA,main=colnames(MA))
   abline(h=0,col="red")
}
dev.off()

tiff(filename="MA_plot.tif")
M <- top.all$logFC #log fold-change
A <- top.all$AveExpr # ave. expr. level
mat <- cbind(A,M)           
plotMA(mat,main="MA plot")                   
abline(h=3,col="red")
dev.off()


## Get differentially expressed genes : Voom

y <- voom(x,design,plot=TRUE)
fit <- eBayes(lmFit(y,design))
ysumm = cbind(y$genes,y$E)
topTable(fit,coef=2)
write.table(ysumm,file="normdata2.txt",row.names=FALSE,sep="\t",quote=FALSE)

##Voom on 

jpeg('plotvar.jpg')
plot(y,xlim=c(-2.5,2.5))
dev.off()

jpeg('plotpca.jpg')
plotMDS(y,xlim=c(-2.5,2.5))
dev.off()

CT <- factor(celltype, levels = c("P6", "P5"))
#design <- model.matrix(~CT)


# colnames(design) <- paste("CT",levels(CT),sep="") 
# cont.matrix <- makeContrasts(P6.P5 = P6 - P5, levels = design)
# fit2 <- contrasts.fit(fit, cont.matrix) 
# fit2 <- eBayes(fit2)
# results <- decideTests(fit2) 
# write.fit(fit2,results=results,"diffexp.txt")

#A  Coef	t	p.value	F	F.p.value	P6.P5	Genes.GeneID	Genes.Chr	Genes.Start	Genes.End	Genes.Strand	Genes.Length
#4.4	0.067	0.52	0.61247	0.27	0.61247	0	27395	chr1;chr1;chr1;chr1;chr1	4773198;4777525;4782568;4783951;4785573	4776801;4777648;4782733;4784105;4785726	-;-;-;-;-	4203

### Write out final results:
finalres = topTable(fit,coef=2,sort="none",n=Inf)
sel.dif <- p.adjust(fit$F.p.value, method = "fdr") < 0.05

FC = 2^fit$coefficients
FC = ifelse(FC<1,-1/FC,FC) 
finalres = cbind(finalres,FC,y$E,sel.dif)
write.table(finalres,file="P6vsP5_final_results_per.txt",row.names=FALSE,sep="\t",quote=FALSE)


ratio2 = 2^(x_rpkm[,8]-x_rpkm[,7])
FC_expt4 = ifelse(ratio2<1,-1/ratio2,ratio2) 
FC_expt4 = cbind(x$genes,x_rpkm[,8],x_rpkm[,7],FC_expt4)
write.table(FC_expt4,file="P6vsP5_expt4_paired2r.txt",row.names=FALSE,sep="\t",quote=FALSE)


v <- voom(fc$counts, design) 
fit <- eBayes(lmFit(v, design)) 
topTable(fit, coef=2)


######Try alignment#####

ref <- "chr_all3.fa"
#buildindex(basename="reference_index",reference=ref)
buildindex(basename="reference_index",reference=ref,gappedIndex=TRUE,indexSplit=TRUE,memory=8000,
           TH_subread=24,colorspace=FALSE)

read1 = read.table("read1.txt")
read1=(as.vector(as.matrix(read1)))
read2 = read.table("read2.txt")
read2=(as.vector(as.matrix(read2)))

#For subread (non-RNASeq)
align(index="reference_index",readfile1=read1,readfile2=read2)
#For RNA-Seq data:
subjunc(index="reference_index",readfile1=read1,readfile2=read2)

##########################################

save.image()
rm(list=ls())
