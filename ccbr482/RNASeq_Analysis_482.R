setwd("/Users/maggiec/Documents/Deblina")
setwd("/Volumes/CCBR/CCRIFXCCR/ccbr482/RNASeq/bam/")
setwd("/Users/maggiec/GitHub/Maggie/ccbr482")
load("Gencode_fc.RData")
load(".RData")
targets <- readTargets()

biocLite("biomaRt")
library(biomaRt)
###Done 

gtf="/data/maggiec/RNASeq/Genomes/hg19/gencode.v22.annotation.gtf"
library(Rsubread)
library(limma)
targets <- readTargets()

fc_paired_r <- featureCounts(files=targets$bam,isGTFAnnotationFile=TRUE,nthreads=32,
            annot.ext=gtf,GTF.attrType="gene_name",strandSpecific=2,isPairedEnd=TRUE)

fc=fc_paired_r
ls()
library(edgeR)
options(digits=3)

x <- DGEList(counts=fc$counts, genes=fc$annotation)

fc1=mat=fc$counts
tfc1=t(fc1)
filter <- apply(fc1, 1, function(x) length(x[x>5])>=1)
fc1filt <- fc1[filter,]
genes <- rownames(fc1filt)

####Histogram: Look at distribution of data before quantile:
library(reshape)
df.m <- melt(as.data.frame(fc1))
jpeg(filename = "disp_plot_raw.jpeg", width = 480, height = 480, 
     units = "px", pointsize = 12,quality = 100)

ggplot(df.m) + geom_density(aes(x = value,
                                colour = variable)) + labs(x = NULL) +
  theme(legend.position='top') + scale_x_log10()

dev.off()

###Do Normalization: cpm is further normalized by quantile
##v1 = scaled by lib.size then quantile norm
x$counts[x$genes$GeneID=="MYCN",]
x$counts[x$genes$GeneID=="RARB",]

x <- DGEList(counts=fc$counts, genes=fc$annotation)
isexpr <- rowSums(cpm(x)>1) >= 2
hasannot <- rowSums(is.na(x$genes))==0
x <- x[isexpr & hasannot,,keep.lib.sizes=FALSE]
dim(x)
[1] 14109     6

x <- calcNormFactors(x)
celltype <- factor(targets$condition)
design <- model.matrix(~celltype)
design

targets
v <- voom(as.matrix(x),design,plot=TRUE,normalize="quantile")

plotMDS(v,top=50,labels=substring(celltype,1,1), 
      col=ifelse(celltype=="C","blue","red"),gene.selection="common")

v2 <- voom(as.matrix(x),design,lib.size =x$samples[,2],
           normalize="none")

fit <- lmFit(v,design) 
fit <- eBayes(fit) 
options(digits=3) 
topTable(fit,coef=2,n=16,sort="p")
volcanoplot(fit,coef=2)  

#####Write out Results ############

finalres=topTable(fit,coef=2,sort="none",n=Inf)

FC = 2^(fit$coefficients[,2])
FC = ifelse(FC<1,-1/FC,FC)

finalres = cbind(v$E,FC,finalres)

finalres[rownames(finalres)=="RARB",] 
finalres[rownames(finalres)=="MYCN",] 


datadir="/Volumes/maggiec/Projects/ccbr482/Results"
setwd(datadir)
write.table(finalres,file="Gencode_FCpval_new.txt",
      row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

#####Look after quantile normalization:

head(v1$weights)
df.n <- melt(as.data.frame(v$E))
df.n <- melt(as.data.frame(v2$E))  #look at non-normalized

jpeg(filename = "disp_plot_nonnorm.jpeg", width = 480, height = 480, 
     units = "px", pointsize = 12,quality = 100)
jpeg(filename = "disp_plot_norm.jpeg", width = 480, height = 480, 
     units = "px", pointsize = 12,quality = 100)
ggplot(df.n) + geom_density(aes(x = value,
                                colour = variable)) + labs(x = NULL) +
  theme(legend.position='top')
dev.off()

datadir="/Volumes/maggiec/Projects/ccbr482/Results"
setwd(datadir)
write.table(v1$E,file="Gencode_voom_Q_Ctrl.txt",
            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(v1$E,file="Gencode_voom_none.txt",
            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

########################################################
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

tfc1=t(fc1)


#####Tried scaling first by total library counts - offset doesn't work as well)
lmfitcoeff=vector()
scale=matrix(nrow=dim(v1$E)[1],ncol=dim(v1$E)[2])
scaleoff=matrix(nrow=dim(v1$E)[1],ncol=dim(v1$E)[2])
for (i in 1:12){
  x=log2(fc1filt[,i]/v1$targets[i,1])
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

scaled_dat=(fc1filt)*(scaleoff)
scaleval=(2^v1$E)*100  # Real data (aiming for)
actscale=(fc1filt/scaleval)

scale_fit=(fc1filt/actscale) ###Need to use this formula to calculate normalized counts
sf = v1$E/log2((fc1filt/colSums(fc1filt))*1000000)
scale_fit2=(fc1filt/sf)

scaled_dat[rownames(scaled_dat)=="Zbtb7b"]
scaleval[rownames(scaleval)=="Zbtb7b"]
scale_fit[rownames(scale_fit)=="Zbtb7b"]
scale_fit2[rownames(scale_fit2)=="Zbtb7b"]
fc1filt[rownames(sumfcfilt)=="Zbtb7b"]

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

write.table(actscale,file="scaling_factors.txt",
            row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)

Group = as.factor(targets$condition)
design=model.matrix(~0+Group)
colnames(design)
v=v1
v=v2

fit <- lmFit(v,design) 
fit <- eBayes(fit) 


plot(fit,coef=2, main="RA vs C", legend="bottomright") 

abline(h=0,col="darkgrey") 


gtf_annot=read.delim("data/gencode.v19.annotation.gtf")

gtf_annot=read.delim("data/gencode.trunc", skip=5,header=FALSE)
gene_annot=as.character(gtf_annot$V9)
gene_list=strsplit(gene_annot, split=";")
 
gene_name=strsplit(gene_list[[1]][5],split=" ")
gene_type=strsplit(gene_list[[1]][3],split=" ")

gene_info=paste(gene_name[[1]][3],gene_type[[1]][3])
