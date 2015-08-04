setwd("/Users/maggiec/Documents/Deblina")
setwd("/Volumes/CCBR/CCRIFXCCR/ccbr482/RNASeq/bam/")
load("Gencode_fc.RData")
load(".RData")
targets <- readTargets()

fc=fc_paired_r
ls()
library(edgeR)
options(digits=2)

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
celltype <- factor(targets$condition)
design <- model.matrix(~celltype)

targets
v1 <- voom(as.matrix(fc1filt),design,lib.size =x$samples[,2],
           normalize="quantile")

dim(v1)
#[1] 11947    18
#[1] 15111    18  (including new lncRNA)

v <- voom(fc$counts, design) 
fit <- eBayes(lmFit(v, design)) 
topTable(fit, coef=2)

v2 <- voom(as.matrix(fc1filt),design,lib.size =x$samples[,2],
           normalize="none")

#####Look after quantile normalization:

head(v1$weights)
df.n <- melt(as.data.frame(v1$E))
df.n <- melt(as.data.frame(v2$E))  #look at non-normalized

jpeg(filename = "disp_plot_nonnorm.jpeg", width = 480, height = 480, 
     units = "px", pointsize = 12,quality = 100)
jpeg(filename = "disp_plot_norm.jpeg", width = 480, height = 480, 
     units = "px", pointsize = 12,quality = 100)
ggplot(df.n) + geom_density(aes(x = value,
                                colour = variable)) + labs(x = NULL) +
  theme(legend.position='top')
dev.off()


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
#v <- voom(fc$counts,design,plot=TRUE)
fit <- lmFit(v,design) 
fit <- eBayes(fit) 
topTable(fit,coef=ncol(design))
head(topTable(fit,coef=2,sort="none",n=Inf))
finalres=topTable(fit,coef=2,sort="none",n=Inf)


FC = 2^(fit$coefficients[,2])
FC = ifelse(FC<1,-1/FC,FC)

finalres = cbind(v$E,FC,finalres)

finalres[rownames(finalres)=="RARB",] 


write.table(finalres,file="Gencode_FCpval.txt",
        row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
