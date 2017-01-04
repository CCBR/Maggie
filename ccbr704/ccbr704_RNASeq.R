setwd("/Users/maggiec/GitHub/maggie/ccbr704/data/")
load("CCBR704.RData")
load(".RData")

DP1=read.table("CPM_TMM_counts.txt",header = TRUE,sep = "\t")
DP2=read.table("CPM_TMM_counts2.txt",header = TRUE,sep = "\t")

colnames(DP2)
head(DP2[,c(1,6,10,13)])
DP2.mat=as.matrix(DP2[,c(6,10,13)])
rownames(DP2.mat)=sapply(strsplit(as.character(DP2$X),split="\\|"),function(x) x[2])
head(DP2.mat)

TEC_counts=as.data.frame(fc$counts[,match(targets1$bam,colnames(fc$counts))])
dim(TEC_counts)
[1] 43280    18

dim(DP2.df)
[1] 15491     4

TEC_DP_counts=merge(DP2.df,TEC_counts,by.x="gene",by.y="gene")

head(TEC_DP_counts)
dim(TEC_DP_counts)
[1] 15491    21

minread=5
minread=(minread/max(colSums(TEC_DP_counts[,2:21])))*1e6  
minread
minsamples=5
filter <- apply(cpm(TEC_DP_counts[,2:21]), 1, 
                function(x) length(x[x>minread])>=minsamples)
TEC_DP_counts.filt=cpm(TEC_DP_counts[filter,2:21])
rownames(TEC_DP_counts.filt)=TEC_DP_counts$gene[filter]

DPcounts=rowMeans(TEC_DP_counts.filt[,1:3])
TECcounts=rowMeans(TEC_DP_counts.filt[,4:20])
filt2=(DPcounts/TECcounts)>100
length(rownames(TEC_DP_counts.filt)[filt2])
filtgenes=rownames(TEC_DP_counts.filt)[filt2]

targets
#bam  Time  Type Rep
1     Sample_10_E15_5_ETP.bam E15.5   ETP   1
2    Sample_11_E17_5_cTEC.bam E17.5  cTEC   1
3  Sample_12_E17_5_mTEC_1.bam E17.5  mTEC   1
4  Sample_13_E17_5_cTEC_2.bam E17.5  cTEC   2
5  Sample_14_E17_5_mTEC_2.bam E17.5  mTEC   2
6     Sample_15_NB_cTEC_1.bam    NB  cTEC   1
7     Sample_16_NB_mTEC_1.bam    NB  mTEC   1
8     Sample_17_NB_cTEC_2.bam    NB  cTEC   2
9     Sample_18_NB_mTEC_2.bam    NB  mTEC   2
10       Sample_19_NB_ETP.bam    NB   ETP   1
11   Sample_1_E13_5_TEC_1.bam E13.5   TEC   1
12      Sample_20_4W_cTEC.bam    4W  cTEC   1  #potential outlier
13      Sample_21_4W_mTEC.bam    4W  mTEC   1
14      Sample_22_5W_cTEC.bam    5W  cTEC   1
15      Sample_23_5W_mTEC.bam    5W  mTEC   1
16   Sample_24_E13_5_CD45.bam E13.5  CD45   1
17   Sample_2_E13_5_TEC_2.bam E13.5   TEC   2
18     Sample_3_E13_5_ETP.bam E13.5   ETP   1
19  Sample_4_E15_5_cTEC_1.bam E15.5  cTEC   1
20  Sample_5_E15_5_mTEC_1.bam E15.5  mTEC   1
21 Sample_6_E15_5_inter_1.bam E15.5 inter   1
22  Sample_7_E15_5_cTEC_2.bam E15.5  cTEC   2
23  Sample_8_E15_5_mTEC_2.bam E15.5  mTEC   2
24 Sample_9_E15_5_inter_2.bam E15.5 inter   2

fc1=mat=fc$counts
tfc1=t(fc1)
filter <- apply(fc1, 1, function(x) length(x[x>5])>=2)
fc1filt <- fc1[filter,] 
#fc1filt2=fc1[rownames(fc1) %in% selectgenes$V1,]

genes <- rownames(fc1filt)
#targets1=targets[n,]
dim(fc1filt)
```

####*Mapping Rate:*

```{r, fig.width=10, fig.height=5, warning=FALSE,echo=FALSE}
#options(scipen=9)
#map=as.numeric(targets$mapped)
#names(map)=sapply(strsplit(targets$bam,split="\\."),function(x) x[1])
#barplot(map, col=as.numeric(as.factor(targets$Phenotype)),las = 2,cex.axis = 0.7)
```


######Bring in new Data:#####
dat=read.table("RawCountFile_genes_filtered.txt",header = TRUE,sep = "\t")
colnames(dat)
[1] "X"                             "X10W_cTEC_2.star.count.txt"   
[3] "X10W_mTEC_1.star.count.txt"    "X10W_mTEC_2.star.count.txt"   
[5] "X10W_mTEC_3.star.count.txt"    "X11M_mTEC_S24.star.count.txt" 
[7] "X4W_cTEC_3.star.count.txt"     "X4W_cTEC.star.count.txt"      
[9] "X4W_mTEC_1.star.count.txt"     "X4W_mTEC_2.star.count.txt"    
[11] "X4W_mTEC_3.star.count.txt"     "X4W_mTEC.star.count.txt"      
[13] "X5W_cTEC.star.count.txt"       "X5W_mTEC.star.count.txt"      
[15] "X8W_cTEC_S12.star.count.txt"   "X8W_mTEC_S11.star.count.txt"  
[17] "X9M_cTEC_S31.star.count.txt"   "X9M_mTEC_S30.star.count.txt"  
[19] "E13_5_CD45.star.count.txt"     "E13_5_ETP.star.count.txt"     
[21] "E13_5_TEC_1.star.count.txt"    "E13_5_TEC_2.star.count.txt"   
[23] "E15_5_cTEC_1.star.count.txt"   "E15_5_cTEC_2.star.count.txt"  
[25] "E15_5_ETP.star.count.txt"      "E15_5_inter_1.star.count.txt" 
[27] "E15_5_inter_2.star.count.txt"  "E15_5_mTEC_1.star.count.txt"  
[29] "E15_5_mTEC_2.star.count.txt"   "E16p5_cTEC_S21.star.count.txt"
[31] "E16p5_mTEC_S20.star.count.txt" "E17_5_cTEC_1.star.count.txt"  
[33] "E17_5_cTEC_2.star.count.txt"   "E17_5_mTEC_1.star.count.txt"  
[35] "E17_5_mTEC_2.star.count.txt"   "Myc_1_cTEC.star.count.txt"    
[37] "Myc_1_mTEC.star.count.txt"     "Myc_2_cTEC.star.count.txt"    
[39] "Myc_2_mTEC.star.count.txt"     "Myc_3_cTEC.star.count.txt"    
[41] "Myc_3_mTEC.star.count.txt"     "Myc_3_mTEC_rep.star.count.txt"
[43] "NB1_mTEC_1.star.count.txt"     "NB1_mTEC_2.star.count.txt"    
[45] "NB_cTEC_1.star.count.txt"      "NB_cTEC_2.star.count.txt"     
[47] "NB_mTEC_1.star.count.txt"   

colnames(dat)=sapply(strsplit(as.character(colnames(dat)),split="\\."),function(x) x[1])
[1] "X"             
[2] "X10W_cTEC_2"   
[3] "X10W_mTEC_1"   
[4] "X10W_mTEC_2"   
[5] "X10W_mTEC_3"   
[6] "X11M_mTEC_S24" 
[7] "X4W_cTEC_3*"    
[8] "X4W_cTEC*"      
[9] "X4W_mTEC_1*"    
[10] "X4W_mTEC_2*"    
[11] "X4W_mTEC_3*"    
[12] "X4W_mTEC*"      
[13] "X5W_cTEC*"      
[14] "X5W_mTEC*"      
[15] "X8W_cTEC_S12"  
[16] "X8W_mTEC_S11"  
[17] "X9M_cTEC_S31"  
[18] "X9M_mTEC_S30"  
[19] "E13_5_CD45*"    
[20] "E13_5_ETP*"     
[21] "E13_5_TEC_1*"   
[22] "E13_5_TEC_2*"   
[23] "E15_5_cTEC_1*"  
[24] "E15_5_cTEC_2*"  
[25] "E15_5_ETP*"     
[26] "E15_5_inter_1*" 
[27] "E15_5_inter_2*" 
[28] "E15_5_mTEC_1*"  
[29] "E15_5_mTEC_2*"  
[30] "E16p5_cTEC_S21"
[31] "E16p5_mTEC_S20"
[32] "E17_5_cTEC_1*"  
[33] "E17_5_cTEC_2*"  
[34] "E17_5_mTEC_1*"  
[35] "E17_5_mTEC_2*"  
[36] "Myc_1_cTEC"    
[37] "Myc_1_mTEC"    
[38] "Myc_2_cTEC"    
[39] "Myc_2_mTEC"    
[40] "Myc_3_cTEC"    
[41] "Myc_3_mTEC"    
[42] "Myc_3_mTEC_rep"
[43] "NB1_mTEC_1*"    
[44] "NB1_mTEC_2*"    
[45] "NB_cTEC_1"     
[46] "NB_cTEC_2"     
[47] "NB_mTEC_1" 

rownames(dat)=sapply(strsplit(as.character(dat$X),split="\\|"),function(x) x[1])

Samp.info=data.frame(sample=colnames(dat)[2:dim(dat)[2]])
Samp.info$cell=sapply(strsplit(as.character(Samp.info$sample),split="\\_"),
                      function(x) x[2])

##Created sampinfo in excel:
Samp.info=read.table("sampinfo.txt",header = TRUE,row.names = 1)


####*QC Check: Look at raw signal distribution and median expression levels*

```{r fig.width=10, fig.height=5, echo=FALSE, warning=FALSE,message=FALSE}

library(ggplot2)
library(reshape)
library(plotly)

fc1filt=as.data.frame(fc1filt)
#fc1filtnames=sapply(strsplit(colnames(fc1filt),split="\\."),function(x) x[1])
#fc1filtnames=sapply(strsplit(fc1filtnames,split="_"),function(x) paste(x[2],x[3],sep = "."))
fc1filtnames=paste(targets$Time,targets$Type,targets$Rep,sep="_")
fc1filtnames
#colnames(fc1filt)=paste(targets$Sample,targets$Batch)
colnames(fc1filt)=fc1filtnames

index=isexpr & hasannot & isexprzero
fc2filt=fc1filt[index,]
df.m <- melt(as.data.frame(fc1filt))

###Plot raw data #####
dat.batch=dat[,colnames(dat) ]
df.m <- melt(as.data.frame(dat[,-1]))


ggplot(df.m, aes(x=value, group=variable)) + 
  geom_density(aes(colour = variable, linetype = variable),
               size=1, alpha=.60, kernel='epanechnikov') +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position='right') + scale_x_log10() + ggtitle("Raw Counts") +
  scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),18))

par(mar=c(10,7,1,1))
boxplot(log(value)~variable,las=2,data=df.m,main="Raw Signal", 
        ylab="Counts",col=as.numeric(as.factor(targets$Type)))
```

###This is for all batches:
Samp.info

cell  Time Batch Rep
X10W_cTEC_2     cTEC   10W     2   2
X10W_mTEC_1     mTEC   10W     2   1
X10W_mTEC_2     mTEC   10W     2   2
X10W_mTEC_3     mTEC   10W     2   3
X11M_mTEC_S24   mTEC   11M     2   1
X4W_cTEC_3      cTEC    4W     2   1
X4W_cTEC        cTEC    4W     1   1
X4W_mTEC_1      mTEC    4W     2   1
X4W_mTEC_2      mTEC    4W     2   2
X4W_mTEC_3      mTEC    4W     2   3
X4W_mTEC        mTEC    4W     1   1
X5W_cTEC        cTEC    4W     1   2
X5W_mTEC        mTEC    4W     1   2
X8W_cTEC_S12    cTEC    8W     2   1
X8W_mTEC_S11    mTEC    8W     2   1
X9M_cTEC_S31    cTEC    9M     2   1
X9M_mTEC_S30    mTEC    9M     2   1
E13_5_CD45      CD45 E13.5     1   1
E13_5_ETP        ETP E13.5     1   1
E13_5_TEC_1      TEC E13.5     1   1
E13_5_TEC_2      TEC E13.5     1   2
E15_5_cTEC_1    cTEC E15.5     1   1
E15_5_cTEC_2    cTEC E15.5     1   2
E15_5_ETP        ETP E15.5     1   1
E15_5_inter_1  inter E15.5     1   1
E15_5_inter_2  inter E15.5     1   2
E15_5_mTEC_1    mTEC E15.5     1   1
E15_5_mTEC_2    mTEC E15.5     1   2
E16p5_cTEC_S21  cTEC E16.5     2   1
E16p5_mTEC_S20  mTEC E16.5     2   1
E17_5_cTEC_1    cTEC E17.5     1   1
E17_5_cTEC_2    cTEC E17.5     1   2
E17_5_mTEC_1    mTEC E17.5     1   1
E17_5_mTEC_2    mTEC E17.5     1   2
Myc_1_cTEC      cTEC   Myc     2   1
Myc_1_mTEC      mTEC   Myc     2   1
Myc_2_cTEC      cTEC   Myc     2   2
Myc_2_mTEC      mTEC   Myc     2   2
Myc_3_cTEC      cTEC   Myc     2   3
Myc_3_mTEC      mTEC   Myc     2   3
Myc_3_mTEC_rep  mTEC   Myc     2   3
NB1_mTEC_1      mTEC    NB     1   1
NB1_mTEC_2      mTEC    NB     1   2
NB_cTEC_1       cTEC    NB     2   1
NB_cTEC_2       cTEC    NB     2   2
NB_mTEC_1       mTEC    NB     2   1

c <- factor(Samp.info$cell) 
t <- factor(Samp.info$Time) 
b <- factor(Samp.info$Batch)
design=model.matrix(~b+c+t)
design

v <- voom(dat[,-1],design=design,normalize="quantile",plot=TRUE)
df.m=melt(as.data.frame(v$E))

####*Data is normalized by TMM: filtered number of genes*

```{r,fig.width=10, fig.height=5, echo=FALSE, warning=FALSE, message=FALSE}

selectgenes=c("Id3","Nfe2","Klf2","Nrarp","Rxrg","Tox2","Nfil3","Rora","Id2",
              "Zbtb16","Rorc","Lef1","Mycn","Gata2","Tcf7","Tox")
#selectgenes=read.table("data/select_genes.txt")
#selectgenes=c(as.character(selectgenes$V1),"Tcf7")
Genecount=fc$counts[rownames(fc$counts) %in% selectgenes,]

#Select Samples for further analysis:
n=c(-9,-15,-16)
targets1=targets[n,]

#x <- DGEList(counts=fc$counts[,n], genes=fc$annotation)
x <- DGEList(counts=fc$counts, genes=fc$annotation)
isexpr <- rowSums(cpm(x)>1.0) >= 2
isexpr[names(isexpr) %in% selectgenes] = TRUE
isexprzero <- rowSums(cpm(x))>0
hasannot <- rowSums(is.na(x$genes))==0
x <- x[isexpr & hasannot & isexprzero,,keep.lib.sizes=FALSE]

#Number of filtered genes 
dim(x)
[1] 20095    24

#x=rpkm(x)
x <- calcNormFactors(x)


par(mfrow=c(2,4))
for (i in 1:dim(x)[2]){
  plotMD(cpm(x, log=TRUE), column=i)
  abline(h=0, col="red", lty=2, lwd=2)
}

```

####*Run Voom*
```{r, echo=FALSE,warning=FALSE,message=FALSE}
#Do analysis for entire group

targets$Time2=targets$Time

#Change Time=5W to Time=4W
targets$Rep[targets$Time2=="5W"]="2"
targets$Time2[targets$Time2=="5W"]="4W"
targets$Typetime = paste(targets$Type,targets$Time2,sep="_")
targets$sample=paste(targets$Typetime,targets$Rep,sep="_")


targets$Time3[targets$Time2=="E13.5"]="0"
targets$Time3[targets$Time2=="E15.5"]="1"
targets$Time3[targets$Time2=="E17.5"]="2"
targets$Time3[targets$Time2=="NB"]="3"
targets$Time3[targets$Time2=="4W"]="4"

targets

bam  Time  Type Rep    Typetime Time2        sample Time3
1     Sample_10_E15_5_ETP.bam E15.5   ETP   1   ETP_E15.5 E15.5   ETP_E15.5_1     1
2    Sample_11_E17_5_cTEC.bam E17.5  cTEC   1  cTEC_E17.5 E17.5  cTEC_E17.5_1     2
3  Sample_12_E17_5_mTEC_1.bam E17.5  mTEC   1  mTEC_E17.5 E17.5  mTEC_E17.5_1     2
4  Sample_13_E17_5_cTEC_2.bam E17.5  cTEC   2  cTEC_E17.5 E17.5  cTEC_E17.5_2     2
5  Sample_14_E17_5_mTEC_2.bam E17.5  mTEC   2  mTEC_E17.5 E17.5  mTEC_E17.5_2     2
6     Sample_15_NB_cTEC_1.bam    NB  cTEC   1     cTEC_NB    NB     cTEC_NB_1     3
7     Sample_16_NB_mTEC_1.bam    NB  mTEC   1     mTEC_NB    NB     mTEC_NB_1     3
8     Sample_17_NB_cTEC_2.bam    NB  cTEC   2     cTEC_NB    NB     cTEC_NB_2     3
9     Sample_18_NB_mTEC_2.bam    NB  mTEC   2     mTEC_NB    NB     mTEC_NB_2     3
10       Sample_19_NB_ETP.bam    NB   ETP   1      ETP_NB    NB      ETP_NB_1     3
11   Sample_1_E13_5_TEC_1.bam E13.5   TEC   1   TEC_E13.5 E13.5   TEC_E13.5_1     0
12      Sample_20_4W_cTEC.bam    4W  cTEC   1     cTEC_4W    4W     cTEC_4W_1     4  #outlier
13      Sample_21_4W_mTEC.bam    4W  mTEC   1     mTEC_4W    4W     mTEC_4W_1     4
14      Sample_22_5W_cTEC.bam    5W  cTEC   2     cTEC_4W    4W     cTEC_4W_2     4
15      Sample_23_5W_mTEC.bam    5W  mTEC   2     mTEC_4W    4W     mTEC_4W_2     4
16   Sample_24_E13_5_CD45.bam E13.5  CD45   1  CD45_E13.5 E13.5  CD45_E13.5_1     0
17   Sample_2_E13_5_TEC_2.bam E13.5   TEC   2   TEC_E13.5 E13.5   TEC_E13.5_2     0
18     Sample_3_E13_5_ETP.bam E13.5   ETP   1   ETP_E13.5 E13.5   ETP_E13.5_1     0
19  Sample_4_E15_5_cTEC_1.bam E15.5  cTEC   1  cTEC_E15.5 E15.5  cTEC_E15.5_1     1
20  Sample_5_E15_5_mTEC_1.bam E15.5  mTEC   1  mTEC_E15.5 E15.5  mTEC_E15.5_1     1
21 Sample_6_E15_5_inter_1.bam E15.5 inter   1 inter_E15.5 E15.5 inter_E15.5_1     1
22  Sample_7_E15_5_cTEC_2.bam E15.5  cTEC   2  cTEC_E15.5 E15.5  cTEC_E15.5_2     1
23  Sample_8_E15_5_mTEC_2.bam E15.5  mTEC   2  mTEC_E15.5 E15.5  mTEC_E15.5_2     1
24 Sample_9_E15_5_inter_2.bam E15.5 inter   2 inter_E15.5 E15.5 inter_E15.5_2     1

table(targets$Type,targets$Time2)
#      4W   E13.5 E15.5 E17.5 NB
CD45   0     1     0     0  0
cTEC   2     0     2     2  2
ETP    0     1     1     0  1
inter  0     0     2     0  0
mTEC   2     0     2     2  2
TEC    0     2     0     0  0

table(targets$Type,targets$Time3)
      0 1 2 3 4
CD45  1 0 0 0 0
cTEC  0 2 2 2 2
ETP   1 1 0 1 0
inter 0 2 0 0 0
mTEC  0 2 2 2 2
TEC   2 0 0 0 0

#celltype <- factor(targets$Type, levels=c("ETP","cTEC","mTEC","TEC","inter"))
#time <- factor(targets$Time, levels=c("E13.5","E15.5","E17.5","NB","4W","5W"))
#design <- model.matrix(~0+celltype+time)

####Original - looking at differences within the different time points####
lev <- c("cTEC_E15.5","cTEC_E17.5","cTEC_NB","cTEC_4W",
         "mTEC_E15.5","mTEC_E17.5","mTEC_NB","mTEC_4W") 
f <- factor(targets$Typetime, levels=lev) 

design <- model.matrix(~0+f) 
colnames(design)=lev

design

#    cTEC_E15.5 cTEC_E17.5 cTEC_NB cTEC_4W mTEC_E15.5 mTEC_E17.5 mTEC_NB mTEC_4W
2           0          1       0       0          0          0       0       0
3           0          0       0       0          0          1       0       0
4           0          1       0       0          0          0       0       0
5           0          0       0       0          0          1       0       0
6           0          0       1       0          0          0       0       0
7           0          0       0       0          0          0       1       0
8           0          0       1       0          0          0       0       0
9           0          0       0       0          0          0       1       0
12          0          0       0       1          0          0       0       0
13          0          0       0       0          0          0       0       1
14          0          0       0       1          0          0       0       0
15          0          0       0       0          0          0       0       1
19          1          0       0       0          0          0       0       0
20          0          0       0       0          1          0       0       0
22          1          0       0       0          0          0       0       0
23          0          0       0       0          1          0       0       0

dev.off()

n=c(-1,-10,-16,-18,-21,-24)

n=c(-1,-10,-16,-18,-21,-24,-12)  #Remove outlier 12

targets1=targets[n,]
targets1

> targets1
               bam  Time Type Rep   Typetime Time2       sample Time3
2    Sample_11_E17_5_cTEC.bam E17.5 cTEC   1 cTEC_E17.5 E17.5 cTEC_E17.5_1     2
3  Sample_12_E17_5_mTEC_1.bam E17.5 mTEC   1 mTEC_E17.5 E17.5 mTEC_E17.5_1     2
4  Sample_13_E17_5_cTEC_2.bam E17.5 cTEC   2 cTEC_E17.5 E17.5 cTEC_E17.5_2     2
5  Sample_14_E17_5_mTEC_2.bam E17.5 mTEC   2 mTEC_E17.5 E17.5 mTEC_E17.5_2     2
6     Sample_15_NB_cTEC_1.bam    NB cTEC   1    cTEC_NB    NB    cTEC_NB_1     3
7     Sample_16_NB_mTEC_1.bam    NB mTEC   1    mTEC_NB    NB    mTEC_NB_1     3
8     Sample_17_NB_cTEC_2.bam    NB cTEC   2    cTEC_NB    NB    cTEC_NB_2     3
9     Sample_18_NB_mTEC_2.bam    NB mTEC   2    mTEC_NB    NB    mTEC_NB_2     3
11   Sample_1_E13_5_TEC_1.bam E13.5  TEC   1  TEC_E13.5 E13.5  TEC_E13.5_1     0
12      Sample_20_4W_cTEC.bam    4W cTEC   1    cTEC_4W    4W    cTEC_4W_1     4   #outlier
13      Sample_21_4W_mTEC.bam    4W mTEC   1    mTEC_4W    4W    mTEC_4W_1     4
14      Sample_22_5W_cTEC.bam    5W cTEC   2    cTEC_4W    4W    cTEC_4W_2     4
15      Sample_23_5W_mTEC.bam    5W mTEC   2    mTEC_4W    4W    mTEC_4W_2     4
17   Sample_2_E13_5_TEC_2.bam E13.5  TEC   2  TEC_E13.5 E13.5  TEC_E13.5_2     0
19  Sample_4_E15_5_cTEC_1.bam E15.5 cTEC   1 cTEC_E15.5 E15.5 cTEC_E15.5_1     1
20  Sample_5_E15_5_mTEC_1.bam E15.5 mTEC   1 mTEC_E15.5 E15.5 mTEC_E15.5_1     1
22  Sample_7_E15_5_cTEC_2.bam E15.5 cTEC   2 cTEC_E15.5 E15.5 cTEC_E15.5_2     1
23  Sample_8_E15_5_mTEC_2.bam E15.5 mTEC   2 mTEC_E15.5 E15.5 mTEC_E15.5_2     1

table(targets1$Type,targets1$Time2)

      4W E13.5 E15.5 E17.5 NB
cTEC  2     0     2     2  2
mTEC  2     0     2     2  2
TEC   0     2     0     0  0

celltype1=targets1$Type
x1=x[,n]
x1$samples

#Another filtration method:
library(HTSFilter)
x_filt=HTSFilter(x1$counts, targets1$Type, s.len=25, plot=TRUE)
hist(log(x_filt$filteredData+1), col="grey", breaks=25, main="",
     xlab="Log(counts+1)")
remdata=rownames(x_filt$removedData)
dim(x_filt$filteredData)

colnames(x1$counts)
dim(x1$counts)
isexprzero <- rowSums(cpm(x1))>0
hasannot <- rowSums(is.na(x1$genes))==0
x1 <- x1[hasannot & isexprzero,,keep.lib.sizes=FALSE]
x_filt=HTSFilter(x1$counts, targets1$Typetime,s.len=25, plot=TRUE)
class(x_filt)
hist(log(x_filt$filteredData+1), col="grey", breaks=25, main="",
     xlab="Log(counts+1)")
remdata=rownames(x_filt$removedData)


lev <- c("TEC_E13.5",
         "cTEC_E15.5","cTEC_E17.5","cTEC_NB","cTEC_4W",
         "mTEC_E15.5","mTEC_E17.5","mTEC_NB","mTEC_4W") 
f <- factor(targets1$Typetime, levels=lev) 

design <- model.matrix(~0+f) 
colnames(design)=lev
dim(design)  #17

design

       TEC_13.5 cTEC_E15.5 cTEC_E17.5 cTEC_NB cTEC_4W mTEC_E15.5 mTEC_E17.5 mTEC_NB mTEC_4W
1         0          0          1       0       0          0          0       0       0
2         0          0          0       0       0          0          1       0       0
3         0          0          1       0       0          0          0       0       0
4         0          0          0       0       0          0          1       0       0
5         0          0          0       1       0          0          0       0       0
6         0          0          0       0       0          0          0       1       0
7         0          0          0       1       0          0          0       0       0
8         0          0          0       0       0          0          0       1       0
10        0          0          0       0       1          0          0       0       0
11        0          0          0       0       0          0          0       0       1
12        0          0          0       0       1          0          0       0       0
13        0          0          0       0       0          0          0       0       1
15        0          1          0       0       0          0          0       0       0
16        0          0          0       0       0          1          0       0       0
17        0          1          0       0       0          0          0       0       0
18        0          0          0       0       0          1          0       0       0
attr(,"assign")

v <- voom(x1,design=design,normalize="quantile",plot=TRUE)
v <- voom(x1,design=design,normalize="none",plot=TRUE)
colnames(v)=targets1$sample

fit <- lmFit(v, design)
colnames(design)
[1] "TEC_E13.5"  "cTEC_E15.5" "cTEC_E17.5" "cTEC_NB"    "cTEC_4W"    "mTEC_E15.5" "mTEC_E17.5" "mTEC_NB"   
[9] "mTEC_4W" 

#Look at Different times within cTEC group
cont <- makeContrasts( "cTEC_E15.5-TEC_E13.5",
                       "cTEC_E17.5-cTEC_E15.5",
                       "cTEC_NB-cTEC_E17.5",
                       "cTEC_4W-cTEC_NB",
                       "mTEC_E15.5-TEC_E13.5",
                       "mTEC_E17.5-mTEC_E15.5",
                       "mTEC_NB-mTEC_E17.5",
                       "mTEC_4W-mTEC_NB",
                       Dif17.5d =(cTEC_E17.5-cTEC_E15.5)-(mTEC_E17.5-mTEC_E15.5),
                       DifNB =(cTEC_NB-cTEC_E15.5)-(mTEC_NB-mTEC_E15.5),
                       Dif4W =(cTEC_4W-cTEC_NB)-(mTEC_4W-mTEC_NB),
                       levels=design) 
fit2 <- contrasts.fit(fit, cont) 
fit2 <- eBayes(fit2) 
topTableF(fit2, adjust="BH")


FC = 2^(fit2$coefficients)
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC) = paste(colnames(FC),"FC",sep=".")
pval = fit2$p.value
colnames(pval)=paste(colnames(pval),"pval",sep=".")
adjpval=apply(fit2$p.value,2,function(x) p.adjust(x, method='BH'))
colnames(adjpval)=paste(colnames(pval),"adj.pval",sep=".")
colorder=targets1[order(targets1$Time3,targets1$Type),"sample"]
head(v$E[,colorder])
vE_order=v$E[,colorder]
finalres=cbind(vE_order,FC,pval,adjpval)

write.table(finalres,file="../results/RNASeq_allresults.txt",sep="\t")

```

####*After Normalization*

```{r, echo=FALSE,webgl=TRUE}

colnames(v$E)=targets1$sample
#colnames(v$E)=paste(fc1filtnames,targets$batch[n])
df.m <- melt(as.data.frame(v$E))

ggplot(df.m, aes(x=value, group=variable)) + 
  geom_density(aes(colour = variable, linetype = variable),
               size=1, alpha=.60, kernel='epanechnikov') +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position='right') + scale_x_log10() + ggtitle("Raw Counts") +
  scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),8))

dev.off()

par(mar=c(10,7,1,1))
boxplot(value~variable,las=2,data=df.m,main="Normalized Signal", 
        ylab="Counts",col=as.numeric(as.factor(celltype1)))

dev.off()

par(mfrow=c(2,4))
for (i in 1:dim(v)[2]){
  plotMD(v, column=i)
  abline(h=0, col="red", lty=2, lwd=2)
}


#v2=normalizeCyclicLoess(v, weights = NULL, 
#          span=0.7, iterations = 3, method = "fast")

edf=as.matrix(v$E)
mTEC=rownames(Samp.info[(Samp.info$cell=="cTEC" & Samp.info$Time != "Myc"),])
edf=as.matrix(v$E[,mTEC])
tedf= t(edf)
pca=prcomp(tedf,scale.=T)
tedf1 = data.frame(tedf)
#Phenotype=targets1$Typetime
#cell_rep=targets1$sample
Phenotype=t
Phenotype=Samp.info[mTEC,]$Time
Phenotype=paste(c,t,sep = ".")
tedf1$group = as.factor(Phenotype)

plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
#pca$x[,1:3]  #look at first 3 PC's
summary(pca)

plot(cumsum(pca$sdev^2/sum(pca$sdev^2)))  #Plot %variance explained

#Another way:
vars <- apply(pca$x, 2, var) 
props <- vars / sum(vars)
plot(cumsum(props))

#Plot 3D PCA:
plot3d(pca$x[,1:3],col = as.integer(tedf1$group),type="s",size=2)
#group.v<-as.vector(cell_rep)
group.v <- as.vector(paste(Samp.info$cell,Samp.info$Time,
                           Samp.info$Batch,sep="_"))
group.v <- as.vector(paste(Samp.info[mTEC,]$cell,Samp.info[mTEC,]$Time,
                           Samp.info[mTEC,]$Batch,sep="_"))
text3d(pca$x[,1], pca$x[,2], pca$x[,3], group.v, cex=1.0, adj = 1.5) 
rgl.postscript("/Users/maggiec/GitHub/Maggie/ccbr704/results/pca3d_norm.pdf","pdf")
rgl.postscript("/Users/maggiec/GitHub/Maggie/ccbr704/results/pca3d_norm2.pdf","pdf")  #remove outlier

```

#####*Similarity Matrix*########

library(amap)
library(lattice) 
library(IDPmisc)
par(mar=c(2,1,1,1))

d=Dist(tedf,method="pearson",diag=TRUE)
m=as.matrix(d)

new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb")
#levelplot(m[1:ncol(m),ncol(m):1],col.regions=new.palette(20))
heatmap(m,symm=TRUE,col=new.palette(20),margins = c(7,7))

cTEC=targets1[order(targets1$Type,targets1$Time3),][c(1:7),]
cTEC_samp=v$E[,colnames(v$E) %in% cTEC$sample]
colnames(cTEC_samp)
[1] "cTEC_E17.5_1"
[2] "cTEC_E17.5_2" 
[3] "cTEC_NB_1"    
[4] "cTEC_NB_2"    
[5] "cTEC_4W_2"    
[6] "cTEC_E15.5_1"
[7] "cTEC_E15.5_2"

ipairs(cTEC_samp, pixs=0.5,cex.diag=0.5)
ipairs(cTEC_samp[,c(1,2,6,7)], pixs=0.5)
```


####*Statistical Analysis of Experimental groups (lmFit)*

```{r, echo=FALSE, warning=FALSE, message=FALSE}
#Run analysis of Experimental group:
targets1
                         bam  Time Type Rep   Typetime Time2       sample Time3
2    Sample_11_E17_5_cTEC.bam E17.5 cTEC   1 cTEC_E17.5 E17.5 cTEC_E17.5_1     2
3  Sample_12_E17_5_mTEC_1.bam E17.5 mTEC   1 mTEC_E17.5 E17.5 mTEC_E17.5_1     2
4  Sample_13_E17_5_cTEC_2.bam E17.5 cTEC   2 cTEC_E17.5 E17.5 cTEC_E17.5_2     2
5  Sample_14_E17_5_mTEC_2.bam E17.5 mTEC   2 mTEC_E17.5 E17.5 mTEC_E17.5_2     2
6     Sample_15_NB_cTEC_1.bam    NB cTEC   1    cTEC_NB    NB    cTEC_NB_1     3
7     Sample_16_NB_mTEC_1.bam    NB mTEC   1    mTEC_NB    NB    mTEC_NB_1     3
8     Sample_17_NB_cTEC_2.bam    NB cTEC   2    cTEC_NB    NB    cTEC_NB_2     3
9     Sample_18_NB_mTEC_2.bam    NB mTEC   2    mTEC_NB    NB    mTEC_NB_2     3
11   Sample_1_E13_5_TEC_1.bam E13.5  TEC   1  TEC_E13.5 E13.5  TEC_E13.5_1     0
13      Sample_21_4W_mTEC.bam    4W mTEC   1    mTEC_4W    4W    mTEC_4W_1     4
14      Sample_22_5W_cTEC.bam    5W cTEC   2    cTEC_4W    4W    cTEC_4W_2     4
15      Sample_23_5W_mTEC.bam    5W mTEC   2    mTEC_4W    4W    mTEC_4W_2     4
17   Sample_2_E13_5_TEC_2.bam E13.5  TEC   2  TEC_E13.5 E13.5  TEC_E13.5_2     0
19  Sample_4_E15_5_cTEC_1.bam E15.5 cTEC   1 cTEC_E15.5 E15.5 cTEC_E15.5_1     1
20  Sample_5_E15_5_mTEC_1.bam E15.5 mTEC   1 mTEC_E15.5 E15.5 mTEC_E15.5_1     1
22  Sample_7_E15_5_cTEC_2.bam E15.5 cTEC   2 cTEC_E15.5 E15.5 cTEC_E15.5_2     1
23  Sample_8_E15_5_mTEC_2.bam E15.5 mTEC   2 mTEC_E15.5 E15.5 mTEC_E15.5_2     1

library(splines) 
class(targets1$Time3) = "numeric"
targets1$Type2=targets1$Type
targets1$Type2[targets1$Type=="TEC"&targets1$Rep=="1"]="mTEC"
targets1$Type2[targets1$Type=="TEC"&targets1$Rep=="2"]="cTEC"

#Select 1 from each group
n1=c(-3,-4,-7,-8,-10,-16,-17)
targets2=targets1[n1,]
targets2

                          bam  Time Type Rep   Typetime Time2       sample Time3 Type2
2    Sample_11_E17_5_cTEC.bam E17.5 cTEC   1 cTEC_E17.5 E17.5 cTEC_E17.5_1     2  cTEC
3  Sample_12_E17_5_mTEC_1.bam E17.5 mTEC   1 mTEC_E17.5 E17.5 mTEC_E17.5_1     2  mTEC
6     Sample_15_NB_cTEC_1.bam    NB cTEC   1    cTEC_NB    NB    cTEC_NB_1     3  cTEC
7     Sample_16_NB_mTEC_1.bam    NB mTEC   1    mTEC_NB    NB    mTEC_NB_1     3  mTEC
11   Sample_1_E13_5_TEC_1.bam E13.5  TEC   1  TEC_E13.5 E13.5  TEC_E13.5_1     0  mTEC
14      Sample_22_5W_cTEC.bam    5W cTEC   2    cTEC_4W    4W    cTEC_4W_2     4  cTEC
15      Sample_23_5W_mTEC.bam    5W mTEC   2    mTEC_4W    4W    mTEC_4W_2     4  mTEC
17   Sample_2_E13_5_TEC_2.bam E13.5  TEC   2  TEC_E13.5 E13.5  TEC_E13.5_2     0  cTEC
19  Sample_4_E15_5_cTEC_1.bam E15.5 cTEC   1 cTEC_E15.5 E15.5 cTEC_E15.5_1     1  cTEC
20  Sample_5_E15_5_mTEC_1.bam E15.5 mTEC   1 mTEC_E15.5 E15.5 mTEC_E15.5_1     1  mTEC

library(splines)
X <- ns(targets2$Time3, df=4)

Group <- factor(targets2$Type2,levels=c("mTEC","cTEC")) 
design <- model.matrix(~Group*X) 
design

design <- model.matrix(~Group:X)
design <- design[,colSums(abs(design))>0]
colnames(design)=make.names(colnames(design))

vfilt=v[,n1]
vfilt$targets
fit <- lmFit(vfilt, design) 
cont.mat <- makeContrasts(cTEC1=GroupcTEC.X2-GroupcTEC.X1, 
                          cTEC2=GroupcTEC.X3-GroupcTEC.X2,
                          cTEC2=GroupcTEC.X4-GroupcTEC.X3,
                          mTEC1=GroupmTEC.X2-GroupmTEC.X1,
                          mTEC2=GroupcTEC.X3-GroupcTEC.X2,
                          mTEC3=GroupcTEC.X4-GroupcTEC.X3,
                          levels=design)

fit2 <- contrasts.fit(fit,cont.mat)
fit2 <- eBayes(fit2)
topTable(fit2)
head(fit2$coeff[order(fit2$coeff[,1]),])

v$targets
         group lib.size norm.factors
cTEC_E17.5_1     1  6109531    0.9547239
mTEC_E17.5_1     1  6657359    0.9431515
cTEC_E17.5_2     1  7403380    1.0067125
mTEC_E17.5_2     1  7512099    0.9733896
cTEC_NB_1        1  7189729    0.9915664
mTEC_NB_1        1  7050544    0.9189894
cTEC_NB_2        1  8229733    0.9918845
mTEC_NB_2        1  7303980    0.9778963
TEC_E13.5_1      1  7147272    1.0296538
mTEC_4W_1        1  6361322    0.9451205
cTEC_4W_2        1  6889820    1.0156583
mTEC_4W_2        1  6697606    0.9321881
TEC_E13.5_2      1  7504929    1.1097570
cTEC_E15.5_1     1  6915244    0.9938911
mTEC_E15.5_1     1  7942478    0.9985528
cTEC_E15.5_2     1  8298503    1.0415438
mTEC_E15.5_2     1  7879525    0.9818616

plot(as.factor(rownames(v$targets)),
 v$targets$lib.size,las=2,cex.axis=0.5,ylim=c(6000000,9000000))

barplot(v$targets$lib.size,main = "library size",
        names.arg = as.factor(rownames(v$targets)),
        cex.names = 0.7, las=2,
        col = "red")

##############
options(digits=3) 

design
fit <- lmFit(v,design) 
colnames(fit)
"cTEC_E15.5" 
"cTEC_E17.5" 
"cTEC_NB"    
"cTEC_4W"    

"mTEC_E15.5" 
"mTEC_E17.5" 
"mTEC_NB"    
"mTEC_4W"   



#Look at Different times within mTEC group:





#Look at Difference of Differences:
cont.dif <- makeContrasts( + Dif6hr =(mu.6hr-mu.0hr)-(wt.6hr-wt.0hr), 
                          + Dif24hr=(mu.24hr-mu.6hr)-(wt.24hr-wt.6hr), 
                          + levels=design) > fit2 <- contrasts.fit(fit, cont.dif) 
fit2 <- eBayes(fit2) 
topTableF(fit2, adjust="BH")





finalres = cbind(,FC[,4:6],pval[,4:6],adjpval[,4:6])

#finalres = cbind(EBR,FC,pval,adjpval)
head(finalres)

write.table(finalres,file="../../results/RNASeq_allbatches_results2.txt",sep="\t")



########Dendrogram ###

library(dendextend)
library(dendextendRcpp)

dist = dist(tedf)
dend = hclust(dist)
#dend = color_branches(dend, k = 2)
dend <- tedf %>% dist %>% hclust(method = "average") %>% as.dendrogram
dend <- dend %>% set("branches_lwd", 0.5)
dend <- dend %>% set("branches_k_color", k = 4) %>%
  plot(main = "Dendrogram for all genes")

###Do K means clustering

library(cluster) 

#Silhouette plot:
kfit <- kmeans(tedf, 4)
dissE <- daisy(tedf) 
dE2   <- dissE^2
sk2   <- silhouette(kfit$cl, dE2)
plot(sk2)

library(factoextra)
library(fpc)
library(NbClust)

#Kmeans plot

fviz_nbclust(nb) + theme_minimal()


plot(tedf, col=kfit$cluster)
points(kfit$centers, col = 1:4, pch = 8)

km.res <- eclust(tedf, "kmeans", k = 4,
                 nstart = 25, graph = FALSE)




library(amap)
kfit  <- Kmeans(tedf, 4, method = "euclidean")
plot(tedf, col = kfit$cluster)

##############

library(GMD)
png(file=paste("heatmap_PC2","png",sep="."),
    width=3.25,height=3.25,units="in",res=3000,pointsize = 5)
main_title="Heatmap of top PC1 genes"
main_title="Heatmap of top PC2 genes"
par(cex.main=1)

ht_global_opt(annotation_legend_title_gp = gpar(fontsize = 5), 
              annotation_legend_labels_gp = gpar(fontsize = 3))
df = data.frame(Sample = targets$Sample[order(targets$order)])
ha2=HeatmapAnnotation(df, 
                      col = list(Sample = c("CHILP" =  "darkorange2", 
                                            "CLP" = "cyan3", "EILP" = "red", "ETP" = "blue")),
                      width = unit(0.2, "mm"),annotation_width = 1.0)

hm=Heatmap(t(x.mat.ord), col = colorRamp2(c(-3, 0, 3), 
                                          c("blue", "white", "red")), 
           cluster_rows = TRUE, cluster_columns = FALSE,
           clustering_distance_rows = "euclidean",
           row_names_gp = gpar(fontsize = 2),
           column_names_gp = gpar(fontsize = 2),
           column_title = main_title, column_title_gp = 
             gpar(fontsize = 5, fontface = "bold"),
           heatmap_legend_param = list(title = "legend",
                                       width = unit(0.1, "mm")),
           bottom_annotation = ha2, 
           bottom_annotation_height = unit(3, "mm"))

draw(hm,show_heatmap_legend = FALSE, show_annotation_legend = TRUE)
dev.off()



###Draw Heatmaps of top TF genes (EILP vs CLP)
edf=EBR[rownames(EBR) %in% selectgenes,] #

finalres=as.data.frame(finalres)

[1,] "CHILP.1"                    
[2,] "CHILP.2"                    
[3,] "CHILP.3"                    
[4,] "CLP.10K.2"                  
[5,] "CLP.300.2"                  
[6,] "CLP.1"                      
[7,] "EILP.1"                     
[8,] "EILP.KO1.3"                 
[9,] "EILP.WT.3"                  
[10,] "EILP.2"                     
[11,] "ETP3.1"                     
[12,] "ETP.WT5.3"                  
[13,] "ETP.WT6.3"                  
[14,] "ETP.WT1.3"                  
[15,] "ETP.WT2.3"                  
[16,] "ETP.WT3.3"                  
[17,] "ETP.WT4.3"                  
[18,] "celltypeETP.FC"             
[19,] "celltypeEILP.FC"            
[20,] "celltypeCHILP.FC"           
[21,] "celltypeETP.pval"           
[22,] "celltypeEILP.pval"          
[23,] "celltypeCHILP.pval"         
[24,] "celltypeETP.pval.adj.pval"  
[25,] "celltypeEILP.pval.adj.pval" 
[26,] "celltypeCHILP.pval.adj.pval"

EILP_TFgenes=finalres[(finalres$celltypeEILP.FC>3) & (finalres$celltypeEILP.pval<0.01),
                      c(4:10,1:3)]
EILP_TFgenes=EILP_TFgenes[rownames(EILP_TFgenes) %in% selectgenes,]
dim(EILP_TFgenes)
[1] 47 10
head(EILP_TFgenes)
EILP_TFgenes=as.matrix(EILP_TFgenes)
y=EILP_TFgenes
y.mat=apply(y, 1, function(y) scale (y,center = TRUE, scale = TRUE))

library(GMD)
library(ComplexHeatmap)
library(circlize)

png(file=paste("heatmap_EILPvsCLP","png",sep="."),
    width=3.25,height=3.25,units="in",res=3000,pointsize = 5)
par(cex.main=1)

ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 10), 
              annotation_legend_labels_gp = gpar(fontsize = 1))


df1 = data.frame(Celltypes=celltype[4:6]) #CLP
df2 = data.frame(Celltypes=celltype[7:10]) #EILP
df3 = data.frame(Celltypes=celltype[1:3]) #CHILP

ha_col1=columnAnnotation(df1, 
                         col = list(Celltypes = c("CLP" =  "darkorange2")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col2=columnAnnotation(df2, 
                         col = list(Celltypes = c("EILP" = "cyan3")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col3=columnAnnotation(df3, 
                         col = list(Celltypes = c("CHILP" = "blue")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)

#Get Row Order
dend = dendsort(hclust(dist(t(y.mat))))

ht1=Heatmap(t(y.mat)[,1:3], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE, cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col1)
#       split = df, combined_name_fun = NULL)
ht2=Heatmap(t(y.mat)[,4:7], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE, cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col2)
ht3=Heatmap(t(y.mat)[,8:10], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            cluster_rows = dend, row_dend_reorder = FALSE,
            row_names_gp = gpar(fontsize = 3, fontface = "bold"),
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col3)

ht_list=ht1+ht2+ht3

draw(ht_list, gap = unit(1, "mm"),show_heatmap_legend = FALSE,show_annotation_legend = FALSE)

#draw(ht_list,show_heatmap_legend = FALSE, show_annotation_legend = FALSE)

dev.off()

###################

EILP_TFgenes2=finalres[(finalres$celltypeEILP.FC>3) & (finalres$celltypeEILP.pval<0.01),c(4:17)]
EILP_TFgenes2=EILP_TFgenes2[rownames(EILP_TFgenes2) %in% selectgenes,]
ETP_TFgenes=finalres[(finalres$celltypeETP.FC>3) & (finalres$celltypeETP.pval<0.01),
                     c(4:17)]
ETP_TFgenes=ETP_TFgenes[rownames(ETP_TFgenes) %in% selectgenes,]
dim(ETP_TFgenes)
[1] 39 14
head(ETP_TFgenes)
EILP_ETP_TFgenes=unique(as.matrix(rbind(EILP_TFgenes2,ETP_TFgenes)))

y=EILP_ETP_TFgenes
y.mat=apply(y, 1, function(y) scale (y,center = TRUE, scale = TRUE))

png(file=paste("heatmap_EILP_ETPvsCLP","png",sep="."),
    width=3.25,height=3.25,units="in",res=3000,pointsize = 5)
par(cex.main=1)

ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 10), 
              annotation_legend_labels_gp = gpar(fontsize = 1))

#If want to split rows:
#pa = pam(mat, k = 3)
#split = paste0("pam", pa$clustering)

df1 = data.frame(Celltypes=celltype[4:6]) #CLP
df2 = data.frame(Celltypes=celltype[7:10]) #EILP
df3 = data.frame(Celltypes=celltype[11:17]) #ETP

ha_col1=columnAnnotation(df1, 
                         col = list(Celltypes = c("CLP" =  "darkorange2")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col2=columnAnnotation(df2, 
                         col = list(Celltypes = c("EILP" = "cyan3")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col3=columnAnnotation(df3, 
                         col = list(Celltypes = c("ETP" = "green")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)

dend = dendsort(hclust(dist(t(y.mat))))

ht1=Heatmap(t(y.mat)[,1:3], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE,cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col1)
#       split = df, combined_name_fun = NULL)
ht2=Heatmap(t(y.mat)[,4:7], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE,cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col2)
ht3=Heatmap(t(y.mat)[,8:13], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            cluster_rows = dend, row_dend_reorder = FALSE,
            row_names_gp = gpar(fontsize = 3, fontface = "bold"),
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col3)

ht_list=ht1+ht2+ht3

draw(ht_list, gap = unit(1, "mm"),show_heatmap_legend = FALSE,show_annotation_legend = FALSE)

#draw(ht_list,show_heatmap_legend = FALSE, show_annotation_legend = FALSE)

dev.off()

################################

#Downregulated Genes in EILP:
neg<-function(x) -x 
finalres=as.data.frame(finalres)
EILP_TFdgenes=finalres[(finalres$celltypeEILP.FC<neg(3)) & (finalres$celltypeEILP.pval<0.01),
                       c(4:17)]
EILP_TFdgenes=EILP_TFdgenes[rownames(EILP_TFdgenes) %in% selectgenes,]
dim(EILP_TFdgenes)
[1] 175  14

head(EILP_TFdgenes)
EILP_TFdgenes=as.matrix(EILP_TFdgenes)
y=EILP_TFdgenes
y.mat=apply(y, 1, function(y) scale (y,center = TRUE, scale = TRUE))

library(GMD)
library(ComplexHeatmap)
library(circlize)

png(file=paste("heatmap_EILPvsCLP_down","png",sep="."),
    width=3.25,height=3.25,units="in",res=3000,pointsize = 5)
par(cex.main=1)

ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 10), 
              annotation_legend_labels_gp = gpar(fontsize = 1))

#If want to split rows:
#pa = pam(mat, k = 3)
#split = paste0("pam", pa$clustering)

df1 = data.frame(Celltypes=celltype[4:6]) #CLP
df2 = data.frame(Celltypes=celltype[7:10]) #EILP
df3 = data.frame(Celltypes=celltype[11:17]) #ETP

ha_col1=columnAnnotation(df1, 
                         col = list(Celltypes = c("CLP" =  "darkorange2")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col2=columnAnnotation(df2, 
                         col = list(Celltypes = c("EILP" = "cyan3")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)
ha_col3=columnAnnotation(df3, 
                         col = list(Celltypes = c("ETP" = "green")),
                         width = unit(1.0, "mm"),annotation_width = 0.5)

dend = dendsort(hclust(dist(t(y.mat))))

ht1=Heatmap(t(y.mat)[,1:3], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE,cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col1)
#       split = df, combined_name_fun = NULL)
ht2=Heatmap(t(y.mat)[,4:7], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            show_row_names = FALSE,cluster_rows = dend, row_dend_reorder = FALSE,
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col2)
ht3=Heatmap(t(y.mat)[,8:13], col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
            cluster_columns = FALSE,
            cluster_rows = dend, row_dend_reorder = FALSE,
            row_names_gp = gpar(fontsize = 2, fontface = "bold"),
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            top_annotation = ha_col3)

ht_list=ht1+ht2+ht3

draw(ht_list, gap = unit(1, "mm"),show_heatmap_legend = FALSE,show_annotation_legend = FALSE)

#draw(ht_list,show_heatmap_legend = FALSE, show_annotation_legend = FALSE)

dev.off()

##############################
#For GO analysis:

finalresgo_EILP=finalres[,c(19,21)]
write.table(finalresgo_EILP,file="../../results/EILPvsCLP.txt",sep="\t",quote = FALSE)

#############################


df = data.frame(Coregroup = cyt_genes_ord_arr$group)
ha_row=rowAnnotation(df, 
                     col = list(Coregroup = c("core1" =  "darkorange2", "core2" = "cyan3")),
                     width = unit(1.0, "mm"),annotation_width = 1.0)
ha2 = HeatmapAnnotation(df = data.frame(cytoexp = cyt_avg_heat),
                        col = list(cytoexp = colorRamp2(c(-1, 8), c("white", "red"))))
ht1=Heatmap(t(x.mat), col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")), 
            cluster_rows = TRUE, cluster_columns = FALSE,
            clustering_distance_rows = "pearson",
            row_names_gp = gpar(fontsize = 2),
            column_title = main_title, column_title_gp = 
              gpar(fontsize = 5, fontface = "bold"),
            heatmap_legend_param = list(title = "legend",width = unit(0.1, "mm")),
            split = df, combined_name_fun = NULL, top_annotation = ha2, 
            top_annotation_height = unit(3, "mm"))

ht_list=ht1+ha_row
draw(ht_list,show_heatmap_legend = FALSE, show_annotation_legend = FALSE)


```

####Calculate Scale Factors ###########

```{r,echo=FALSE}

#Get scale factors:

scaleval=(2^EBR)*100  # Real data (aiming for)
#actscale=(fc1filt/scaleval)
actscale=x$counts/scaleval
write.table(actscale,file="scaling_factors.txt", row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)

palette <- colorRampPalette(c('blue','yellow'))(12)
palette2 <- rainbow(20)

palette
color=as.numeric(as.factor(targets1$cell[n]))
f.cell=as.factor(targets1$cell[n])
topPC1=base::sort(abs(pca$rotation[,"PC1"]),decreasing=TRUE)[1:500]
v1=v$E[rownames(v$E) %in% names(topPC1),]
v1=sort(v1,)
head(v1)
dim(v1)
#hc.rows <- hclust(dist(v$E))
dev.off()

diff=v$E[,14]-v$E[,15]

topPC1=base::sort(diff,decreasing=TRUE)[1:100]
topPC1[1:50]
v1=v1[order(match(rownames(v1),names(topPC1))),]

vselect=v$E[rownames(v$E) %in% selectgenes,]
colnames(vselect)=cell_rep
dd <- as.dist((1 - cor(t(vselect))/2))
hc.rows <- hclust(dd)

heatmap(as.matrix(vselect), 
        ColSideColors = palette2[color],
        col=palette, labCol = targets1$cell,
        Rowv=as.dendrogram(hc.rows),
        #labRow=NA,
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(vselect)),
                                  method="average")))
legend("topright",legend=levels(f.cell), fill=palette2,title="Group",
       cex=0.5)

d3heatmap(vselect, Rowv = TRUE, Colv = TRUE, scale = "none", dendrogram = "none", main = "title", color = "YlOrRd",cexRow=1,cexCol=0.5,main="topgenes")

```


####*Heatmap: Top genes*

```{r, echo=FALSE, message=FALSE,warning=FALSE}
library(d3heatmap)
library(reshape2)

topgenes=topTable(fit,coef="celltypeETP",n=50,sort.by="p")
topgenes=topTable(fit,coef="celltypeEILP",n=50,sort.by="p")
topgenes=topTable(fit,coef="celltypeCHILP",n=50,sort.by="p")

topgenes=selectgenes
topgenes_data=v$E[rownames(v$E) %in% selectgenes,]
colnames(topgenes_data)=colnames(v$E)
topgenes_data=topgenes_data[match(topgenes,rownames(topgenes_data)),]

#Original Heatmap (notice batch effect)
topgenes_data=v$E[rownames(v$E) %in% rownames(topgenes),]
colnames(topgenes_data)=colnames(v$E)
topgenes_data=topgenes_data[match(rownames(topgenes),rownames(topgenes_data)),]
#colnames(topgenes_data)=targets$Phenotype[match(colnames(topgenes_data),targets$bam)]
#colnames(topgenes_data)=paste(colnames(topgenes_data),targets$Replicate,sep=".")
#cc=colsplit(string=colnames(topgenes_data), pattern="_", names=c("Part1", "Part2"))
#dd=colsplit(string=cc$Part2, pattern="\\.", names=c("Part1", "Part2"))
#colnames(topgenes_data)=dd$Part1
d3heatmap(topgenes_data, scale = "row", dendrogram = "both", main = "title",
          color = "RdYlBu",cexRow=0.8,main="topgenes")
```


#Look at genes##

bed=read.table("mm10-TCRsIgs2.bed",header=FALSE)
colnames(bed)=c("chr","start","stop","name")
bed[bed$start %in% fc$annotation$Start,]

Gene=fc$annotation[fc$annotation$Start %in% bed$start,]
Genecount=fc$counts[rownames(fc$counts) %in% Gene,]
merge(fc$annotation,bed,by.x=fc$annotation$Start,by.y=bed$start)

#####Run TSCAN analysis ########
colnames(v$E)
[1] "cTEC_E17.5_1" "mTEC_E17.5_1" "cTEC_E17.5_2" "mTEC_E17.5_2" "cTEC_NB_1"    "mTEC_NB_1"   
[7] "cTEC_NB_2"    "mTEC_NB_2"    "TEC_E13.5_1"  "mTEC_4W_1"    "cTEC_4W_2"    "mTEC_4W_2"   
[13] "TEC_E13.5_2"  "cTEC_E15.5_1" "mTEC_E15.5_1" "cTEC_E15.5_2" "mTEC_E15.5_2"

library(TSCAN)
library(raster)
ALL_LP=v$E[,c(9,13,14,16,1,3,5,7,11)] #cTEC
ALL_LP=v$E[,c(9,13,14,16,3,5,7,11)] #cTEC minus #E17.5_1
head(ALL_LP)

dim(ALL_LP[rownames(ALL_LP) %in% remdata=="FALSE",])
[1] 17142     9
> dim(x_filt$filteredData)
[1] 17465    17
ALL_LP_filt=ALL_LP[rownames(ALL_LP) %in% remdata=="FALSE",]
dim(ALL_LP_filt)

ALL_LP2=v$E[,c(9,13,15,17,2,4,6,8,10,12)]  #mTEC
head(ALL_LP2)
dim(ALL_LP2[rownames(ALL_LP2) %in% remdata=="FALSE",])
[1] 17142     9
dim(x_filt$filteredData)
[1] 17465    17
ALL_LP2_filt=ALL_LP2[rownames(ALL_LP2) %in% remdata=="FALSE",]


procdata <- preprocess(ALL_LP)

logdata=ALL_LP
logdata=ALL_LP_filt    #cTEC

logdata=ALL_LP2
logdata=ALL_LP2_filt   #mTEC

logdatafilt=apply(logdata,1,function(x) length(x[x>=1])>=1)
filtdata=logdata[logdatafilt,]
dim(filtdata)
cvlogdatafilt=apply(filtdata,1,function(x) cv(2^x))
hist(cvlogdatafilt,breaks=100)
quantile(cvlogdatafilt, c(.50, .75, .90))
topcv=quantile(cvlogdatafilt, .25)
dim(filtdata[cvlogdatafilt>topcv,])
procdata=filtdata[cvlogdatafilt>topcv,]
#procdata=ALL_LP
dim(procdata)

library(rgl)
pca=prcomp(t(procdata),scale.=T) 
group.v=sapply(strsplit(colnames(procdata),split="\\."),function(x) x[1])

plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
plot3d(pca$x[,1:3],col = as.integer(as.factor(group.v)), type="s",size=2)
text3d(pca$x[,1], pca$x[,2], pca$x[,3], colnames(procdata), cex=1.0, adj = 1.2)
rgl.postscript("/Users/maggiec/GitHub/Maggie/ccbr546/results/pca3d_ALL_LP.pdf","pdf")


lpmclust <- exprmclust(procdata,clusternum = 5)
lpmclust <- exprmclust(procdata[1:1000,],clusternum = 3)
lpmclust
plotmclust(lpmclust)
lporder <- TSCANorder(lpmclust)
lporder <- TSCANorder(lpmclust,flip=T,listbranch = T)
lporder

##Hack into exprmclust:

library(mclust)
library(igraph)

data=logdata # use HTSFilt
clusternum=9
clusternum=8 #cTEC

data=procdata
clusternum=8

data=filtdata
clusternum=8

#function (data, clusternum = 2:9, modelNames = "VVV", reduce = T) 
#{
  set.seed(12345)
#  if (reduce) {
    sdev <- prcomp(t(data), scale = T)$sdev[1:20]
    x <- 1:20
    optpoint <- which.min(sapply(2:10, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(sdev ~ x + x2)$residuals^2)
    }))
    pcadim = optpoint + 1
    tmpdata <- t(apply(data, 1, scale))
    colnames(tmpdata) <- colnames(data)
    tmppc <- prcomp(t(tmpdata), scale = T)
    pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
 # }
#  else {
#    pcareduceres <- t(data)
#  }
  clusternum <- clusternum[clusternum > 1]
#  res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, 
#                                 modelNames = "VVV"))  #didn't work
  res <- suppressWarnings(Mclust(pcareduceres, G = clusternum))
  clusterid <- apply(res$z, 1, which.max)
  clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = res$G)
  for (cid in 1:res$G) {
    clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == 
                cid]), , drop = F])
  }
  dp <- as.matrix(dist(clucenter))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, 
       clucenter = clucenter)
#}
  igraph::distances(gp)  # for mTEC:
  1        2        3        4        5        6        7        8        9
  1   0.0000 151.6859 154.4178 190.0598 187.7585 203.2636 220.7357 197.2951 211.8080
  2 151.6859   0.0000 142.8621 171.6360 174.2993 184.7742 202.8030 183.4527 194.6346
  3 154.4178 142.8621   0.0000 171.2698 168.6025 181.5835 202.5401 180.1024 190.7491
  4 190.0598 171.6360 171.2698   0.0000 170.4632 170.3392 188.8360 175.4498 179.8285
  5 187.7585 174.2993 168.6025 170.4632   0.0000 174.6412 189.1884 174.4176 182.4038
  6 203.2636 184.7742 181.5835 170.3392 174.6412   0.0000 164.4923 157.4724 160.5678
  7 220.7357 202.8030 202.5401 188.8360 189.1884 164.4923   0.0000 159.0483 166.5861
  8 197.2951 183.4527 180.1024 175.4498 174.4176 157.4724 159.0483   0.0000 146.2239
  9 211.8080 194.6346 190.7491 179.8285 182.4038 160.5678 166.5861 146.2239   0.0000 


lpmclust=list(pcareduceres = pcareduceres, 
        MSTtree = dp_mst, clusterid = clusterid,
        clucenter = clucenter)
lpmclust
plotmclust(lpmclust,cell_name_size = 5)
plot(lpmclust$MSTtree)

mclustobj=lpmclust


function (mclustobj, MSTorder = NULL, orderonly = T, flip = F, 
          listbranch = F) 
{
  if (!is.null(MSTorder) & length(MSTorder) == 1) {
    stop("MSTorder is not a path!")
  }
  set.seed(12345)
  clucenter <- mclustobj$clucenter
  row.names(clucenter) <- paste0("clu", 1:nrow(clucenter))
  clusterid <- mclustobj$clusterid
  pcareduceres <- mclustobj$pcareduceres
  adjmat <- as_adjacency_matrix(mclustobj$MSTtree, sparse = FALSE)
  if (is.null(MSTorder)) {
    orderinMST <- 1
    clutable <- table(mclustobj$clusterid)
    alldeg <- degree(mclustobj$MSTtree)
    allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg == 1]), 
                           as.numeric(names(alldeg)[alldeg == 1]))
    allcomb <- allcomb[allcomb[, 1] < allcomb[, 2], ]
    numres <- t(apply(allcomb, 1, function(i) {
      tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree, 
                                          i[1], i[2])$vpath[[1]])
      c(length(tmp), sum(clutable[tmp]))
    }))
    optcomb <- allcomb[order(numres[, 1], numres[, 2], decreasing = T)[1], 
                       ]
    branchcomb <- allcomb[-order(numres[, 1], numres[, 2], 
                                 decreasing = T)[1], ]
    MSTorder <- get.shortest.paths(mclustobj$MSTtree, optcomb[1], 
                                   optcomb[2])$vpath[[1]]
    if (flip) {
      MSTorder <- rev(MSTorder)
    }
  }
  else {
    edgeinMST <- sapply(1:(length(MSTorder) - 1), function(i) {
      adjmat[MSTorder[i], MSTorder[i + 1]]
    })
    if (sum(edgeinMST == 0) > 0) {
      orderinMST <- 0
    }
    else {
      orderinMST <- 1
    }
  }

#Set input variables:  
  MSTinout=1
  internalorder=MSTorder
  
  internalorderfunc <- function(internalorder, MSTinout) {
    TSCANorder <- NULL
    for (i in 1:(length(internalorder) - 1)) {
      currentcluid <- internalorder[i]
      #> currentcluid
      #+ 1/3 vertex, named:
      #  [1] 1
      nextcluid <- internalorder[i + 1]
      #+ 1/3 vertex, named:
      #  [1] 2
      currentclucenter <- clucenter[currentcluid, ]
      #[1] -41.6935672  61.8574257  13.6194363  -7.4981269  -0.3347489   0.2865341   0.1581459
      nextclucenter <- clucenter[nextcluid, ]
      #[1] -16.6812933 -58.4265349  -1.6693466  -1.1513396  -0.9064566  -0.8055364  -0.7048275
      currentreduceres <- pcareduceres[clusterid == currentcluid, ]
      #currentreduceres
      #PC1      PC2       PC3       PC4       PC5       PC6        PC7
      #CLP.1 -45.62629 62.12603  7.678426 -12.58850  78.84035 -18.70177 -11.148273
      #CLP.2 -52.12689 64.22814 31.937031  57.99704 -54.65701  -9.96409   1.341105
      #CLP.3 -27.32753 59.21810  1.242852 -67.90292 -25.18758  29.52546  10.281607
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[currentcluid, ] == 1)))
      }
      #connectcluid
      #[1] 2
      else {
        if (i == 1) {
          connectcluid <- nextcluid
        }
        else {
          connectcluid <- c(nextcluid, internalorder[i - 1])
        }
      }
      cludist <- sapply(connectcluid, function(x) {
        rowSums(sweep(currentreduceres, 2, clucenter[x, ], "-")^2)
      })
      #x=connectcluid
      #rowSums(sweep(currentreduceres, 2, clucenter[x, ], "-")^2)
      CLP.1    CLP.2    CLP.3 
      22377.82 23905.67 20048.10
      TEC_E13.5_1 TEC_E13.5_2 
      26370.09    28859.79 
      #cludist
      #[,1]
      #CLP.1 22377.82
      #CLP.2 23905.67
      #CLP.3 20048.10
      #cludist
      #[,1]
      #TEC_E13.5_1 26370.09
      #TEC_E13.5_2 28859.79
      mindistid <- apply(cludist, 1, which.min)
      #mindistid
      #CLP.1 CLP.2 CLP.3 
      #1     1     1
      #TEC_E13.5_1 TEC_E13.5_2 
      #1           1 
      edgecell <- names(which(mindistid == which(connectcluid == nextcluid)))
      #edgecell
      #[1] "CLP.1" "CLP.2" "CLP.3"
      
      difvec <- nextclucenter - currentclucenter
      #difvec
      #[1]   25.0122740 -120.2839607  -15.2887829    6.3467873   -0.5717077   -1.0920705   -0.8629735
      tmppos <- pcareduceres[edgecell, ] %*% difvec
      #tmppos
      #PC1      PC2       PC3       PC4       PC5       PC6        PC7
      #CLP.1 -45.62629 62.12603  7.678426 -12.58850  78.84035 -18.70177 -11.148273
      #CLP.2 -52.12689 64.22814 31.937031  57.99704 -54.65701  -9.96409   1.341105
      #CLP.3 -27.32753 59.21810  1.242852 -67.90292 -25.18758  29.52546  10.281607
      [,1]
      TEC_E13.5_1 -11831.07
      TEC_E13.5_2 -13075.92
      
      pos <- as.vector(tmppos)
      #pos
      #CLP.1     CLP.2     CLP.3 
      #-8826.302 -9108.639 -8283.195 
      #TEC_E13.5_1 TEC_E13.5_2 
      #-11831.07   -13075.92 
      names(pos) <- row.names(tmppos)
      #[1] "CLP.1" "CLP.2" "CLP.3"
      [1] "TEC_E13.5_1" "TEC_E13.5_2"
      TSCANorder <- c(TSCANorder, names(sort(pos)))
      #TSCANorder
      #[1] "CLP.2" "CLP.1" "CLP.3"
      [1] "TEC_E13.5_2" "TEC_E13.5_1"
      nextreduceres <- pcareduceres[clusterid == nextcluid,]
      #nextreduceres
      #PC1       PC2        PC3        PC4        PC5       PC6       PC7
      #EILP.1  -4.179801 -70.06032  74.554485  29.105486  26.315740  43.68861  14.00482
      #EILP.2 -14.511074 -54.93478  -9.406074  -4.025238 -14.056735 -17.81326 -74.65124
      #EILP.3 -15.307204 -57.63337  -8.073171 -10.462428  -5.070884 -65.17300  45.84909
      #EILP.4 -32.727095 -51.07767 -63.752626 -19.223179 -10.813948  36.07550  11.97803
      PC1        PC2        PC3        PC4        PC5        PC6        PC7        PC8        PC9 
      -50.418252 -29.201625 -26.565607  20.082463  55.873599 -17.099986 -20.045339 -66.585049  -6.888034 
      
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[nextcluid, ] == 1)))
      }
      else {
        if (i == length(internalorder) - 1) {
          connectcluid <- currentcluid
        }
        else {
          connectcluid <- c(currentcluid, internalorder[i + 2])
        }
      }
      cludist <- sapply(connectcluid, function(x) {
        rowSums(sweep(nextreduceres, 2, clucenter[x,], "-")^2)
      })
      connectcluid
      [1] 1 3
      #clucenter[x,]
      #[,1]     [,2]      [,3]      [,4]       [,5]      [,6]      [,7]
      #clu1 -41.69357 61.85743  13.61944 -7.498127 -0.3347489 0.2865341 0.1581459
      #clu3  95.90294 24.06693 -17.09046 13.549870  2.3150366 1.1812716 1.1724362
      clucenter[x,]
      [,1]      [,2]       [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
      clu1 -109.18498  44.18544   6.716776 -14.91665 -21.80898  6.425824 -1.927514  3.534904 -4.977807
      clu3  -44.36303 -35.77149 -13.849342  38.11689  48.87898 -2.626823 13.387579 65.826009 23.672822
      
      rowSums(sweep(nextreduceres, 2, clucenter[x,], "-")^2)
      EILP.1   EILP.2   EILP.3   EILP.4 
      33455.36 10228.95 35859.39 11012.67 
      Warning message:
        In sweep(nextreduceres, 2, clucenter[x, ], "-") :
        STATS is longer than the extent of 'dim(x)[MARGIN]'
      #cludist
      #[,1]     [,2]
      #EILP.1 26648.20 30064.84
      #EILP.2 21033.88 25178.52
      #EILP.3 21848.69 26153.98
      #EILP.4 20489.02 26950.57
      
      mindistid <- apply(cludist, 1, which.min)
      edgecell <- names(which(mindistid == which(connectcluid == 
              currentcluid)))
      difvec <- nextclucenter - currentclucenter
      tmppos <- pcareduceres[edgecell, ] %*% difvec
      pos <- as.vector(tmppos)
      names(pos) <- row.names(tmppos)
      TSCANorder <- c(TSCANorder, names(sort(pos)))
    }
    if (orderonly) {
      TSCANorder
    }
    else {
      data.frame(sample_name = TSCANorder, State = clusterid[TSCANorder], 
                 Pseudotime = 1:length(TSCANorder), stringsAsFactors = F)
    }
  }
  if (!orderinMST) {
    internalorderfunc(MSTorder, 0)
  }
  else {
    if (exists("branchcomb") & listbranch) {
      allres <- list()
      allres[[paste("backbone", paste(MSTorder, collapse = ","))]] <- 
        internalorderfunc(MSTorder, 1)
      for (tmpcombid in 1:nrow(branchcomb)) {
        tmporder <- get.shortest.paths(mclustobj$MSTtree, 
          branchcomb[tmpcombid, 1], branchcomb[tmpcombid, 2])$vpath[[1]]
        allres[[paste("branch:", paste(tmporder, collapse = ","))]] <- 
          internalorderfunc(tmporder, 1)
      }
      allres
    }
    else {
      internalorderfunc(MSTorder, 1)
    }
  }
}
<environment: namespace:TSCAN>function (mclustobj, MSTorder = NULL, orderonly = T, flip = F, 
                                        listbranch = F) 

function (x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) 
  {
    FUN <- match.fun(FUN)
    dims <- dim(x)
    if (check.margin) {
      dimmargin <- dims[MARGIN]
      dimstats <- dim(STATS)
      lstats <- length(STATS)
      if (lstats > prod(dimmargin)) {
        warning("STATS is longer than the extent of 'dim(x)[MARGIN]'")
      }
      else if (is.null(dimstats)) {
        cumDim <- c(1L, cumprod(dimmargin))
        upper <- min(cumDim[cumDim >= lstats])
        lower <- max(cumDim[cumDim <= lstats])
        if (lstats && (upper%%lstats != 0L || lstats%%lower != 
                       0L)) 
          warning("STATS does not recycle exactly across MARGIN")
      }
      else {
        dimmargin <- dimmargin[dimmargin > 1L]
        dimstats <- dimstats[dimstats > 1L]
        if (length(dimstats) != length(dimmargin) || any(dimstats != 
                                                         dimmargin)) 
          warning("length(STATS) or dim(STATS) do not match dim(x)[MARGIN]")
      }
    }
    perm <- c(MARGIN, seq_along(dims)[-MARGIN])
    FUN(x, aperm(array(STATS, dims[perm]), order(perm)), ...)
  }
<bytecode: 0x11d586a68>
  <environment: namespace:base>
  
  
lporder <- TSCAN::TSCANorder(lpmclust, MSTorder = c(2,3,5))
lpmclust$MSTtree
ord=c(2,1,3,4,6,5,7,9)  #order for mTEC
ord=c(2,1,3,6,7,8,9)  #order for cTEC
lporder=colnames(ALL_LP2_filt)[ord]  #mTEC
lporder=colnames(ALL_LP_filt)[ord]
lporder
[1] "TEC_E13.5_2"  "TEC_E13.5_1"  "cTEC_E15.5_1" "cTEC_E17.5_2" "cTEC_NB_1"    "cTEC_NB_2"   
[7] "cTEC_4W_2"   
class(lporder)
[1] "character"


diffval=difftest(procdata, TSCANorder(lpmclust))
diffval=difftest(procdata, TSCANorder(lpmclust,MSTorder = c(1,3,2))) #if sidebranch
diffval=difftest(ALL_LP2, TSCANorder(lpmclust))
diffval=difftest(ALL_LP2_filt, TSCANorder(lpmclust,MSTorder = c(1,3,2)))
diffval=diffval[order(diffval$pval)[diffval$qval < 0.05],]

diffgene=row.names(diffval)[diffval$qval < 0.05]
diffval[1:100,]

###start here:
procdata=data  #use HTSFilt data
diffval <- difftest(procdata,lporder)
diffval_cTEC <- difftest(procdata,lporder)
#diffval=difftest(procdata, TSCANorder(lpmclust,MSTorder = c(1,3,2)))
head(diffval[order(diffval$pval),])
diffval=diffval[order(diffval$pval),]
diffval_cTED=diffval_cTEC[order(diffval_cTEC$pval),]
head(row.names(diffval)[diffval$qval < 0.05])
plot(lpmclust$MSTtree)

class(synth_gex) #matrix

gex_cleaned <- as_data_frame(t(synth_gex)) 
%>% mutate(Pseudotime = ex_pseudotime) %>%
  gather(Gene, Expression, -Pseudotime)
class(gex_cleaned)
tail(gex_cleaned)

sde <- switchde(synth_gex, ex_pseudotime)


##Running Sigmoidal Differential Expression Models#####
cTEC_order=(2^ALL_LP_filt[,lporder])+1
dim(cTEC_order)
[1] 17142     7
cTEC_order.filt=cTEC_order[!rownames(cTEC_order) %in% filtgenes,]
dim(cTEC_order.filt)
[1] 17076     7
ptime=c(1:7)
sde_cTEC <- switchde(cTEC_order.filt, ptime)
corrgene=arrange(sde_cTEC, qval)
early=filter(corrgene,t0<=3 & k<0 & pval<0.05)
late=filter(corrgene,t0>=5 & k>0 & pval<0.05)

early_GO=early[,c(1,5,2)]
early_GO$k="2"
late_GO=late[,c(1,5,2)]
late_GO$k="-2"
other_GO=corrgene[!corrgene$gene %in% early$gene,c(1,5,2)]
other_GO=other_GO[!other_GO$gene %in% late$gene,]
other_GO$k="0"
all_GO=rbind(early_GO,late_GO,other_GO)

#setwd("SDE_cTEC")  #For all high correl genes
setwd("../SDE_early_cTEC")  #For early genes
setwd("../SDE_late_cTEC")  #For late genes
l=dim(late)[1]
for (i in 1:l){
  #gene=corrgene$gene[i]
  #gene=early$gene[i]
  gene=late$gene[i]
  pars <- extract_pars(sde_cTEC, gene)
  plot.new()
  png(file=paste(gene,"sde.png",sep="."))
  pval=signif(subset(sde_cTEC$pval,sde_cTEC$gene==gene),2)
  qval=signif(subset(sde_cTEC$qval,sde_cTEC$gene==gene),2)
  p=switchplot(cTEC_order[gene, ], ptime, pars, pval, qval, gene)
  print(p)
  dev.off()
}

switchplot <- function(x, pseudotime, pars, pval, qval, gene) {
  ggplot(data_frame(Expression = x, Pseudotime = pseudotime), 
         aes_string(x = "Pseudotime", y = "Expression")) +
    geom_point(alpha = 0.5, fill = "grey", colour = "black", shape = 21) + theme_bw() +
    stat_function(fun = sigmoid, args = list(params = pars), color = 'red') + 
    annotate("text", x = max(pseudotime)-2, y = 0.85*max(x), label = gene) +
    annotate("text", x = max(pseudotime)-2, y = 0.8*max(x), label = paste0("p=",pval)) +
    annotate("text", x = max(pseudotime)-2, y = 0.75*max(x), label = paste0("q=",qval)) 
    }

sigmoid <- function(pst, params) {
  mu0 <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- 2 * mu0 / (1 + exp(-k*(pst - t_0)))
  return(mu)
}

mTEC_order=log(2^mTEC_order+1)
corrgene=rownames(diffval[,])
sde_mTEC2 <- switchde(mTEC_order[corrgene,], as.array(ptime))
arrange(sde_mTEC, qval)

gene <- sde_mTEC$gene[which.min(sde_mTEC$qval)]

corrgene2=arrange(sde_mTEC, qval)
setwd("SDE2")
for (i in 1:500){
gene=corrgene2$gene[i]
#gene=corrgene[i]
pars <- extract_pars(sde_mTEC2, gene)
plot.new()
png(file=paste(gene,"sde2.png",sep="."))
p=switchplot(mTEC_order[gene, ], ptime, pars)
print(p)
dev.off()
}

for (i in 1:10){
pars <- extract_pars(sde_mTEC, corrgene[i])
print(pars)
png(file=paste(corrgene[i],"sde.png",sep="."),
width=3.25,height=3.25,units="in",res=3000,pointsize = 5)
switchplot(mTEC_order[corrgene[i], ], ptime, pars)
dev.off()
}

write.table(diffval,"genefit_mTEC.txt",col.names=NA,quote=FALSE,sep="\t")
write.table(diffval_cTEC,"genefit_cTEC.txt",col.names=NA,quote=FALSE,sep="\t")


#For mTEC:
q="Skint10"
q="Sphk1"   #500   very good!
q="Fabp12"   #1000  ok  pval=3.352013e-12

#For cTEC:
q="Lgals9"
q="Haus3"    #500    very good!
q="Gm15243"   #1000   ok

Expr <- ALL_LP_filt[q,ord]  #cTEC
Expr <- ALL_LP2_filt[q,ord]   #mTEC

Expr <- log2(ALL_LP["Id2",]+1)

Expr <- log2(ALL_LP["Gata3",]+1)
Expr <- log2(ALL_LP["Tox",]+1)
Expr <- log2(ALL_LP["Tcf7",]+1)
Expr <- log2(ALL_LP2["Bcl11b",]+1)
Expr <- log2(ALL_LP["Cd122",]+1)
Expr[ptime]

TFgene=c("Id2","Gata3","Tox","Tcf7","Bcl11b")

ptime=TSCANorder(lpmclust,flip=T,orderonly=FALSE)$sample_name
#for (i in 1:length(TFgene)){

for (i in 1:length(diffgene)){
  # png(file=paste("ptime.231",diffgene[i],"png",sep="."))
  file=paste("ptime.231",diffgene[i],"png",sep=".")
  Expr <- log2(ALL_LP[diffgene[i],]+1)
  Expr
  singlegeneplot(Expr, TSCANorder(lpmclust,orderonly=FALSE))
  singlegeneplot(Expr, lporder)
  plot(as.factor(lporder),Expr)
  #plot(Expr[ptime])
  singlegeneplot(Expr, TSCANorder(lpmclust,flip=T,orderonly=FALSE))
  ggsave(file, width = 5, height = 5)
}


function (geneexpr, TSCANorder, cell_size = 2) 
{
  Pseudotime <- c(1,2,3,4,5,6,7,8)
  Pseudotime <- c(1,2,3,4,5,6,7) #cTEC
  #geneexpr <- geneexpr[TSCANorder[, 1]]
  geneexpr <- Expr[lporder]
  #exprdata <- cbind(TSCANorder, geneexpr)
  exprdata <- as.data.frame(cbind(Pseudotime, geneexpr))
  #exprdata$State <- factor(exprdata$State)
  exprdata$State=as.factor(targets$Typetime[match(rownames(exprdata),targets$sample)])
  exprdata$predict <- fitted.values(mgcv::gam(geneexpr ~ s(Pseudotime, 
              k = 3), data = exprdata))
  q <- ggplot(aes(Pseudotime, geneexpr), data = exprdata)
  q <- q + geom_point(aes_string(color = "State"), size = I(cell_size))
  q <- q + geom_line(aes(Pseudotime, predict), data = exprdata)
  q <- q + ylab("Expression") + xlab("Pseudotime")
  q <- q + theme(strip.background = element_rect(colour = "white", 
       fill = "white")) + theme(panel.border = element_blank(), 
       axis.line = element_line()) + 
       theme(panel.grid.minor.x = element_blank(), 
       panel.grid.minor.y = element_blank()) + 
       theme(panel.grid.major.x = element_blank(), 
       panel.grid.major.y = element_blank()) + 
       theme(panel.background = element_rect(fill = "white"))
  q
}
<environment: namespace:TSCAN>

  

```

###*Provenance:*

```{r, echo=FALSE}

sessionInfo()

```