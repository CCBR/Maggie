Ghany

source("https://bioconductor.org/biocLite.R")

### Run in Biowulf ################
library(affy)
library(oligo)
library(genefilter)
library(limma)
library("arrayQualityMetrics")
library(stringr)
library(R2HTML)
library(pd.hugene.2.0.st)

setwd("/data/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
setwd("/Users//maggiec//Documents/Ghany/")

setwd("/Users/maggiec/GitHub/Maggie/ccbr601/data")
load(".RData")

celFiles=read.table("covdesc2",header=TRUE,fill=TRUE)[,1]

celFiles=read.table("covdesc",header=TRUE)[,1]

############################################################

# Use oligo package to read in cel files for Gene arrays

#celFiles <- list.celfiles()  

pd<-read.AnnotatedDataFrame("covdesc2",header=TRUE,sep="\t")
x <- varMetadata(pd)
x <- data.frame(x, channel = "_ALL_")
varMetadata(pd) <- x

affyData = read.celfiles(as.character(celFiles))
phenoData(affyData) = pd

colnames(affyData)
colnames(affyData) = str_match(celFiles,"(.*).CEL")[,2]

eset = rma(affyData)
e <- exprs(eset)
te <- t(e)
tedf = as.data.frame(te)
tedf=cbind(pData(eset),tedf)
tedf2 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","Timing",
                       "Failed","HCV.RNA.at.F.U.Bx.Value"), 
            function(x) colMeans(x[,24:ncol(x)]))
unique(tedf2[1:4])

tedf2_cirr_6hr=tedf2[tedf2$Treatment == "6h-IFN" & tedf2$Baseline.Fibrosis.Group == "F_5_6",]
tedf2_cirr_Base=tedf2[tedf2$Treatment == "Base" & tedf2$Baseline.Fibrosis.Group == "F_5_6" & 
                        tedf2$Timing == "6hrs post IFN",]
tedf2_noncirr_Base=tedf2[tedf2$Treatment == "Base" & 
      (tedf2$Baseline.Fibrosis.Group == "F_0_2" | tedf2$Baseline.Fibrosis.Group == "F_3_4")
      & tedf2$Timing == "6hrs post IFN",]
tedf2_noncirr_6hr=tedf2[tedf2$Treatment == "6h-IFN" & 
      (tedf2$Baseline.Fibrosis.Group == "F_0_2" | tedf2$Baseline.Fibrosis.Group == "F_3_4")
      & tedf2$Timing == "6hrs post IFN",]


tedf2_noncirr_6h_Base_merge=merge(tedf2_noncirr_Base,tedf2_noncirr_6hr,by="Patient.ID")


temat_noncirr_Base = data.matrix(tedf2_noncirr_Base)
temat_noncirr_6hr = data.matrix(tedf2_noncirr_6hr)
tedf2_noncirr_6h_Base_ratio = temat_noncirr_6hr[,5:ncol(temat_noncirr_6hr)]-
  temat_noncirr_Base[,5:ncol(temat_noncirr_Base)]

temat_cirr_6hr = data.matrix(tedf2_cirr_6hr)
temat_cirr_Base = data.matrix(tedf2_cirr_Base)
tedf2_cirr_6h_Base_ratio = temat_cirr_6hr[,5:ncol(temat_cirr_6hr)]-
  temat_cirr_Base[,5:ncol(temat_cirr_Base)]



tedf2_cirr_noncirr_Base_merge=rbind(tedf2_noncirr_Base,tedf2_cirr_Base)

tedf2_cirr_noncirr_6h_merge=rbind(tedf2_noncirr_6h_Base_ratio,tedf2_cirr_6h_Base_ratio)
tedf2_pheno=rbind(tedf2_noncirr_6hr[,1:4],tedf2_cirr_6hr[,1:4])
tedf2_cirr_noncirr_6h_merge = cbind(tedf2_pheno,tedf2_cirr_noncirr_6h_merge)


tedf2_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group[which(tedf2_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group=="F_3_4")]<-"F_0_2"
tedf2_cirr_noncirr_6h_merge$Baseline.Fibrosis.Group[which(tedf2_cirr_noncirr_6h_merge$Baseline.Fibrosis.Group=="F_3_4")]<-"F_0_2"

BFG = factor(tedf2_cirr_noncirr_6h_merge$Baseline.Fibrosis.Group)
results = colttests(data.matrix(tedf2_cirr_noncirr_6h_merge[,5:ncol(tedf2_cirr_noncirr_6h_merge)]), 
                    BFG, tstatOnly = FALSE)
results_base = colttests(data.matrix(tedf2_cirr_noncirr_Base_merge[,5:ncol(tedf2_cirr_noncirr_Base_merge)]), 
                    BFG, tstatOnly = FALSE)


colmeans_cirr <- apply(tedf2_cirr_6h_Base_ratio, 2, mean)
colmeans_noncirr <- apply(tedf2_noncirr_6h_Base_ratio, 2, mean)
cirr_noncirr_ratio=colmeans_cirr-colmeans_noncirr
FC_cirr_noncirr= 2^cirr_noncirr_ratio
FC_cirr_noncirr = ifelse(FC_cirr_noncirr<1,-1/FC_cirr_noncirr,FC_cirr_noncirr)
cirr_noncirr=cbind(FC_cirr_noncirr,results)
sig=cirr_noncirr[(cirr_noncirr$FC_cirr_noncirr < -1.5 | 
    cirr_noncirr$FC_cirr_noncirr > 1.5) & cirr_noncirr$p.value < 0.01,]

#dim(sig)
[1] 315   4

colmeans_noncirr_base <- apply(temat_noncirr_Base[,5:ncol(temat_noncirr_Base)], 2, mean)
colmeans_cirr_base <- apply(temat_cirr_Base[,5:ncol(temat_cirr_Base)], 2, mean)
cirr_noncirr_base_ratio=colmeans_cirr_base-colmeans_noncirr_base
FC_cirr_noncirr_base= 2^cirr_noncirr_base_ratio
FC_cirr_noncirr_base = ifelse(FC_cirr_noncirr_base<1,-1/FC_cirr_noncirr_base,FC_cirr_noncirr_base)
cirr_noncirr_base=cbind(FC_cirr_noncirr_base,results_base)
sig_base=cirr_noncirr_base[(cirr_noncirr_base$FC_cirr_noncirr_base < -1.5 | 
      cirr_noncirr_base$FC_cirr_noncirr_base > 1.5) & cirr_noncirr_base$p.value < 0.01,]

# dim(sig_base)
[1] 98  4

###########Look at 2 and 4 day ##########################

tedf2_DUAL_noncirr_Base=tedf2[tedf2$Treatment == "Base" & 
    (tedf2$Baseline.Fibrosis.Group == "F_0_2" | tedf2$Baseline.Fibrosis.Group == "F_3_4")
     & (tedf2$Timing == "Week 2" | tedf2$Timing == "Week 4"),]

tedf2_DUAL_cirr_Base=tedf2[tedf2$Treatment == "Base" & 
    tedf2$Baseline.Fibrosis.Group == "F_5_6" & (tedf2$Timing == "Week 2" | 
          tedf2$Timing == "Week 4"),]

tedf2_DUAL_noncirr_wk4=tedf2[(tedf2$Treatment == "Wk2" | tedf2$Treatment == "Wk4") & 
    (tedf2$Baseline.Fibrosis.Group == "F_0_2" | tedf2$Baseline.Fibrosis.Group == "F_3_4")
    & (tedf2$Timing == "Week 2" | tedf2$Timing == "Week 4"),]

tedf2_DUAL_cirr_wk4=tedf2[(tedf2$Treatment == "Wk2" | tedf2$Treatment == "Wk4") & 
  tedf2$Baseline.Fibrosis.Group == "F_5_6" & (tedf2$Timing == "Week 2" | 
  tedf2$Timing == "Week 4"),]

temat_DUAL_noncirr_Base = data.matrix(tedf2_DUAL_noncirr_Base)
temat_DUAL_cirr_Base = data.matrix(tedf2_DUAL_cirr_Base)

temat_DUAL_noncirr_wk4 = data.matrix(tedf2_DUAL_noncirr_wk4)
temat_DUAL_cirr_wk4 = data.matrix(tedf2_DUAL_cirr_wk4)

tedf2_DUAL_noncirr_wk4_Base_ratio = temat_DUAL_noncirr_wk4[,5:ncol(temat_DUAL_noncirr_wk4)]-
  temat_DUAL_noncirr_Base[,5:ncol(temat_DUAL_noncirr_Base)]
tedf2_DUAL_cirr_wk4_Base_ratio = temat_DUAL_cirr_wk4[,5:ncol(temat_DUAL_cirr_wk4)]-
  temat_DUAL_cirr_Base[,5:ncol(temat_DUAL_cirr_Base)]


tedf2_DUAL_cirr_noncirr_Base_merge=rbind(tedf2_DUAL_noncirr_Base,tedf2_DUAL_cirr_Base)
tedf2_DUAL_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group[which(tedf2_DUAL_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group=="F_3_4")]<-"F_0_2"

tedf2_DUAL_cirr_noncirr_wk4_merge=rbind(tedf2_DUAL_noncirr_wk4_Base_ratio,
                                        tedf2_DUAL_cirr_wk4_Base_ratio)
tedf2_DUAL_pheno=rbind(tedf2_DUAL_noncirr_wk4[,1:4],tedf2_DUAL_cirr_wk4[,1:4])
tedf2_DUAL_cirr_noncirr_wk4_merge = cbind(tedf2_DUAL_pheno,tedf2_DUAL_cirr_noncirr_wk4_merge)


tedf2_DUAL_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group[which(tedf2_DUAL_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group=="F_3_4")]<-"F_0_2"
BFG = factor(tedf2_DUAL_cirr_noncirr_Base_merge$Baseline.Fibrosis.Group)
DUAL_results = colttests(data.matrix(tedf2_DUAL_cirr_noncirr_wk4_merge[,5:ncol(tedf2_DUAL_cirr_noncirr_wk4_merge)]), 
                    BFG, tstatOnly = FALSE)
DUAL_results_base = colttests(data.matrix(tedf2_DUAL_cirr_noncirr_Base_merge[,5:ncol(tedf2_DUAL_cirr_noncirr_Base_merge)]), 
                    BFG, tstatOnly = FALSE)


colmeans_DUAL_noncirr_base <- apply(temat_DUAL_noncirr_Base[,5:ncol(temat_DUAL_noncirr_Base)], 2, mean)
colmeans_DUAL_cirr_base <- apply(temat_DUAL_cirr_Base[,5:ncol(temat_DUAL_cirr_Base)], 2, mean)
DUAL_cirr_noncirr_base_ratio=colmeans_DUAL_cirr_base-colmeans_DUAL_noncirr_base
FC_DUAL_cirr_noncirr_base= 2^DUAL_cirr_noncirr_base_ratio
FC_DUAL_cirr_noncirr_base = ifelse(FC_DUAL_cirr_noncirr_base<1,-1/FC_DUAL_cirr_noncirr_base,FC_DUAL_cirr_noncirr_base)
DUAL_cirr_noncirr_base=cbind(FC_DUAL_cirr_noncirr_base,DUAL_results_base)
DUAL_sig_base=cirr_noncirr_base[(DUAL_cirr_noncirr_base$FC_DUAL_cirr_noncirr_base < -1.5 | 
            DUAL_cirr_noncirr_base$FC_DUAL_cirr_noncirr_base > 1.5) & DUAL_cirr_noncirr_base$p.value < 0.01,]
#dim(DUAL_sig_base)
#[1] 105   4

#dim(DUAL_sig_base)
#[1] 324   4

colmeans_DUAL_cirr <- apply(tedf2_DUAL_cirr_wk4_Base_ratio, 2, mean)
colmeans_DUAL_noncirr <- apply(tedf2_DUAL_noncirr_wk4_Base_ratio, 2, mean)
DUAL_cirr_noncirr_ratio=colmeans_DUAL_cirr-colmeans_DUAL_noncirr
FC_DUAL_cirr_noncirr= 2^DUAL_cirr_noncirr_ratio
FC_DUAL_cirr_noncirr = ifelse(FC_DUAL_cirr_noncirr<1,-1/FC_DUAL_cirr_noncirr,FC_DUAL_cirr_noncirr)
DUAL_cirr_noncirr=cbind(FC_DUAL_cirr_noncirr,DUAL_results)
DUAL_sig=DUAL_cirr_noncirr[(DUAL_cirr_noncirr$FC_DUAL_cirr_noncirr < -1.5 | 
   DUAL_cirr_noncirr$FC_DUAL_cirr_noncirr > 1.5) & DUAL_cirr_noncirr$p.value < 0.01,]

#dim(DUAL_sig)
#[1] 274   4

##############Responders vs Non-responders#################################


##Fisher's exact test for number of responders vs non-responders (with or without cirr):
> tab
[1] 10  4
tab=cbind(tab,c(15,3))

tab   
Cirrhosis  Non-Cirr
[1,]  10 15 #resp
[2,]   4  3 #non-responders

fisher.test(tab)
Fishers Exact Test for Count Data
data:  tab
p-value = 0.6691  # not significant
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.06105214 3.75213036
sample estimates:
odds ratio 
0.5112225 

tedf3 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","Timing","Failed"), 
               function(x) colMeans(x[,24:ncol(x)]))
unique(tedf2[1:4])



##### ttest in a loop - too long ####################################################

fc_pvalmat=matrix(nrow=ncol(tedf2_cirr_noncirr_6h_merge),ncol=4)
rownames(fc_pvalmat)=names(tedf2_cirr_noncirr_6h_merge[,5:ncol(tedf2_cirr_noncirr_6h_merge[,5])]) 

for (i in 5:ncol(tedf2_cirr_noncirr_6h_merge)){
  n=i-4
ttest_cirr_noncirr_6hr=t.test(tedf2_cirr_noncirr_6h_merge[,i] 
  ~ Baseline.Fibrosis.Group, data=tedf2_cirr_noncirr_6h_merge,var.equal=TRUE)
fc_pvalmat[n,1]= ttest_cirr_noncirr_6hr$estimate[1]
fc_pvalmat[n,2]= ttest_cirr_noncirr_6hr$estimate[2]
fc_pvalmat[n,3]= ttest_cirr_noncirr_6hr$estimate[2]-ttest_cirr_noncirr_6hr$estimate[1]
fc_pvalmat[n,4]=ttest_cirr_noncirr_6hr$p.value
}

#### Group together Cirrhosis vs Non-Cirrhosis and reanalyze ############


tedf2 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","Timing",
                       "Failed"), 
               function(x) colMeans(x[,24:ncol(x)]))

tedf2$Cirrhosis.Group[tedf2$Baseline.Fibrosis.Group=="F_0_2"]="Non.Cirr"
tedf2$Cirrhosis.Group[tedf2$Baseline.Fibrosis.Group=="F_3_4"]="Non.Cirr"
tedf2$Cirrhosis.Group[tedf2$Baseline.Fibrosis.Group=="F_5_6"]="Cirr"
Cirrhosis.Group = as.factor(tedf2$Cirrhosis.Group)


tedf2$Treatment.Group[tedf2$Timing=="Pre-IFN"]="Gt1a"
tedf2$Treatment.Group[tedf2$Timing=="6hrs post IFN"]="Gt1a"
tedf2$Treatment.Group[tedf2$Timing=="Week 2"]="Gt1b"
tedf2$Treatment.Group[tedf2$Timing=="Week 4"]="Gt1b"
Treatment.Group = as.factor(tedf2$Treatment.Group)
table(Treatment.Group,Cirrhosis.Group)

Cirrhosis.Group
Treatment.Group Cirr Non-Cirr
Gt1a   20       22
Gt1b    8       14

tedf2$Treatment.Time[tedf2$Treatment=="Wk4"]="Treated"
tedf2$Treatment.Time[tedf2$Treatment=="Wk2"]="Treated"
tedf2$Treatment.Time[tedf2$Treatment=="6h-IFN"]="Treated"
tedf2$Treatment.Time[tedf2$Treatment=="Pre-IFN"]="Treated"
tedf2$Treatment.Time[tedf2$Treatment=="Base"]="Untreated"
Treatment.Time = as.factor(tedf2$Treatment.Time)

table(Treatment.Time,Cirrhosis.Group,Treatment.Group)
, , Treatment.Group = Gt1a

Cirrhosis.Group
Treatment.Time Cirr Non-Cirr
Treated     10       11
Untreated   10       11

, , Treatment.Group = Gt1b

Cirrhosis.Group
Treatment.Time Cirr Non-Cirr
Treated      4        7
Untreated    4        7

##### 1)  First analyze Gt1a ######

dim(tedf2[tedf2$Treatment.Group=="Gt1a",])
[1]    42 53625

Gt1a=tedf2[tedf2$Treatment.Group=="Gt1a",]
Gt1a <- Gt1a[ ,c(1:5,53623:53625,6:53622)]
#Gt1a$Cirrhosis.Group[Gt1a$Cirrhosis.Group=="Non-Cirr"]="Non.Cirr"

Treatment.Time.Gt1a=as.factor(Gt1a$Treatment.Time)
Cirrhosis.Group.Gt1a=as.factor(Gt1a$Cirrhosis.Group)
table(Treatment.Time.Gt1a,Cirrhosis.Group.Gt1a)
                   Cirrhosis.Group.Gt1a
Treatment.Time.Gt1a Cirr Non-Cirr
          Treated     10       11
          Untreated   10       11

Treat <- factor(paste(Gt1a$Cirrhosis.Group,Gt1a$Treatment.Time,sep=".")) 

#Cirr.Untreated     Non-Cirr.Untreated
#Cirr.Treated       Non-Cirr.Treated   

design <- model.matrix(~0+Treat) 
colnames(design) <- levels(Treat)
Gt1a.mat=t(as.matrix(Gt1a[, 9:53625]))
corfit <- duplicateCorrelation(Gt1a.mat,design,block=Gt1a$Patient.ID)
fit <- lmFit(Gt1a.mat,design,block=Gt1a$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(DiseasedvsNormalForUntreated = Cirr.Untreated-Non.Cirr.Untreated,
              DiseasedvsNormalForTreated = Cirr.Treated-Non.Cirr.Treated,
              TreatedvsUntreatedForNormal = Non.Cirr.Treated-Non.Cirr.Untreated,
              TreatedvsUntreatedForDiseased = Cirr.Treated-Cirr.Untreated,
              levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="DiseasedvsNormalForUntreated")
finalres=topTable(fit2,sort="none",n=Inf)

results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
    counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^finalres[,1:4]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

DiseasedvsNormalForUntreated_pval=p.adjust(fit2$p.value[,1], method='BH')
DiseasedvsNormalForTreated=p.adjust(fit2$p.value[,2], method='BH')
TreatedvsUntreatedForNormal=p.adjust(fit2$p.value[,3], method='BH')
TreatedvsUntreatedForDiseased=p.adjust(fit2$p.value[,4], method='BH')
pvaladjall=cbind(DiseasedvsNormalForUntreated_pval,DiseasedvsNormalForTreated,
              TreatedvsUntreatedForNormal,TreatedvsUntreatedForDiseased)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))

dim(finalres)

######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1a.mat2=(as.matrix(Gt1a[, 9:53625]))
pca=prcomp(Gt1a.mat2)
Treat <- factor(paste(Gt1a$Cirrhosis.Group,Gt1a$Treatment.Time,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Treat)), 
       type="s",size=2)
group.v<-as.vector(Treat)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Qindiv.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

for (n in 1:4) {
group=colnames(fit2$contrasts)[n]
png(paste("Volcano_plot_",n,".png",sep=""), width=1000, height = 700) 
lod = fit2[["p.value"]][,n]
lod = -log10(lod)
mtstat = fit2[["t"]][,n]
m = fit2[["coefficients"]][,n] 
o1 = order(abs(m),decreasing=TRUE)[1:100]
o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
o=union(o1,o2)
smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
abline(h=2,v=c(-1,1)) 
dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
          columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$PROBEID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres2 = merge(finalres,annot,by.x="PROBEID",by.y="PROBEID")

colnames(finalres2)
[1] "ID"                                        "DiseasedvsNormalForUntreated_FC"          
[3] "DiseasedvsNormalForTreated_FC"             "TreatedvsUntreatedForNormal_FC"           
[5] "TreatedvsUntreatedForDiseased_FC"          "DiseasedvsNormalForUntreated_pval"        
[7] "DiseasedvsNormalForTreated_pval"           "TreatedvsUntreatedForNormal_pval"         
[9] "TreatedvsUntreatedForDiseased_pval"        "DiseasedvsNormalForUntreated_pval_adjpval"
[11] "DiseasedvsNormalForTreated_adjpval"        "TreatedvsUntreatedForNormal_adjpval"      
[13] "TreatedvsUntreatedForDiseased_adjpval"     "SYMBOL"                                   
[15] "REFSEQ"                                    "GENE"

head(finalres2[order(finalres2[,6],finalres2[,2],decreasing=FALSE),])

#DiseasedvsNormalForUntreated_pval
DvsN_Untx=finalres2[(finalres2[,6]<=0.05 & finalres2[,2] > 1.5),]
cbind(DvsN_Untx$SYMBOL,DvsN_Untx$GENE)

write.table(finalres2,file="All_results_Gt1a.txt",row.names=FALSE,sep="\t",quote=FALSE)

############  2) Merge Gt1a and Gt1b ###################################
  
  dim(tedf2[tedf2$Treatment.Group=="Gt1a" | tedf2$Treatment.Group=="Gt1b",] )
[1]    64 53625

Gt1ab=tedf2[tedf2$Treatment.Group=="Gt1a" | tedf2$Treatment.Group=="Gt1b",]
Gt1ab <- Gt1ab[ ,c(1:5,53623:53625,6:53622)]
Gt1ab[,1:8]
#Gt1a$Cirrhosis.Group[Gt1a$Cirrhosis.Group=="Non-Cirr"]="Non.Cirr"

Treatment.Time.Gt1ab=as.factor(Gt1ab$Treatment.Time)
Cirrhosis.Group.Gt1ab=as.factor(Gt1ab$Cirrhosis.Group)
table(Treatment.Time.Gt1ab,Cirrhosis.Group.Gt1ab)

Cirrhosis.Group.Gt1ab
Treatment.Time.Gt1ab Cirr Non.Cirr
Treated     14       18
Untreated   14       18

Treat <- factor(paste(Gt1ab$Cirrhosis.Group,Gt1ab$Treatment.Time,sep=".")) 

#Cirr.Untreated     Non-Cirr.Untreated
#Cirr.Treated       Non-Cirr.Treated   

design <- model.matrix(~0+Treat) 
colnames(design) <- levels(Treat)
Gt1ab.mat=t(as.matrix(Gt1ab[, 9:53625]))
corfit <- duplicateCorrelation(Gt1ab.mat,design,block=Gt1ab$Patient.ID)
fit <- lmFit(Gt1ab.mat,design,block=Gt1ab$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(DiseasedvsNormalForUntreated = Cirr.Untreated-Non.Cirr.Untreated,
                    DiseasedvsNormalForTreated = Cirr.Treated-Non.Cirr.Treated,
                    TreatedvsUntreatedForNormal = Non.Cirr.Treated-Non.Cirr.Untreated,
                    TreatedvsUntreatedForDiseased = Cirr.Treated-Cirr.Untreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="DiseasedvsNormalForUntreated")
finalres=topTable(fit2,sort="none",n=Inf)

results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^finalres[,1:4]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

DiseasedvsNormalForUntreated_pval=p.adjust(fit2$p.value[,1], method='BH')
DiseasedvsNormalForTreated=p.adjust(fit2$p.value[,2], method='BH')
TreatedvsUntreatedForNormal=p.adjust(fit2$p.value[,3], method='BH')
TreatedvsUntreatedForDiseased=p.adjust(fit2$p.value[,4], method='BH')
pvaladjall=cbind(DiseasedvsNormalForUntreated_pval,DiseasedvsNormalForTreated,
                 TreatedvsUntreatedForNormal,TreatedvsUntreatedForDiseased)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))

dim(finalres)

######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1a.mat2=(as.matrix(Gt1a[, 9:53625]))
pca=prcomp(Gt1a.mat2)
Treat <- factor(paste(Gt1a$Cirrhosis.Group,Gt1a$Treatment.Time,sep=".")) 

Gt1ab.mat2=(as.matrix(Gt1ab[, 9:53625]))
pca=prcomp(Gt1ab.mat2)
Treat <- factor(paste(Gt1ab$Cirrhosis.Group,Gt1ab$Treatment.Time,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Treat)), 
       type="s",size=2)
group.v<-as.vector(Treat)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Gt1ab.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

for (n in 1:4) {
  group=colnames(fit2$contrasts)[n]
  png(paste("Volcano_plot_",n,".png",sep=""), width=1000, height = 700) 
  lod = fit2[["p.value"]][,n]
  lod = -log10(lod)
  mtstat = fit2[["t"]][,n]
  m = fit2[["coefficients"]][,n] 
  o1 = order(abs(m),decreasing=TRUE)[1:100]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
  o=union(o1,o2)
  smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
  points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$ID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres3 = merge(finalres,annot,by.x="ID",by.y="PROBEID") #Gt1ab

colnames(finalres3)
[1] "ID"                                        "DiseasedvsNormalForUntreated_FC"          
[3] "DiseasedvsNormalForTreated_FC"             "TreatedvsUntreatedForNormal_FC"           
[5] "TreatedvsUntreatedForDiseased_FC"          "DiseasedvsNormalForUntreated_pval"        
[7] "DiseasedvsNormalForTreated_pval"           "TreatedvsUntreatedForNormal_pval"         
[9] "TreatedvsUntreatedForDiseased_pval"        "DiseasedvsNormalForUntreated_pval_adjpval"
[11] "DiseasedvsNormalForTreated_adjpval"        "TreatedvsUntreatedForNormal_adjpval"      
[13] "TreatedvsUntreatedForDiseased_adjpval"     "SYMBOL"                                   
[15] "REFSEQ"                                    "GENE"

head(finalres2[order(finalres2[,6],finalres2[,2],decreasing=FALSE),])
head(finalres2[order(finalres3[,6],finalres3[,2],decreasing=FALSE),])



#DiseasedvsNormalForUntreated_pval
DvsN_Untx=finalres2[(finalres2[,6]<=0.05 & finalres2[,2] > 1.5),]
cbind(DvsN_Untx$SYMBOL,DvsN_Untx$GENE)

write.table(finalres3,file="All_results_Gt1ab.txt",row.names=FALSE,sep="\t",quote=FALSE)

###### 3) Look at Failed vs Success ################

Gt1ab$Failed[Gt1ab$Failed=="Withdrew"]="1"

Gt1ab$Failure[Gt1ab$Failed=="1"]="Failed"
Gt1ab$Failure[Gt1ab$Failed=="0"]="Success"

Failure.Gt1ab=as.factor(Gt1ab$Failure)
table(Treatment.Time.Gt1ab,Failure.Gt1ab)

Failure.Gt1ab
Treatment.Time.Gt1ab Failed Success
Treated        7      25
Untreated      7      25

Failure <- factor(paste(Gt1ab$Failure,Gt1ab$Treatment.Time,sep=".")) 
unique(Failure)

[1] Success.Untreated Success.Treated   Failed.Untreated  Failed.Treated   
Levels: Failed.Treated Failed.Untreated Success.Treated Success.Untreated

design <- model.matrix(~0+Failure) 
colnames(design) <- levels(Failure)
Gt1ab.mat=t(as.matrix(Gt1ab[, 9:53625]))
corfit <- duplicateCorrelation(Gt1ab.mat,design,block=Gt1ab$Patient.ID)
fit <- lmFit(Gt1ab.mat,design,block=Gt1ab$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(FailedvsSuccessUntreated = Failed.Untreated-Success.Untreated,
                    FailedvsSuccessForTreated = Failed.Treated-Success.Treated,
                    TreatedvsUntreatedForSuccess = Success.Treated-Success.Untreated,
                    TreatedvsUntreatedForFailed = Failed.Treated-Failed.Untreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="FailedvsSuccessUntreated")
finalres=topTable(fit2,sort="none",n=Inf)

results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^finalres[,1:4]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

DiseasedvsNormalForUntreated_pval=p.adjust(fit2$p.value[,1], method='BH')
DiseasedvsNormalForTreated=p.adjust(fit2$p.value[,2], method='BH')
TreatedvsUntreatedForNormal=p.adjust(fit2$p.value[,3], method='BH')
TreatedvsUntreatedForDiseased=p.adjust(fit2$p.value[,4], method='BH')
pvaladjall=cbind(DiseasedvsNormalForUntreated_pval,DiseasedvsNormalForTreated,
                 TreatedvsUntreatedForNormal,TreatedvsUntreatedForDiseased)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))

dim(finalres)

######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1a.mat2=(as.matrix(Gt1a[, 9:53625]))
pca=prcomp(Gt1a.mat2)
Treat <- factor(paste(Gt1a$Cirrhosis.Group,Gt1a$Treatment.Time,sep=".")) 

Gt1ab.mat2=(as.matrix(Gt1ab[, 9:53625]))
pca=prcomp(Gt1ab.mat2)
Treat <- factor(paste(Gt1ab$Failure,Gt1ab$Treatment.Time,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Treat)), 
       type="s",size=2)
group.v<-as.vector(Treat)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Gt1ab.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

for (n in 1:4) {
  group=colnames(fit2$contrasts)[n]
  png(paste("Volcano_plot_",n,".png",sep=""), width=1000, height = 700) 
  lod = fit2[["p.value"]][,n]
  lod = -log10(lod)
  mtstat = fit2[["t"]][,n]
  m = fit2[["coefficients"]][,n] 
  o1 = order(abs(m),decreasing=TRUE)[1:100]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
  o=union(o1,o2)
  smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
  points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$ID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres3 = merge(finalres,annot,by.x="ID",by.y="PROBEID") #Gt1ab

colnames(finalres3)
[1] "ID"                                        "DiseasedvsNormalForUntreated_FC"          
[3] "DiseasedvsNormalForTreated_FC"             "TreatedvsUntreatedForNormal_FC"           
[5] "TreatedvsUntreatedForDiseased_FC"          "DiseasedvsNormalForUntreated_pval"        
[7] "DiseasedvsNormalForTreated_pval"           "TreatedvsUntreatedForNormal_pval"         
[9] "TreatedvsUntreatedForDiseased_pval"        "DiseasedvsNormalForUntreated_pval_adjpval"
[11] "DiseasedvsNormalForTreated_adjpval"        "TreatedvsUntreatedForNormal_adjpval"      
[13] "TreatedvsUntreatedForDiseased_adjpval"     "SYMBOL"                                   
[15] "REFSEQ"                                    "GENE"

head(finalres2[order(finalres2[,6],finalres2[,2],decreasing=FALSE),])
head(finalres2[order(finalres3[,6],finalres3[,2],decreasing=FALSE),])

############  4) Look at Gt1b only ###################################


dim(tedf2[tedf2$Treatment.Group=="Gt1b",] )
[1]    22 53625
dim(Gt1b)
Gt1b=tedf2[tedf2$Treatment.Group=="Gt1b",]
Gt1b <- Gt1b[ ,c(1:5,53623:53625,6:53622)]
Gt1b[,1:8]
Gt1b <- Gt1b[Gt1b$Patient.ID!="772",] # Remove Withdrawn Patient

#Gt1a$Cirrhosis.Group[Gt1a$Cirrhosis.Group=="Non-Cirr"]="Non.Cirr"

Treatment=as.factor(Gt1b$Treatment.Time)
#Cirrhosis.Group.Gt1b=as.factor(Gt1b$Cirrhosis.Group)
#table(Treatment.Time.Gt1b,Cirrhosis.Group.Gt1b)
Patient <- factor(Gt1b$Patient.ID) 

table(Patient,Treatment)
#Treatment
#Patient Treated Untreated
110        1         1
513        1         1
814        1         1
975        1         1
1059       1         1
1154       1         1
1230       1         1
1505       1         1
2168       1         1
2783       1         1
2957       1         1

#design <- model.matrix(~Patient+Treatment) 
Gt1b.mat=t(as.matrix(Gt1b[, 9:53625]))

#######PCA Plot #################

colnames(Gt1b.mat)=paste(Gt1b$Patient,Gt1b$Treatment.Time,Gt1b$Failed,sep=".")
Treatment <- relevel(Treatment, ref="Untreated")
design <- model.matrix(~0+Treatment) 

corfit <- duplicateCorrelation(Gt1b.mat,design,block=Gt1b$Patient.ID)
corfit$consensus
fit <- lmFit(Gt1b.mat,design,block=Gt1b$Patient.ID,correlation=corfit$consensus)

#fit <- lmFit(Gt1b, design) 
#fit <- eBayes(fit) 
#topTable(fit, coef="TreatmentUntreated")  
#finalres=topTable(fit,coef="TreatmentUntreated",sort="none",n=Inf)

cm <- makeContrasts(Treatment = TreatmentTreated-TreatmentUntreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1)
volcanoplot(fit2,coef=1)
finalres=topTable(fit2,sort="none",n=Inf)

FC = 2^(fit2$coefficients)
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(fit2$coefficient),"FC",sep="_")
pval=fit2$p.value
colnames(pval)=paste(colnames(fit2$p.value),"pval",sep="_")

finalres = as.data.frame(cbind(Gt1b.mat,FC,pval,finalres))

dim(finalres)
head(finalres)

#key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$ID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres5 = merge(finalres,annot,by.x="ID",by.y="PROBEID") #Gt1b
write.table(finalres5,file="all_results_Gt1b.txt",row.names=FALSE,sep="\t",quote=FALSE)

####Heatmap of all genes differentially expressed ###########
#Global heat map of all ISGs (2000 of them? Though 469 were differentially expressed)

head(finalres5)
finalres5_filt=filter(finalres5,abs(Treatment_FC) >= 1.2,Treatment_pval<=0.01,SYMBOL != "NA")
finalres5_filt=arrange(finalres5_filt,Treatment_pval)
finalres5_filt=finalres5_filt[!duplicated(finalres5_filt$SYMBOL), ]
dim(finalres5_filt)
[1] 386  12
head(finalres5_filt)
colnames(finalres5_filt)
finalres_heat=finalres5_filt[,c(2:23)]
colnames(finalres_heat)=paste(Gt1b$Treatment.Time,Gt1b$Failed,sep=".")
rownames(finalres_heat)=finalres5_filt$SYMBOL

palette <- colorRampPalette(c('blue','red'))(12)
palette2 <- c('green','yellow','red','blue')

palette
color=as.numeric(as.factor(colnames(finalres_heat)))
f.patient=as.factor(colnames(finalres_heat))

hc.rows <- hclust(dist(finalres_heat))

heatmap(as.matrix(finalres_heat), 
        ColSideColors = palette2[color],
        col=palette, labCol= NA, labRow=NA,
        Rowv=as.dendrogram(hc.rows),
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(finalres_heat)),
              method="average")))

legend("topright",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.8)

##########  Look at Treatment Failure within Gt1b ############


dim(Gt1b)
Gt1b[1:10,1:10]
#Gt1b <- Gt1b[ ,c(1:5,53623:53625,6:53622)]
Gt1b$Failed[Gt1b$Failed=="1"]="Failed"
Gt1b$Failed[Gt1b$Failed=="0"]="Success"
Gt1b_Tx=Gt1b[Gt1b$Treatment.Time=="Treated",]

Patient <- factor(Gt1b$Patient.ID) 
Treatment=as.factor(Gt1b$Treatment.Time)
Failure=as.factor(Gt1b$Failed)

table(Failure,Treatment)
Failure
Treatment   Failed Success
Treated        4       7
Untreated      4       7

Failure <- factor(paste(Gt1b$Failed,Gt1b$Treatment.Time,sep=".")) 
unique(Failure)

# [1] Success.Untreated Success.Treated   Failed.Untreated  Failed.Treated   
# Levels: Failed.Treated Failed.Untreated Success.Treated Success.Untreated

design <- model.matrix(~0+Failure) 
colnames(design) <- levels(Failure)
Gt1b.mat=t(as.matrix(Gt1b[, 9:53625]))
corfit <- duplicateCorrelation(Gt1b.mat,design,block=Gt1b$Patient.ID)
corfit$consensus
fit <- lmFit(Gt1b.mat,design,block=Gt1b$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(FailedvsSuccessUntreated = Failed.Untreated-Success.Untreated,
                    FailedvsSuccessForTreated = Failed.Treated-Success.Treated,
                    TreatedvsUntreatedForSuccess = Success.Treated-Success.Untreated,
                    TreatedvsUntreatedForFailed = Failed.Treated-Failed.Untreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="FailedvsSuccessUntreated")
volcanoplot(fit2,coef=1)
finalres=topTable(fit2,sort="none",n=Inf)



results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^(fit2$coefficients[,1:4])
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

FailedvsSuccessUntreated_pval=p.adjust(fit2$p.value[,1], method='BH')
FailedvsSuccessForTreated=p.adjust(fit2$p.value[,2], method='BH')
TreatedvsUntreatedForSuccess=p.adjust(fit2$p.value[,3], method='BH')
TreatedvsUntreatedForFailed=p.adjust(fit2$p.value[,4], method='BH')
pvaladjall=cbind(FailedvsSuccessUntreated_pval,FailedvsSuccessForTreated,
                 TreatedvsUntreatedForSuccess,TreatedvsUntreatedForFailed)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

colnames(Gt1b.mat)=paste(Failure,Gt1b$Patient.ID,sep=".")
finalres = as.data.frame(cbind(Gt1b.mat,FC,pvalall,pvaladjall))

FailedvsSuccessUntx=finalres[(finalres[,1]>=1.2 | finalres[,1]<=-1.2) && finalres[,3] < 0.1]

TxvsUntx_Success=finalres[((finalres$TreatedvsUntreatedForSuccess_FC>=1.2 | 
                             finalres$TreatedvsUntreatedForSuccess_FC<=-1.2) &
                             finalres$TreatedvsUntreatedForSuccess_pval < 0.1),]


TxvsUntx_Failed=finalres[((finalres[,4]>=1.2 | finalres[,4]<=-1.2) & finalres[,8] < 0.1),]

dim(finalres)
[1] 53617    34


######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1b.mat2=(as.matrix(Gt1b[, 9:53625]))
pca=prcomp(Gt1b.mat2)
#Treat <- factor(paste(Gt1b$Cirrhosis.Group,Gt1b$Treatment.Time,sep=".")) 
plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
pca$x[,1:3] 
plot(pca,type="lines") 
cell_rep=paste(Patient,Treatment,Failure,sep=".")
#cell_rep=paste(Phenotype,Failure,Treatment,sep=".")
library(rgl)
open3d() 
Gt1b.mat2$group = as.factor(Patient)
plot3d(pca$x[,1:3],col=as.integer(Gt1b.mat2$group), 
       type="s",size=2)
group.v<-as.vector(cell_rep)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1.5) 
rgl.postscript("pca3d_Gt1ab.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

n=12
  group=colnames(fit$contrasts)[n]
  png(paste("Volcano_plot_",n,".png",sep=""), width=1000, height = 700) 
  lod = fit[["p.value"]][,n]
  lod = -log10(lod)
  mtstat = fit[["t"]][,n]
  m = fit[["coefficients"]][,n] 
  o1 = order(abs(m),decreasing=TRUE)[1:100]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
  o=union(o1,o2)
  smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
  points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$ID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres6 = merge(finalres,annot,by.x="ID",by.y="PROBEID") #Gt1b
write.table(finalres6,file="failvsnonfail_results_Gt1b.txt",row.names=FALSE,sep="\t",quote=FALSE)

###### Look at ISG Expression in Gt1b ######

ISG=c("ACADSB","ACSM2A","ACSM3","ADH6","ADIPOQ","ANXA9","APOA5",
      "APOF","BDH1","CCL19","CCL2","CCL20","CCL4","CRIP1","CXCL10",
      "CXCL11","CXCL9","CYP1A1","CYP4A11","DDX58","DDX60","EIF2AK2",
      "FCGR1A","FDPS","HCP5","HERC5","HMGCS2","IFI27","IFI35","IFI44L",
      "IFI6","IFIH1","IFIT1","IFIT2","IFIT3","IFITM1","IRF9","ISG15",
      "ISG20","LCAT","LEP","LGALS3","LIPC","LSS","MUC13","MX1","MX2",
      "NPC1L1","OAS1","OAS2","OAS3","OASL","PLSCR1","PPARA","PRKAB2",
      "RSAD2","SLFN5","STAT1","TAP1","UBE2L6","XAF1")

ISG46=c("CXCL10","CXCL9","CXCL11","CCL4","CCL2","CCL19","CCL20","TAP1",
        "OASL","ISG15","MX1","OAS3","OAS2","IFIT1","IFIT3","STAT1","OAS1",
        "IFIT2","IFITM1","ISG20","BDH1","PRKAB2","CYP4A11","HMGCS2","PPARA",
        "APOA5","CYP1A1","ACADSB","ADH6","ANXA9","ACSM3","LIPC","NPC1L1","FDPS",
        "LSS","LCAT","APOF","DDX60","PLSCR1","EIF2AK2","IRF9","IFI6","IFI44L",
        "IFI27","IFIH1","DDX58","RSAD2")

#ISG=c("FCGR1B","ADIPOQ","RSAD2")
#ISG=c("CXCL10","CXCL11","RSAD2")

ISG_annotsub=subset(annot,annot$SYMBOL %in% ISG46 == TRUE)
ISG_sub=Gt1b[,match(ISG_annotsub$PROBEID,colnames(Gt1b))]

ISG_sub=cbind(Gt1b[,1:8],ISG_sub)
colnames(ISG_sub)=ISG_annotsub$SYMBOL[match(colnames(ISG_sub),ISG_annotsub$PROBEID)]
names(ISG_sub)[1:8] <- names(Gt1b)[1:8]
dim(ISG_sub)
[1] 22 62
numcol=dim(ISG_sub)[2]

FC.df=data.frame(matrix(ncol = numcol, nrow = 11))
colnames(FC.df)=colnames(ISG_sub)

for (i in 0:10){
  j=i*2+1
  k=j+1
  l=i+1
  FC.df[l,1:8]=ISG_sub[k,1:8]
  FC.df[l,9:numcol]=ISG_sub[k,9:numcol]-ISG_sub[j,9:numcol]
}

library(ggplot2)
# Use single fill color

length = dim(ISG_sub)[2]-8
for(i in 1:length){
  j=i+8
  filename=paste(names(ISG_sub)[j],"_",i,sep="")
  #png(paste(filename,"png",sep="."))
  ymin=floor(min(ISG_sub[,j]))
  ymax=ceiling(max(ISG_sub[,j]))
  size=(ymax-ymin)*0.05
plots=ggplot(ISG_sub, aes(x=Treatment, y=ISG_sub[,j], fill=Failed)) + 
  geom_dotplot(binaxis='y', stackdir='center',binwidth = size) +
  ylab(names(ISG_sub)[j]) + ylim(c(ymin,ymax)) + theme_bw()  
ggsave(plots,filename=paste(filename,"png",sep="."))
dev.off()
}

length = dim(FC.df)[2]-8
for(i in 1:length){
  j=i+8
  filename=paste(names(FC.df)[j],"_FC_",i,sep="")
  #png(paste(filename,"png",sep="."))
  ymin=floor(min(FC.df[,j]))
  ymax=ceiling(max(FC.df[,j]))
  size=(ymax-ymin)*0.05
  plots=ggplot(FC.df, aes(x=Failed, y=FC.df[,j], fill=Failed)) + 
    geom_boxplot(notch = FALSE) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = size) +
    ylab(names(FC.df)[j]) + ylim(c(ymin,ymax)) + theme_bw()  
  ggsave(plots,filename=paste(filename,"png",sep="."))
  dev.off()
}

save.image()

#########Draw heatmap #######
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#limit to highest gene signal 
tISG_sub=t(ISG_sub)

colnames(tISG_sub)=paste(tISG_sub["Failed",],tISG_sub["Treatment.Time",],
                         sep=".")
colnames(tISG_sub)

#Filter more
myvars <- colnames(tISG_sub) %in% c("Success.Treated","Failed.Treated")
colnames(tISG_sub[,!myvars])
tISG_sub <- tISG_sub[,!myvars] 

dim(tISG_sub)
[1] 65 22

tISG_sub[1:14,]
tISG_sub.mat=as.matrix(tISG_sub[9:dim(tISG_sub)[1],])
rownames(tISG_sub.mat)=rownames(tISG_sub[9:dim(tISG_sub)[1],])
class(tISG_sub.mat) <- "numeric"
tISG_sub.mat <- cbind(tISG_sub.mat,apply(tISG_sub.mat,1,mean))
dim(tISG_sub.mat)
[1] 73 23

dimnames(tISG_sub.mat)[[2]][12]<- "Avg"
dimnames(tISG_sub.mat)[[2]][23]<- "Avg" # For unfiltered
head(tISG_sub.mat)
tISG_sub.mat.ordered=tISG_sub.mat[order(tISG_sub.mat[,23]),] 
tISG_sub.mat.unique=unique(tISG_sub.mat.ordered,fromLast=TRUE)
tISG_sub.mat.unique.scale=tISG_sub.mat.unique[,1:22]
hc.rows <- hclust(dist(tISG_sub.mat.unique.scale))

save.image()
palette <- colorRampPalette(c('blue','red'))(12)
palette2 <- c('green','orange','red','blue')

palette
color=as.numeric(as.factor(colnames(tISG_sub.mat.unique.scale)))
f.patient=as.factor(colnames(tISG_sub.mat.unique.scale))

heatmap(tISG_sub.mat.unique.scale, 
    ColSideColors = palette2[color],
    col=palette, labCol= NA,
    Rowv=as.dendrogram(hc.rows),
    #Colv = NA)
    Colv=as.dendrogram(hclust(dist(t(tISG_sub.mat.unique[,1:22]),)),
                       method="average"))

legend("topright",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.8)

#######Heatmap of all genes ########

#limit to highest gene signal

colnames(Gt1b)=annot$SYMBOL[match(colnames(Gt1b),annot$PROBEID)]
Gt1b_sub=Gt1b[(colnames(Gt1b) != "NA")==TRUE,]
dim(Gt1b_sub)
tGt1b_sub=t(Gt1b_sub)
colnames(tISG_sub)=paste(tISG_sub["Failed",],tISG_sub["Treatment.Time",],
                         sep=".")
tISG_sub.mat=as.matrix(tISG_sub[9:dim(tISG_sub)[1],])
rownames(tISG_sub.mat)=rownames(tISG_sub[9:dim(tISG_sub)[1],])
class(tISG_sub.mat) <- "numeric"
tISG_sub.mat <- cbind(tISG_sub.mat,apply(tISG_sub.mat,1,mean))
dim(tISG_sub.mat)
[1] 73 23
dimnames(tISG_sub.mat)[[2]][23]<- "Avg"
head(tISG_sub.mat)
tISG_sub.mat.ordered=tISG_sub.mat[order(tISG_sub.mat[,23]),] 
tISG_sub.mat.unique=unique(tISG_sub.mat.ordered,fromLast=TRUE)
tISG_sub.mat.unique.scale=tISG_sub.mat.unique[,1:22]
hc.rows <- hclust(dist(tISG_sub.mat.unique.scale))

save.image()
palette <- colorRampPalette(c('red','blue'))(12)
palette2 <- c('red','blue','green','yellow')

palette
color=as.numeric(as.factor(colnames(tISG_sub.mat.unique.scale)))
f.patient=as.factor(colnames(tISG_sub.mat.unique.scale))

heatmap(tISG_sub.mat.unique.scale, 
        ColSideColors = palette2[color],
        col=palette, labCol= NA,
        Rowv=as.dendrogram(hc.rows),
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(tISG_sub.mat.unique[,1:22])))))

legend("topright",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.7)


################Re-do finalres6 - just treated vs untreated################

colnames(finalres5)
[1] "ID"        "FC"        "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val"
[8] "B"         "PROBEID"   "SYMBOL"    "REFSEQ"    "GENE"  

write.table(finalres5,file="all_results_Gt1b.txt",row.names=FALSE,sep="\t",quote=FALSE)


#########################
colnames(finalres6)
[1] "ID"                                    "FailedvsSuccessUntreated_FC"          
[3] "FailedvsSuccessForTreated_FC"          "TreatedvsUntreatedForSuccess_FC"      
[5] "TreatedvsUntreatedForFailed_FC"        "FailedvsSuccessUntreated_pval"        
[7] "FailedvsSuccessForTreated_pval"        "TreatedvsUntreatedForSuccess_pval"    
[9] "TreatedvsUntreatedForFailed_pval"      "FailedvsSuccessUntreated_pval_adjpval"
[11] "FailedvsSuccessForTreated_adjpval"     "TreatedvsUntreatedForSuccess_adjpval" 
[13] "TreatedvsUntreatedForFailed_adjpval"   "SYMBOL"                               
[15] "REFSEQ"                                    "GENE"

head(finalres6[order(finalres2[,6],finalres2[,2],decreasing=FALSE),])
write.table(finalres4.GO,file="Gt1b_results.txt",row.names=FALSE,sep="\t",quote=FALSE)

#DiseasedvsNormalForUntreated_pval
DvsN_Untx=finalres2[(finalres2[,6]<=0.05 & finalres2[,2] > 1.5),]
cbind(DvsN_Untx$SYMBOL,DvsN_Untx$GENE)

write.table(finalres6,file="Fail_success_results_Gt1b.txt",row.names=FALSE,sep="\t",quote=FALSE)

finalres4.GO=cbind(finalres4$ID,finalres4$FC,finalres4$P.Value)
colnames(finalres4.GO)=c("ID","FC","pval")
write.table(finalres4.GO,file="Gt1b_results.txt",row.names=FALSE,sep="\t",quote=FALSE)

####### 5) Look at Gt1a Pre vs Post Analysis ######################

dim(tedf2[tedf2$Treatment.Group=="Gt1a",])
[1]    42 53625

Gt1a=tedf2[tedf2$Treatment.Group=="Gt1a",]
Gt1a <- Gt1a[ ,c(1:5,53623:53625,6:53622)]
#Gt1a$Cirrhosis.Group[Gt1a$Cirrhosis.Group=="Non-Cirr"]="Non.Cirr"

Gt1a$Timing[Gt1a$Timing=="6h_post_IFN"]="Post_IFN"
Gt1a$Timing[Gt1a$Timing=="Pre-IFN"]="Pre_IFN"
Gt1a$Treatment[Gt1a$Treatment=="Pre-IFN"]="Pre_IFN"
Gt1a$Treatment[Gt1a$Treatment=="6h-IFN"]="Post_IFN"

IFN.Time.Gt1a = as.factor(Gt1a$Treatment)
Treatment.Timing.Gt1a=as.factor(Gt1a$Timing)

table(Treatment.Timing.Gt1a,IFN.Time.Gt1a)
IFN.Time.Gt1a
Treatment.Timing.Gt1a Base Post_IFN Pre_IFN
Post_IFN   11       11       0
Pre_IFN    10        0      10

Treat <- factor(paste(Gt1a$Treatment,Gt1a$Timing,sep=".")) 
unique(Treat)

[1] Base.Pre_IFN      Pre_IFN.Pre_IFN   Post_IFN.Post_IFN Base.Post_IFN    
Levels: Base.Post_IFN Base.Pre_IFN Post_IFN.Post_IFN Pre_IFN.Pre_IFN

design <- model.matrix(~0+Treat) 
colnames(design) <- levels(Treat)
Gt1a.mat=t(as.matrix(Gt1a[, 9:53625]))
corfit <- duplicateCorrelation(Gt1a.mat,design,block=Gt1a$Patient.ID)
fit <- lmFit(Gt1a.mat,design,block=Gt1a$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(IFN_6hrs_vs_base = Post_IFN.Post_IFN-Base.Post_IFN,
                    Pre_IFN_vs_base = Pre_IFN.Pre_IFN-Base.Pre_IFN,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="IFN_6hrs_vs_base")
finalres=topTable(fit2,sort="none",n=Inf)

results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^finalres[,1:2]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

IFN_6hrs_vs_base_pval=p.adjust(fit2$p.value[,1], method='BH')
Pre_IFN_vs_base_pval=p.adjust(fit2$p.value[,2], method='BH')

pvaladjall=cbind(IFN_6hrs_vs_base_pval,Pre_IFN_vs_base_pval)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))
head(finalres)
dim(finalres)

######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1a.mat2=(as.matrix(Gt1a[, 9:53625]))
pca=prcomp(Gt1a.mat2)
Treat <- factor(paste(Gt1a$Treatment,Gt1a$Timing,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Treat)), 
       type="s",size=2)
group.v<-as.vector(Treat)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Gt1a_TreatTiming.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

for (n in 1:2) {
  group=colnames(fit2$contrasts)[n]
  png(paste("Volcano_plot_Gt1a_TT",n,".png",sep=""), width=1000, height = 700) 
  lod = fit2[["p.value"]][,n]
  lod = -log10(lod)
  mtstat = fit2[["t"]][,n]
  m = fit2[["coefficients"]][,n] 
  o1 = order(abs(m),decreasing=TRUE)[1:100]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
  o=union(o1,o2)
  smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
  points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$PROBEID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres5 = merge(finalres,annot,by.x="PROBEID",by.y="PROBEID")

colnames(finalres5)
[1] "PROBEID"                       "IFN_6hrs_vs_base_FC"          
[3] "Pre_IFN_vs_base_FC"            "IFN_6hrs_vs_base_pval"        
[5] "Pre_IFN_vs_base_pval"          "IFN_6hrs_vs_base_pval_adjpval"
[7] "Pre_IFN_vs_base_pval_adjpval"  "SYMBOL"                       
[9] "REFSEQ"                                   "GENE"

head(finalres5[order(finalres2[,4],finalres2[,2],decreasing=FALSE),])

#DiseasedvsNormalForUntreated_pval
Post_IFN_Results=finalres5[,c(1,2,4)]
Pre_IFN_Results=finalres5[,c(1,3,5)]

write.table(finalres5,file="Gt1a_all_TT.txt",row.names=FALSE,sep="\t",quote=FALSE)
write.table(Pre_IFN_Results,file="Gt1a_pre_TT_Results.txt",row.names=FALSE,sep="\t",quote=FALSE)
write.table(Post_IFN_Results,file="Gt1a_post_TT_Results.txt",row.names=FALSE,sep="\t",quote=FALSE)

####### 5) Look at Gt1a Pre-Base vs Post-Base Analysis ######################

dim(tedf2[tedf2$Treatment.Group=="Gt1a",])
[1]    42 53625

Gt1a=tedf2[tedf2$Treatment.Group=="Gt1a",]
Gt1a <- Gt1a[ ,c(1:5,53623:53625,6:53622)]
#Gt1a$Cirrhosis.Group[Gt1a$Cirrhosis.Group=="Non-Cirr"]="Non.Cirr"

Gt1a$Timing[Gt1a$Timing=="6h_post_IFN"]="Post_IFN"
Gt1a$Timing[Gt1a$Timing=="Pre-IFN"]="Pre_IFN"
Gt1a$Treatment[Gt1a$Treatment=="Pre-IFN"]="Pre_IFN"
Gt1a$Treatment[Gt1a$Treatment=="6h-IFN"]="Post_IFN"

IFN.Time.Gt1a = as.factor(Gt1a$Treatment)
Treatment.Timing.Gt1a=as.factor(Gt1a$Timing)

Pre_Group = Gt1a[Gt1a$Timing=="Pre_IFN",]
Post_Group = Gt1a[Gt1a$Timing=="Post_IFN",]

Treat <- factor(paste(Gt1a$Treatment,Gt1a$Timing,sep=".")) 
unique(Treat)

[1] Base.Pre_IFN      Pre_IFN.Pre_IFN   Post_IFN.Post_IFN Base.Post_IFN    
Levels: Base.Post_IFN Base.Pre_IFN Post_IFN.Post_IFN Pre_IFN.Pre_IFN

design <- model.matrix(~0+Treat) 
colnames(design) <- levels(Treat)
Gt1a.mat=t(as.matrix(Gt1a[, 9:53625]))
corfit <- duplicateCorrelation(Gt1a.mat,design,block=Gt1a$Patient.ID)
fit <- lmFit(Gt1a.mat,design,block=Gt1a$Patient.ID,correlation=corfit$consensus)

cm <- makeContrasts(IFN_6hrs_vs_base = Post_IFN.Post_IFN-Base.Post_IFN,
                    Pre_IFN_vs_base = Pre_IFN.Pre_IFN-Base.Pre_IFN,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="IFN_6hrs_vs_base")
finalres=topTable(fit2,sort="none",n=Inf)

results <- decideTests(fit2,adjust.method="none",p.value=0.05,lfc=0.58)
vennDiagram(results[,1:2],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)
vennDiagram(results[,3:4],include=c("up","down"),
            counts.col=c("red","green"),cex=1.0)

FC = 2^finalres[,1:2]
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

IFN_6hrs_vs_base_pval=p.adjust(fit2$p.value[,1], method='BH')
Pre_IFN_vs_base_pval=p.adjust(fit2$p.value[,2], method='BH')

pvaladjall=cbind(IFN_6hrs_vs_base_pval,Pre_IFN_vs_base_pval)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))
head(finalres)
dim(finalres)

######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1a.mat2=(as.matrix(Gt1a[, 9:53625]))
pca=prcomp(Gt1a.mat2)
Treat <- factor(paste(Gt1a$Treatment,Gt1a$Timing,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Treat)), 
       type="s",size=2)
group.v<-as.vector(Treat)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("pca3d_Gt1a_TreatTiming.pdf","pdf")

####### Volcano Plot ################
#draw volcano plot using moderated t-statistics

for (n in 1:2) {
  group=colnames(fit2$contrasts)[n]
  png(paste("Volcano_plot_Gt1a_TT",n,".png",sep=""), width=1000, height = 700) 
  lod = fit2[["p.value"]][,n]
  lod = -log10(lod)
  mtstat = fit2[["t"]][,n]
  m = fit2[["coefficients"]][,n] 
  o1 = order(abs(m),decreasing=TRUE)[1:100]
  o2 = order(abs(mtstat),decreasing=TRUE)[1:100]
  o=union(o1,o2)
  smoothScatter(m,lod,main=group, xlab="Log-ratio",ylab="LOD")
  points(m[o1],lod[o1],pch=18,col="blue")  # most significant fold change (abs)
  points(m[o2],lod[o2],pch=1,col="red")  # most significant p-val
  abline(h=2,v=c(-1,1)) 
  dev.off()
}

###### Annotate ##########

#attr(affyData,"annotation")
[1] "pd.hugene.2.0.st"

library("pd.hugene.2.0.st")
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
ls("package:hugene20stprobeset.db") #probe level
ls("package:hugene20sttranscriptcluster.db")

key=keys(hugene20stprobeset.db) #probe level
key=keys(hugene20sttranscriptcluster.db) #Probeset level
annot = AnnotationDbi::select(hugene20sttranscriptcluster.db, keys=key, 
                              columns=c("SYMBOL","REFSEQ","GENENAME"),keytype="PROBEID")
annot=annot[!duplicated(annot$PROBEID), ]

finalres$PROBEID = rownames(finalres)
colnames(annot) = c("PROBEID","SYMBOL","REFSEQ","GENE")

finalres5 = merge(finalres,annot,by.x="PROBEID",by.y="PROBEID")

colnames(finalres5)
[1] "PROBEID"                       "IFN_6hrs_vs_base_FC"          
[3] "Pre_IFN_vs_base_FC"            "IFN_6hrs_vs_base_pval"        
[5] "Pre_IFN_vs_base_pval"          "IFN_6hrs_vs_base_pval_adjpval"
[7] "Pre_IFN_vs_base_pval_adjpval"  "SYMBOL"                       
[9] "REFSEQ"                                   "GENE"

head(finalres5[order(finalres2[,4],finalres2[,2],decreasing=FALSE),])

#DiseasedvsNormalForUntreated_pval
Post_IFN_Results=finalres5[,c(1,2,4)]
Pre_IFN_Results=finalres5[,c(1,3,5)]

write.table(finalres5,file="Gt1a_all_TT.txt",row.names=FALSE,sep="\t",quote=FALSE)
write.table(Pre_IFN_Results,file="Gt1a_pre_TT_Results.txt",row.names=FALSE,sep="\t",quote=FALSE)
write.table(Post_IFN_Results,file="Gt1a_post_TT_Results.txt",row.names=FALSE,sep="\t",quote=FALSE)

###########  6) Add Viral Load into Equation ##################

##Look at original phenotype##
tedf[1:22,1:23]

colnames(tedf[,1:30])
[1] "Patient.ID"                   "Replicate"                   
[3] "Treatment"                    "Patient"                     
[5] "Patient.ID.1"                 "Age"                         
[7] "Sex"                          "Race"                        
[9] "Baseline.Ishak.Score"         "Baseline.Fibrosis.Group"     
[11] "HCV.Genotype"                 "Timing"                      
[13] "Basline.Bx.Code"              "Date.of.Baseline.Bx"         
[15] "F.U.Bx.Code"                  "Date.of.F.U.Bx"              
[17] "Days.Base.to.F.U.Bx"          "F.U.Bx.Week"                 
[19] "HCV.RNA.at.F.U.Bx..Neg..Pos." "HCV.RNA.at.F.U.Bx.Value"     
[21] "SVR12"                        "Failed"                      
[23] "Cause.of.Failure"             "16650001"  


tedf3 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","HCV.Genotype",
                       "Failed","HCV.RNA.at.F.U.Bx.Value"), 
               function(x) colMeans(x[,24:ncol(x)]))
c("tedf3$HCV.RNA.at.F.U.Bx.Value","tedf3$Failed")


dim(tedf3[tedf3$HCV.Genotype=="1b",] )
[1]    22 53625

Gt1b=tedf3[tedf3$HCV.Genotype=="1b",]
Gt1b[,1:6]

Gt1b$Virus="1"
Gt1b$Virus[(Gt1b$HCV.RNA.at.F.U.Bx.Value=="<43") | 
              (Gt1b$HCV.RNA.at.F.U.Bx.Value== "Neg") ]="0"
Gt1b$Treatment.Group="Treated"
Gt1b$Treatment.Group[(Gt1b$Treatment=="Base")]="Untreated" 

dim(Gt1b)
[1]    22 53625
Gt1b <- Gt1b[ ,c(1:6,53624:53625,7:53623)]
Gt1b[1:22,1:8]

Patient.ID Treatment Baseline.Fibrosis.Group HCV.Genotype Failed HCV.RNA.at.F.U.Bx.Value
1         110      Base                   F_5_6           1b      0                     <43
2         110       Wk4                   F_5_6           1b      0                     <43
7         513      Base                   F_0_2           1b      0                     Neg
8         513       Wk4                   F_0_2           1b      0                     Neg
11        814      Base                   F_0_2           1b      1                     249
12        814       Wk2                   F_0_2           1b      1                     249
17        975      Base                   F_3_4           1b      0                      46
18        975       Wk2                   F_3_4           1b      0                      46
25       1059      Base                   F_5_6           1b      1                     <43
26       1059       Wk2                   F_5_6           1b      1                     <43
27       1154      Base                   F_3_4           1b      0                     Neg
28       1154       Wk4                   F_3_4           1b      0                     Neg
29       1230      Base                   F_0_2           1b      0                     <43
30       1230       Wk4                   F_0_2           1b      0                     <43
35       1505      Base                   F_3_4           1b      0                     <43
36       1505       Wk4                   F_3_4           1b      0                     <43
49       2168      Base                   F_5_6           1b      1                     <43
50       2168       Wk2                   F_5_6           1b      1                     <43
55       2783      Base                   F_5_6           1b      1                     101
56       2783       Wk2                   F_5_6           1b      1                     101
63       2957      Base                   F_3_4           1b      0                     Neg
64       2957       Wk4                   F_3_4           1b      0                     Neg
Virus Treatment.Group
1      0       Untreated
2      0         Treated
7      0       Untreated
8      0         Treated
11     1       Untreated
12     1         Treated
17     1       Untreated
18     1         Treated
25     0       Untreated
26     0         Treated
27     0       Untreated
28     0         Treated
29     0       Untreated
30     0         Treated
35     0       Untreated
36     0         Treated
49     0       Untreated
50     0         Treated
55     1       Untreated
56     1         Treated
63     0       Untreated
64     0         Treated


Gt1b_Base=Gt1b[Gt1b$Treatment=="Base",]
dim(Gt1b_Base)
[1]    11 53625

Failed=as.factor(Gt1b_Base$Failed)
Virus=as.factor(Gt1b_Base$Virus)
table(Failed,Virus) ##  Looks like 2 with virus failed and 1 with virus was successful

### 3 with virus and 9 without - 2 failed and 1 succeeded.
Virus
Failed 0 1
0 6 1
1 2 2

Gt1b.mat=(as.matrix(Gt1b[, 9:53625]))
pca=prcomp(Gt1b.mat)

TreatVirus <- factor(paste(Gt1b$Virus,Gt1b$Treatment.Group,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(TreatVirus)), 
       type="s",size=2)
group.v<-as.vector(TreatVirus)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
rgl.postscript("figures/pca3d_Gt1b_TreatVirus.pdf","pdf")

#####miRNA Analysis ##################

#Global miRNA heatmap

#Pathway for miRNAs

miRdat=read.csv("miR_rawdata.txt",sep="\t",header = TRUE)
Pheno=read.csv("miR_Pheno.txt",sep="\t",strip.white = TRUE,header=TRUE)

negs=miRdat[miRdat$Class.Name=="Negative",4:21]
negMean=colMeans(negs)
meannegMean=mean(negMean)
meannegMean=meannegMean+2*sd(negMean)

tmiRdat = t(miRdat)

dim(tmiRdat)
dimnames(tmiRdat)[[2]]<-tmiRdat[1,]

myval= tmiRdat[3,] == "Endogenous1"
tmiRdat_sub=tmiRdat[,myval==TRUE]

tmiRdat_sub=as.data.frame(tmiRdat_sub[4:21,])
tmiRdat_sub$name=rownames(tmiRdat_sub)
i <- sapply(Pheno, is.factor)
Pheno[i] <- lapply(Pheno[i], as.character)

tmiRdat_sub_merge=merge(Pheno,tmiRdat_sub,
      by.x="Sample",by.y="name")
dim(tmiRdat_sub_merge)
tmiRdat_sub_merge$SampleID=paste(tmiRdat_sub_merge$Patient_ID,
  tmiRdat_sub_merge$Response,tmiRdat_sub_merge$Treatment,sep=".")

tmiRdat_sub_merge.mat=t(tmiRdat_sub_merge[,6:803])
class(tmiRdat_sub_merge.mat) <- "numeric"
filter <- apply(tmiRdat_sub_merge.mat, 1, function(x) length(x[x>meannegMean])>=2)
tmiRdat_sub_merge.mat.filt <- tmiRdat_sub_merge.mat[filter,]
dim(tmiRdat_sub_merge.mat.filt)
[1] 301  18

#Draw density plot

miRdat_sub_merge=tmiRdat_sub_merge.mat.filt 
head(miRdat_sub_merge)
dimnames(miRdat_sub_merge)[[2]] <- tmiRdat_sub_merge[,804]
class(miRdat_sub_merge) <- "numeric"
miRdat_sub_merge=log(miRdat_sub_merge,2)
miRdat_sub_merge_norm=normalizeQuantiles(miRdat_sub_merge)

#miRdat_sub_merge_norm=normalizeMedianValues(miRdat_sub_merge)

df.m <- melt.data.frame(as.data.frame(miRdat_sub_merge_norm))

ggplot(df.m) + 
  geom_density(aes(x = value, colour = variable)) + labs(x = NULL) +
  theme(legend.position='right') + scale_x_log10() + ggtitle("Raw Counts")

par(mar=c(10,7,1,1))
boxplot(log(value)~variable,las=2,data=df.m,main="Raw Signal", 
        ylab="Counts",col=c(2,2,3,3,4,4))


Treatment=as.factor(tmiRdat_sub_merge$Treatment)
#Cirrhosis.Group.Gt1b=as.factor(Gt1b$Cirrhosis.Group)
#table(Treatment.Time.Gt1b,Cirrhosis.Group.Gt1b)
Patient <- factor(tmiRdat_sub_merge$Patient_ID) 

table(Patient,Treatment)

#Treatment
#Patient Treated Untreated
110        1         1
814        1         1
975        1         1
1154       1         1
1230       1         1
1505       1         1
2168       1         1
2783       1         1
2957       1         1


Failure=as.factor(tmiRdat_sub_merge$Response)
table(Treatment,Failure)
Failure
Treatment   NR R
Treated    3 6
Untreated  3 6

Treatment.Fail=as.factor(paste(Failure,Treatment,sep="."))

Treatment.Fail
[1] R.Untreated  R.Treated    R.Untreated  R.Treated    R.Untreated  R.Treated    R.Untreated 
[8] R.Treated    NR.Treated   NR.Untreated R.Treated    R.Untreated  NR.Treated   NR.Untreated
[15] R.Untreated  R.Treated    NR.Untreated NR.Treated  
Levels: NR.Treated NR.Untreated R.Treated R.Untreated

tedf= t(miRdat_sub_merge_norm)
pca=prcomp(tedf,scale.=T)
tedf1 = data.frame(tedf)
Phenotype=Patient
Phenotype=Treatment.Fail
plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
pca$x[,1:3] 

Pheno$Treatment
cell_rep=paste(Phenotype,Patient,sep=".")
cell_rep=paste(Phenotype,Failure,Treatment,sep=".")
tedf1$group = as.factor(Patient)
plot3d(pca$x[,1:3],col = as.integer(tedf1$group),type="s",size=2)
group.v<-as.vector(cell_rep)
text3d(pca$x[,1:3], pca$y, pca$z, group.v, cex=1.0, adj = 1.2) 
rgl.postscript("pca3d_miR.pdf","pdf")


#design <- model.matrix(~Patient+Treatment) 

dim(miRdat_sub_merge_norm)
#tmiRdat_sub_merge.mat=t(as.matrix(tmiRdat_sub_merge[, 6:803]))
#dim(tmiRdat_sub_merge.mat)
[1] 332  18
[1] 404  18  # Second try - used a less stringent negative control threshold

class(tmiRdat_sub_merge.mat) <- "numeric"


design <- model.matrix(~0+Treatment.Fail) 
colnames(design)=levels(Treatment.Fail)
corfit <- duplicateCorrelation(miRdat_sub_merge_norm,design,
                               block=tmiRdat_sub_merge$Patient_ID)
corfit$consensus
fit <- lmFit(miRdat_sub_merge_norm,design,
            block=tmiRdat_sub_merge$Patient_ID,
            correlation=corfit$consensus)


cm <- makeContrasts(NRvsR.Untreated = NR.Untreated-R.Untreated,
                    NRvsR.Treated = NR.Treated-R.Treated,
                    NR.TreatedvsUntreated = NR.Treated-NR.Untreated,
                    R.TreatedvsUntreated = R.Treated-R.Untreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="R.TreatedvsUntreated")
volcanoplot(fit2,coef="R.TreatedvsUntreated")
finalresmir=topTable(fit2,sort="none",n=Inf)

FC = 2^(fit2$coefficients)
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")


miR_select=c("hsa-miR-3144-5p","hsa-miR-543","hsa-miR-149-5p","hsa-miR-654-5p","hsa-miR-539-5p","hsa-miR-375","hsa-miR-147a","hsa-miR-140-3p","hsa-miR-193a-3p","hsa-miR-548q","hsa-miR-423-3p","hsa-miR-2682-5p","hsa-miR-107","hsa-miR-199b-5p","hsa-miR-575","hsa-miR-9-5p","hsa-miR-4286","hsa-miR-455-3p","hsa-miR-513a-5p","hsa-miR-874-3p","hsa-miR-185-5p","hsa-miR-135b-5p","hsa-miR-423-5p","hsa-miR-99b-5p","hsa-miR-520f-3p","hsa-miR-497-5p","hsa-miR-331-3p","hsa-miR-210-3p","hsa-miR-4536-5p","hsa-miR-642a-5p","hsa-miR-24-3p","hsa-miR-483-3p","hsa-miR-93-5p","hsa-miR-22-3p","hsa-miR-28-5p")
miR_select_D=c("hsa-miR-140-3p","hsa-miR-2682-5p","hsa-miR-199b-5p","hsa-miR-455-3p","hsa-miR-513a-5p","hsa-miR-185-5p")
finalresmir=cbind(FC,pvalall,finalresmir)
write.table(finalresmir,file="finalres_mir_FC.txt",col.names = TRUE,
            row.names = TRUE, sep = "\t")


NR.TxvsUntx=rownames(finalresmir[finalresmir$NR.TreatedvsUntreated_pval<0.05 & finalresmir$NR.TreatedvsUntreated_FC>0,])
R.TxvsUntx=rownames(finalresmir[finalresmir$R.TreatedvsUntreated_pval<0.05 & finalresmir$R.TreatedvsUntreated_FC>0,])
NRvsR.Tx=rownames(finalresmir[finalresmir$R.TreatedvsUntreated_pval<0.05 & finalresmir$NR.TreatedvsUntreated_FC>0,])
NRvsR.Untx=rownames(finalresmir[finalresmir$NRvsR.Untreated_pval<0.05 & finalresmir$NRvsR.Untreated_FC>0,])

#Selected miRNA heat map (ones that were significant)

palette <- colorRampPalette(c('blue','red'))(12)
palette2 <- c('green','yellow','red','blue')

palette

sigmiR=c(NRvsR.Untx,NR.TxvsUntx,R.TxvsUntx,NRvsR.Tx)
sigmiRdat=miRdat_sub_merge_norm[rownames(miRdat_sub_merge_norm) %in% sigmiR,]
sigmiRdat=miRdat_sub_merge_norm[rownames(miRdat_sub_merge_norm) %in% R.TxvsUntx,]

colnames(miRdat_sub_merge_norm)=Treatment.Fail 
colnames(sigmiRdat)=Treatment.Fail

miRdat_sub_merge_norm=sigmiRdat   # Trade for the new one
color=as.numeric(as.factor(colnames(miRdat_sub_merge_norm)))
f.patient=as.factor(colnames(miRdat_sub_merge_norm))

hc.rows <- hclust(dist(miRdat_sub_merge_norm))

heatmap(as.matrix(miRdat_sub_merge_norm), 
        ColSideColors = palette2[color],
        col=palette, labCol= NA, 
        Rowv=as.dendrogram(hc.rows),
        #labRow=NA,
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(miRdat_sub_merge_norm)),
                                  method="average")))
legend("topright",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.7)

miR_mRNA=read.table("all_miRNA_filt.txt",header=TRUE,sep="\t",fill=TRUE)
head(miR_mRNA)
miRNA_count = group_by(miR_mRNA, miRNA)
miRNA_count_stat = summarise(miRNA_count, count = n())

NR.TxvsUntx_miR=miR_mRNA[miR_mRNA$miR %in% NR.TxvsUntx,]
R.TxvsUntx_miR=miR_mRNA[miR_mRNA$miR %in% R.TxvsUntx,]


R.TxvsUntx_mRNA=finalres[(finalres$TreatedvsUntreatedForSuccess_pval<0.05 &
                finalres$TreatedvsUntreatedForSuccess_FC<0),c(25,29,36)]

R.TxvsUntx_miR_mRNA=inner_join(x = R.TxvsUntx_mRNA, y = NR.TxvsUntx_miR, by=c("Gene" = "mRNA"))

test=R.TxvsUntx_miR_mRNA %>% group_by(miRNA) %>% mutate(count = n())

by_miRNA <- group_by(R.TxvsUntx_miR_mRNA, miRNA)
miRNA_count_R.TxvsUntx <- summarise(by_miRNA, mean(TreatedvsUntreatedForSuccess_FC),count = n())
         
x=data.frame(miRNA_count_stat)[data.frame(miRNA_count_stat)$miRNA %in% 
              data.frame(miRNA_count_R.TxvsUntx)$miRNA==TRUE,]
y=data.frame(miRNA_count_R.TxvsUntx)

######NanoString mRNA Data #######################


NanomRNA=read.csv("mRNA_normdata.txt",sep="\t",header = TRUE)
mRNA_Pheno=read.csv("mRNA_Pheno.txt",sep="\t",header=TRUE)



#Heat map of failed vs. success untreated from the 46 selected genes (from nano string mRNA chip data)

negs_mRNA=NanomRNA[NanomRNA$Class.Name=="Negative",4:24]
negMean=colMeans(negs_mRNA)
meannegMean=mean(negMean)
meannegMean=meannegMean+2*sd(negMean)

tNanomRNA = t(NanomRNA)

dim(tNanomRNA)
dimnames(tNanomRNA)[[2]]<-tNanomRNA[1,]

myval= tNanomRNA[3,] == "Endogenous"
tNanomRNA_sub=tNanomRNA[,myval==TRUE]
tNanomRNA_sub=as.data.frame(tNanomRNA_sub[4:24,])
tNanomRNA_sub$name=rownames(tNanomRNA_sub)
mRNA_Pheno$Sample=paste("X",mRNA_Pheno$Sample,sep="")

tNanomRNA_sub_merge=merge(mRNA_Pheno,tNanomRNA_sub,
                        by.x="Sample",by.y="name")

tNanomRNA_sub_merge$SampleID=paste(tNanomRNA_sub_merge$Patient_ID,
        tNanomRNA_sub_merge$Response,tNanomRNA_sub_merge$Treatment,sep=".")
head(tNanomRNA_sub_merge)
dim(tNanomRNA_sub_merge)
tNanomRNA_sub_merge.mat=t(tNanomRNA_sub_merge[,6:52])
class(tNanomRNA_sub_merge.mat) <- "numeric"
tail(tNanomRNA_sub_merge.mat)

#Draw density plot

NanomRNA_sub_merge=tNanomRNA_sub_merge.mat
dimnames(NanomRNA_sub_merge)[[2]] <- tNanomRNA_sub_merge[,53]
class(NanomRNA_sub_merge) <- "numeric"
NanomRNA_sub_merge=log(NanomRNA_sub_merge,2)

NanomRNA_sub_merge_norm=normalizeQuantiles(NanomRNA_sub_merge)

#miRdat_sub_merge_norm=normalizeMedianValues(miRdat_sub_merge)

df.m <- melt.data.frame(as.data.frame(NanomRNA_sub_merge))

ggplot(df.m) + 
  geom_density(aes(x = value, colour = variable)) + labs(x = NULL) +
  theme(legend.position='right') + scale_x_log10() + ggtitle("Raw Counts")

par(mar=c(10,7,1,1))
boxplot(log(value)~variable,las=2,data=df.m,main="Raw Signal", 
        ylab="Counts",col=c(2,2,2,3,3,3))

Treatment=as.factor(tNanomRNA_sub_merge$Treatment)
#Cirrhosis.Group.Gt1b=as.factor(Gt1b$Cirrhosis.Group)
#table(Treatment.Time.Gt1b,Cirrhosis.Group.Gt1b)
Patient <- factor(tNanomRNA_sub_merge$Patient_ID) 
Failure <- factor(tNanomRNA_sub_merge$Response)

table(Patient,Treatment)
Treatment
Patient Treated Untreated
110        1         1
513        1         1
814        1         0
975        1         1
1059       1         1
1154       1         1
1230       1         1
1505       1         1
2168       1         1
2783       1         1
2957       1         1
#design <- model.matrix(~Patient+Treatment) 

tedf= t(NanomRNA_sub_merge)
pca=prcomp(tedf,scale.=T)
tedf1 = data.frame(tedf)
Phenotype=Patient
#Phenotype=Treatment.Fail
plot(pca,type="lines")  #Decide how many PC's are relevant for plotting
pca$x[,1:3] 

cell_rep=paste(Patient,Failure,Treatment,sep=".")
Phenotype=Patient
tedf1$group = as.factor(Phenotype)
plot3d(pca$x[,1:3],col = as.integer(tedf1$group),type="s",size=2)
group.v<-as.vector(cell_rep)
text3d(pca$x[,1:3], pca$y, pca$z, group.v, cex=1.0, adj = 1) 


#Heatmap of treated and untreated Subjects
palette <- colorRampPalette(c('blue','red'))(12)
palette2 <- c('green','yellow','red','blue')

palette
TreatFail=paste(tNanomRNA_sub_merge$Response,tNanomRNA_sub_merge$Treatment,sep=".")
color=as.numeric(as.factor(TreatFail))
f.patient=as.factor(TreatFail)

hc.rows <- hclust(dist(NanomRNA_sub_merge))

heatmap(as.matrix(NanomRNA_sub_merge), 
        ColSideColors = palette2[color],
        col=palette, labCol= NA, labRow=NA,
        Rowv=as.dendrogram(hc.rows),
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(NanomRNA_sub_merge)),
                                  method="average")))

legend("topright",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.7)

#Heatmap of untreated Subjects only
palette <- colorRampPalette(c('blue','red'))(12)
palette2 <- c('orange','mediumblue','red','green')

palette

tNanomRNA_sub_merge_Untx=filter(tNanomRNA_sub_merge,
                        tNanomRNA_sub_merge$Treatment=="Untreated")
NanomRNA_sub_merge_Untx=t(tNanomRNA_sub_merge_Untx[,6:52])
class(NanomRNA_sub_merge_Untx) <- "numeric"

TreatFail=paste(tNanomRNA_sub_merge_Untx$Response,tNanomRNA_sub_merge_Untx$Treatment,sep=".")
TreatFail
color=as.numeric(as.factor(TreatFail))
f.patient=as.character(TreatFail)
f.patient[f.patient=="R.Untreated"]="Responder Baseline"
f.patient[f.patient=="NR.Untreated"]="Failure Baseline"
f.patient = as.factor(f.patient)


hc.rows <- hclust(dist(NanomRNA_sub_merge_Untx))

heatmap(as.matrix(NanomRNA_sub_merge_Untx), 
        ColSideColors = palette2[color],
        col=palette, labCol= NA, 
        Rowv=as.dendrogram(hc.rows),
        #Colv = NA)
        Colv=as.dendrogram(hclust(dist(t(NanomRNA_sub_merge_Untx)),
                                  method="average")))

legend("right",legend=levels(f.patient), fill=palette2,title="Group",
       cex=0.8)

####Confirmation dotplots###############
#mRNA confirmation plots of baseline, wk2 and wk 4 for failed and success (one slide with 6 genes (CXCL10, CXCL11, RSAD2, ISG15, DDX60, IRF9)
options(digits=9)
ISG_6=c("CXCL10","CXCL11","RSAD2","ISG15","DDX60","IRF9")
Treatment=as.factor(mRNA_Pheno$Treatment.Time[2:22])
Failed=as.factor(mRNA_Pheno$Response)
#ISG_annotsub=subset(annot,annot$SYMBOL %in% ISG == TRUE)
#ISG_sub=Gt1b[,match(ISG_annotsub$PROBEID,colnames(Gt1b))]
tNanomRNA_sub_merge
NanomRNA6=as.data.frame(tNanomRNA_sub_merge[,colnames(tNanomRNA_sub_merge) %in% ISG_6])

cols = c(1:6);    
NanomRNA6[,cols] = apply(NanomRNA6[,cols], 2, function(x) as.numeric(as.character(x)))

NanomRNA6=cbind(tNanomRNA_sub_merge[,1:6],NanomRNA6)
NanomRNA6
#NanomRNA6=NanomRNA6[NanomRNA6$Treatment=="Untreated",]  # Only Untreated

NanomRNA6$Response2[NanomRNA6$Response=="R"]="Responder"
NanomRNA6$Response2[NanomRNA6$Response=="NR"]="Failure"
NanomRNA6$Treatment2[NanomRNA6$Treatment=="Untreated"]="Baseline"
NanomRNA6$Treatment2[NanomRNA6$Treatment=="Treated"]="On-Tx"
length = dim(NanomRNA6)[2]-8
for(i in 1:length){
  j=i+6
  filename=paste(names(NanomRNA6)[j],"nano3",sep="_")
  #png(paste(filename,"png",sep="."))
  ymin=floor(min(NanomRNA6[,j]))
  ymax=ceiling(max(NanomRNA6[,j]))
  size=(ymax-ymin)*0.05
  plots=ggplot(NanomRNA6, aes(x=Treatment2, y=NanomRNA6[,j])) + 
    geom_dotplot(aes(fill=factor(Response2)),
          binaxis='y', stackdir='center', binwidth = size) +
    scale_fill_manual(values=c("orange", "blue")) +
    xlab("") + ylab(names(NanomRNA6)[j]) + ylim(c(ymin,ymax)) + 
    theme(text = element_text(size=26)) +
    theme_bw() + theme(legend.text = element_text(size = 16)) +
    theme(legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=16)) +
    theme(axis.title.y  = element_text(size=16))
  plots
  ggsave(plots,filename=paste(filename,"png",sep="."))
  dev.off()
}

NanomRNA6[,7:12]=apply(NanomRNA6[,7:12], 2, scale) 
df.m=melt(NanomRNA6[,7:13], id="Response2")
plots=ggplot(df.m, aes(x=variable, y=value)) + 
  geom_dotplot(aes(fill=factor(Response2)),
    binaxis='y', stackdir='center') +
  scale_fill_manual(values=c("orange", "blue")) + theme_bw() +
  theme(text = element_text(size=20)) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x  = element_text(size=16))
plots

length(Treatment)

dim(miRdat_sub_merge_norm)
#tmiRdat_sub_merge.mat=t(as.matrix(tmiRdat_sub_merge[, 6:803]))
#dim(tmiRdat_sub_merge.mat)
[1] 332  18
[1] 404  18  # Second try - used a less stringent negative control threshold

class(tmiRdat_sub_merge.mat) <- "numeric"
Failure=as.factor(tmiRdat_sub_merge$Response)
table(Treatment,Failure)
Failure
Treatment   NR R
Treated    3 6
Untreated  3 6

Treatment.Fail=as.factor(paste(Failure,Treatment,sep="."))

Treatment.Fail
[1] R.Untreated  R.Treated    R.Untreated  R.Treated    R.Untreated  R.Treated    R.Untreated 
[8] R.Treated    NR.Treated   NR.Untreated R.Treated    R.Untreated  NR.Treated   NR.Untreated
[15] R.Untreated  R.Treated    NR.Untreated NR.Treated  
Levels: NR.Treated NR.Untreated R.Treated R.Untreated

design <- model.matrix(~0+Treatment.Fail) 
colnames(design)=levels(Treatment.Fail)
corfit <- duplicateCorrelation(miRdat_sub_merge_norm,design,
                               block=tmiRdat_sub_merge$Patient_ID)
corfit$consensus
fit <- lmFit(miRdat_sub_merge_norm,design
             block=tmiRdat_sub_merge$Patient_ID,
             correlation=corfit$consensus)


cm <- makeContrasts(NRvsR.Untreated = NR.Untreated-R.Untreated,
                    NRvsR.Treated = NR.Treated-R.Treated,
                    NR.TreatedvsUntreated = NR.Treated-NR.Untreated,
                    R.TreatedvsUntreated = R.Treated-R.Untreated,
                    levels=design) 
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="R.TreatedvsUntreated")
volcanoplot(fit2,coef="R.TreatedvsUntreated")
finalresmir=topTable(fit2,sort="none",n=Inf)

FC = 2^(fit2$coefficients)
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(colnames(FC),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")

finalresmir=cbind(FC,pvalall,finalresmir)
write.table(finalresmir,file="finalres_mir_FC.txt",col.names = TRUE,
            row.names = TRUE, sep = "\t")


NR.TxvsUntx=rownames(finalresmir[finalresmir$NR.TreatedvsUntreated_pval<0.05 & finalresmir$NR.TreatedvsUntreated_FC>0,])
R.TxvsUntx=rownames(finalresmir[finalresmir$R.TreatedvsUntreated_pval<0.05 & finalresmir$R.TreatedvsUntreated_FC>0,])
NRvsR.Tx=rownames(finalresmir[finalresmir$R.TreatedvsUntreated_pval<0.05 & finalresmir$NR.TreatedvsUntreated_FC>0,])
NRvsR.Untx=rownames(finalresmir[finalresmir$NRvsR.Untreated_pval<0.05 & finalresmir$NRvsR.Untreated_FC>0,])





miR_mRNA=read.table("all_miRNA_filt.txt",header=TRUE,sep="\t",fill=TRUE)
head(miR_mRNA)
miRNA_count = group_by(miR_mRNA, miRNA)
miRNA_count_stat = summarise(miRNA_count, count = n())

NR.TxvsUntx_miR=miR_mRNA[miR_mRNA$miR %in% NR.TxvsUntx,]
R.TxvsUntx_miR=miR_mRNA[miR_mRNA$miR %in% R.TxvsUntx,]


R.TxvsUntx_mRNA=finalres[(finalres$TreatedvsUntreatedForSuccess_pval<0.05 &
                            finalres$TreatedvsUntreatedForSuccess_FC<0),c(25,29,36)]

R.TxvsUntx_miR_mRNA=inner_join(x = R.TxvsUntx_mRNA, y = NR.TxvsUntx_miR, by=c("Gene" = "mRNA"))

test=R.TxvsUntx_miR_mRNA %>% group_by(miRNA) %>% mutate(count = n())

by_miRNA <- group_by(R.TxvsUntx_miR_mRNA, miRNA)
miRNA_count_R.TxvsUntx <- summarise(by_miRNA, mean(TreatedvsUntreatedForSuccess_FC),count = n())

x=data.frame(miRNA_count_stat)[data.frame(miRNA_count_stat)$miRNA %in% 
                                 data.frame(miRNA_count_R.TxvsUntx)$miRNA==TRUE,]
y=data.frame(miRNA_count_R.TxvsUntx)
