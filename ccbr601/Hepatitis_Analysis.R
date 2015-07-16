Ghany

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
tedf2 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","Timing"), 
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


tedf2 <- ddply(tedf, c("Patient.ID", "Treatment","Baseline.Fibrosis.Group","Timing","Failed"), 
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

design <- model.matrix(~Patient+Treatment) 
Gt1b.mat=t(as.matrix(Gt1b[, 9:53625]))
fit <- lmFit(Gt1b.mat, design) 
fit <- eBayes(fit) 
topTable(fit, coef="TreatmentUntreated")  
finalres=topTable(fit,coef="TreatmentUntreated",sort="none",n=Inf)

FC = 2^finalres[,1]
FC = 1/FC
FC = ifelse(FC<1,-1/FC,FC)

finalres = as.data.frame(cbind(FC,finalres))

dim(finalres)

##########  Look at Treatment Failure within Gt1b ############
dim(Gt1b)
Gt1b <- Gt1b[ ,c(1:5,53623:53625,6:53622)]
Gt1b[1:10,1:10]
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
topTable(fit2, coef="TreatedvsUntreatedForSuccess")

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

FailedvsSuccessUntreated_pval=p.adjust(fit2$p.value[,1], method='BH')
FailedvsSuccessForTreated=p.adjust(fit2$p.value[,2], method='BH')
TreatedvsUntreatedForSuccess=p.adjust(fit2$p.value[,3], method='BH')
TreatedvsUntreatedForFailed=p.adjust(fit2$p.value[,4], method='BH')
pvaladjall=cbind(FailedvsSuccessUntreated_pval,FailedvsSuccessForTreated,
                 TreatedvsUntreatedForSuccess,TreatedvsUntreatedForFailed)
colnames(pvaladjall)=paste(colnames(pvaladjall),"adjpval",sep="_")

finalres = as.data.frame(cbind(FC,pvalall,pvaladjall))

FailedvsSuccessUntx=finalres[(finalres[,1]>=1.2 | finalres[,1]<=-1.2) && finalres[,3] < 0.1]

TxvsUntx_Success=finalres[((finalres$TreatedvsUntreatedForSuccess_FC>=1.2 | 
                             finalres$TreatedvsUntreatedForSuccess_FC<=-1.2) &
                             finalres$TreatedvsUntreatedForSuccess_pval < 0.1),]


TxvsUntx_Failed=finalres[((finalres[,4]>=1.2 | finalres[,4]<=-1.2) & finalres[,8] < 0.1),]



dim(finalres)





######### Draw PCA Plot ###############

setwd("/Volumes/CCBR/projects/ccbr601/Niddk4MGHany-5-6-15/04-24-14-MGhany_HA-HuGene2_ST/")
load("Gt1a.RData")

Gt1b.mat2=(as.matrix(Gt1b[, 9:53625]))
pca=prcomp(Gt1b.mat2)
Treat <- factor(paste(Gt1b$Cirrhosis.Group,Gt1b$Treatment.Time,sep=".")) 

library(rgl)
open3d() 

plot3d(pca$x[,1:3],col=as.numeric(as.factor(Failure)), 
       type="s",size=2)
group.v<-as.vector(Failure)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj = 1) 
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


################Re-do finalres6 - just treated vs untreated################

colnames(finalres6)
[1] "ID"        "FC"        "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val"
[8] "B"         "PROBEID"   "SYMBOL"    "REFSEQ"    "GENE"  

write.table(finalres6,file="all_results_Gt1b.txt",row.names=FALSE,sep="\t",quote=FALSE)


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





