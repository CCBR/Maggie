

#Change directory
setwd("/Users/maggiec/GitHub/Maggie/ccbr540/")

####BRCA ####################
BRCA_blacklist=read.table("data/BRCA_blacklist_samples.txt")
BRCA_blacklist[] <- lapply(BRCA_blacklist, as.character)
BRCA_vcf=read.table("data/BRCA.vcf.txt",header=TRUE,comment.char="",check.names=FALSE)
dim(BRCA_vcf)
index=match(colnames(BRCA_vcf),BRCA_blacklist$V1)
BRCA_sub=BRCA_vcf[,is.na(index)]
BRCA_sub[] <- lapply(BRCA_sub, as.character)
dim(BRCA_sub)
x=dim(BRCA_sub)[1]
y=dim(BRCA_sub)[2]
BRCA_mat=matrix(data=NA,nrow=x,ncol=7)
for (i in 1:x){
no_gt=stringi::stri_detect_fixed(BRCA_sub[i,10:y],"./.")
hom_ref=stringi::stri_detect_fixed(BRCA_sub[i,10:y],"0/0")
het_alt=stringi::stri_detect_fixed(BRCA_sub[i,10:y],"0/1")
hom_alt=stringi::stri_detect_fixed(BRCA_sub[i,10:y],"1/1")
BRCA_mat[i,1]=sum(no_gt == TRUE)
BRCA_mat[i,2]=sum(hom_ref == TRUE)
BRCA_mat[i,3]=sum(het_alt == TRUE)
BRCA_mat[i,4]=sum(hom_alt == TRUE)
}

BRCA_mat[,5]=rowSums(BRCA_mat[,1:4])
BRCA_mat[,6]=apply(BRCA_mat[,3:4],1,sum)
BRCA_mat[,7]=(BRCA_mat[,6]/BRCA_mat[,5])*100
colnames(BRCA_mat)=c("no_gt","hom_ref","het_alt","hom_alt","total","alt","AF")
BRCA.df=cbind(BRCA_sub[,c(1:4,8)],BRCA_mat)
BRCA.df[BRCA.df$POS=="115348046",]
write.table(BRCA.df,file="Results/BRCA_AF.txt",sep="\t")

##### Other Cancers ###################

type="READ"
blacklist=read.table(paste("data/",type,"_blacklist_samples.txt",sep=""))
blacklist[] <- lapply(blacklist, as.character)
vcf=read.table(paste("data/",type,".vcf.txt",sep=""),header=TRUE,
               comment.char="",check.names=FALSE)
dim(vcf)
index=match(colnames(vcf),blacklist$V1)
sub=vcf[,is.na(index)]
sub[] <- lapply(sub, as.character)
dim(sub)
x=dim(sub)[1]
y=dim(sub)[2]
mat=matrix(data=NA,nrow=x,ncol=9)
for (i in 1:x){
  no_gt=stringi::stri_detect_fixed(sub[i,10:y],"./.")
  hom_ref=stringi::stri_detect_fixed(sub[i,10:y],"0/0")
  het_alt=stringi::stri_detect_fixed(sub[i,10:y],"0/1")
  hom_alt=stringi::stri_detect_fixed(sub[i,10:y],"1/1")
  mat[i,1]=sum(no_gt == TRUE)
  mat[i,2]=sum(hom_ref == TRUE)
  mat[i,3]=sum(het_alt == TRUE)
  mat[i,4]=sum(hom_alt == TRUE)
}

mat[,5]=rowSums(mat[,1:4])
mat[,6]=apply(mat[,3:4],1,sum)
mat[,7]=(mat[,6]/mat[,5])*100
mat[,8]=(mat[,3]/mat[,5])*100
mat[,9]=(mat[,4]/mat[,5])*100
colnames(mat)=c("no_gt","hom_ref","het_alt","hom_alt","total","alt","AF","%het_alt","%hom_alt")
df=cbind(sub[,c(1:4,8)],mat)
df$ID=type
HABP2_snp=df[df$POS=="115348046",c("ID","total","%hom_alt","%het_alt")]
HABP2_snp
write.csv()
