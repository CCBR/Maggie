
#BRCA, COAD, LUAD, LUSC, READ, SKCM, THCA
type="THCA"

p <- function(..., sep='') {
  paste(..., sep=sep, collapse=sep)
}

library(stringr)
setwd(p("/Users/maggiec/GitHub/Maggie/ccbr540/data/gdac.broadinstitute.org_",type,
        ".Merge_Clinical.Level_1.2015060100.0.0"))

#clin = read.csv("tp_clindata.txt", na.strings = "",header=TRUE, sep="\t")
#clin = read.csv("tp_BRCA_clinical.txt", na.strings = "",header=TRUE, sep="\t")
clin = read.csv(p("tp_",type,"_clinical.txt"), na.strings = "",header=TRUE, sep="\t")

# birth=clin$patient.days_to_birth
# relative=clin$patient.first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_types.first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_type
# patient=clin$patient.bcr_patient_barcode
# survival=clin$patient.clinical_cqcf.days_to_death
# ethnicity=clin$patient.ethnicity
# status=clin$patient.extrathyroid_carcinoma_present_extension_status
# dead=clin$patient.vital_status
# age=clin$patient.age_at_initial_pathologic_diagnosis
# gender=clin$patient.gender
# hist=clin$patient.histological_type
# metastasis=clin$patient.malignant_neoplasm_metastatic_involvement_site
# focal=clin$patient.clinical_cqcf.tumor_focality
# state=clin$patient.follow_ups.follow_up.person_neoplasm_cancer_status
# race=clin$patient.race
# stage=clin$patient.stage_event.pathologic_stage

#snpdat = read.table("/data/CCBR/projects/ccbr540/HABP2bams/HABP2_patient_snpstatus.txt", header=FALSE, fill=TRUE, sep="\t")
snpdat = read.table(p("../",type,"_sample.txt"), header=TRUE, fill=TRUE, sep="\t")
#colnames(snpdat)=c("pat_id","snp")


id=str_split_fixed(snpdat$pat_id, "-", 4)
snpdat$patid2=tolower(paste(id[,1],id[,2],id[,3],sep="-"))
# patdata=as.data.frame(cbind(as.character(patient),as.character(relative),
#                             survival,
#                             as.character(ethnicity),as.character(status),
#                             as.character(dead),age,as.character(gender),
#                             as.character(hist),as.character(metastasis),
#                             as.character(focal),as.character(state),as.character(race),
#                             as.character(stage)))
# colnames(patdata)=c("pat","relative","survival","ethnicity","status","dead",
#                     "age","gender","hist","metastasis","focal","state","race","stage")
# mergedat=merge(patdata,snpdat,by.x="pat",by.y="patid2")
# 
# snpdat$snp[snpdat[]="432614"
# #snpid <- substr(mergedat2$snp, 1, 3)

snpdat$snpid="0/0"

for (i in 2:7){
snpgt=strsplit(as.character(snpdat[,i]), split=":")
snpdat$snp=sapply(snpgt,function(x) x[1])
snpdat$snpid[snpdat$snp=="0/1"]="0/1"
}

snpdat=snpdat[!duplicated(snpdat$patid2), ]

dim(snpdat)
[1] 1057   10  #BRCA
[1] 299  10  #COAD
[1] 573  10   #LUAD -need to redo
[1] 484  10  #LUSC
[1] 114  10  #READ
[1] 470  10 #SKCM
[1] 373  10  #THCA

dim(snpdat[snpdat$snpid=="0/1",])
[1] 87 10  # BRCA
[1] 27 10 # COAD
[1] 82 10  #LUAD
[1] 46 10  #LUSC
[1]  6 10  #READ
[1] 53 10  #SKCM
[1] 29 10  #THCA

dim(snpdat[snpdat$snpid=="0/0",])
[1] 1027   10  #BRCA
[1] 284  10  #COAD
[1] 597  10  #LUAD
[1] 438  10  #LUSC
[1] 108  10  #READ
[1] 417  10  #SKCM
[1] 344  10  #THCA

mergedat2=merge(clin,snpdat,by.x="patient.bcr_patient_barcode",by.y="patid2")

dim(clin)
[1] 1085 3670  #BRCA
[1]  453 3151  #COAD
[1]  521 3011  #LUAD
[1]  495 2694  #LUSC
[1]  171 2742  #READ
[1]  469 1877  #SKCM
[1] 501 569  #THCA

dim(mergedat2)
[1] 1044 3679  #BRCA
[1]  292 3160  #COAD
[1]  513 3020  #LUAD
[1]  476 2703 #LUSC
[1]  114 2751  #READ
[1]  469 1886 #SKCM
[1] 364 578  #THCA

#snpid <- substr(mergedat$snpid, 1, 3)
#mergedat$SNP <- as.factor(ifelse(snpid == "0/0",0,1))
#mergedat$age2 <- as.factor(ifelse(as.matrix(mergedat$age) <= 60,0,1))

snpid=mergedat2$snpid
mergedat2$SNP <- as.factor(ifelse(snpid == "0/0",0,1))
#mergedat$age2 <- as.factor(ifelse(as.matrix(mergedat$age) <= 60,0,1))

mergedat2$age2 <- as.factor(ifelse(mergedat2$patient.age_at_initial_pathologic_diagnosis <= 
                              60,0,1))
#mergedat2$gt <- as.factor(ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA", 
#                                     0,1))
#mergedat2$gt2 <- as.factor(ifelse(mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="not identified" | 
#                      mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="not identifed",
#                      0,ifelse(mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="NA","NA",1)))


#mergedat2$gt <- as.factor(ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="not identified" |
#                              mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="negative", 
#                                 0,ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA","NA",1)))

#mergedat2$stage <- as.factor(ifelse(mergedat2$patient.stage_event.pathologic_stage=="stage i" |
#                                   mergedat2$patient.stage_event.pathologic_stage=="stage ii", 
#                                 "low",ifelse(mergedat2$patient.stage_event.pathologic_stage=="stage iii" |
#                                    mergedat2$patient.stage_event.pathologic_stage=="stage iva" |
#                                      mergedat2$patient.stage_event.pathologic_stage=="stage ivc",
#                                   "high","NA")))

#mergedat2$gt[mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA"]="0" 

# snprel=table(mergedat$relative,mergedat$SNP)
# snpstage=table(mergedat$stage,mergedat$SNP)
# metstage=table(mergedat$met,mergedat$SNP)
# genstage=table(mergedat$gen,mergedat$SNP)
# statstage=table(mergedat$status,mergedat$SNP)
# racestage=table(mergedat$race,mergedat$SNP)
# survstage=table(mergedat$surv,mergedat$SNP)
# tab=table(mergedat2$gt,mergedat2$SNP)
# 
# fisher.test(as.matrix(racestage))


result=vector()

for (i in 2:ncol(mergedat2)) {
#for (i in 2:100) {
  if (length(levels(mergedat2[,i])) <=10 & length(levels(mergedat2[,i])) > 1){
  tab=table(mergedat2[,i],mergedat2$SNP)
  tab
  tabname=colnames(mergedat2)[i]
  tabname
  if (tab[1,1]=="0" & tab[1,2]=="0" | tab[2,1]=="0" & tab[2,2]=="0") {
    fishtest$p.value=1}
  else{fishtest=fisher.test(tab,simulate.p.value=TRUE)}
  categ=paste(levels(mergedat2[,i]),",",collapse="")
  res=c(i,tabname,fishtest$p.value,categ)
  result=rbind(result,res)
}
}

#######For THCA: fishers exact test=0.000999500249875062 ##########
> i=486
> tab=table(mergedat2[,i],mergedat2$SNP)
> tab

0   1
NA                  319  20
nodular hyperplasia  14   7
other, specify        4   0

result[result[,3]<=0.05,1:3]  
#result_HABP=subset(mergedat2,grepl("0/1", snp))
result_ALL=mergedat2

#write.table(result,"snp_assoc_clin_BRCA.txt",row.names=FALSE,sep="\t")
#write.table(result_HABP,"snp_assoc_clin_HABP.txt",row.names=FALSE,sep="\t")

write.table(result,p("../../results/snp_assoc_clin_",type,".txt"),
                     row.names=FALSE,sep="\t")
write.table(result_ALL,p("../../results/all_clin_",type,".txt"),
            row.names=FALSE,sep="\t")

#############################END########

tab=table(mergedat2[,i],mergedat2$SNP)
fishtest=fisher.test(tab)

dim(subset(mergedat2,grepl("0/1", snp)))
selected<-c("parent","child","sibling")
dim(mergedat2[mergedat2$relative %in% selected,])
mergedat2[mergedat2$dead =="dead",]

getwd()
/gpfs/gsfs4/users/CCBR/projects/ccbr540/gdac.broadinstitute.org_THCA.Merge_Clinical.Level_1.2015040200.0.0

setwd("/Volumes/CCBR/projects/ccbr540/")
setwd("gdac.broadinstitute.org_THCA.Merge_Clinical.Level_1.2015040200.0.0")
load(".RData")

head(patdata)
