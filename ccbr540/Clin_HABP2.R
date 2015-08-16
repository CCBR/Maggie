

clin = read.csv("tp_clindata.txt", na.strings = "",header=TRUE, sep="\t")

birth=clin$patient.days_to_birth
relative=clin$patient.first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_types.first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_type
patient=clin$patient.bcr_patient_barcode
survival=clin$patient.clinical_cqcf.days_to_death
ethnicity=clin$patient.ethnicity
status=clin$patient.extrathyroid_carcinoma_present_extension_status
dead=clin$patient.vital_status
age=clin$patient.age_at_initial_pathologic_diagnosis
gender=clin$patient.gender
hist=clin$patient.histological_type
metastasis=clin$patient.malignant_neoplasm_metastatic_involvement_site
focal=clin$patient.clinical_cqcf.tumor_focality
state=clin$patient.follow_ups.follow_up.person_neoplasm_cancer_status
race=clin$patient.race
stage=clin$patient.stage_event.pathologic_stage

snpdat = read.table("/data/CCBR/projects/ccbr540/HABP2bams/HABP2_patient_snpstatus.txt", header=FALSE, fill=TRUE, sep="\t")
colnames(snpdat)=c("pat_id","snp")
library(stringr)

id=str_split_fixed(snpdat$pat_id, "-", 4)
snpdat$patid2=tolower(paste(id[,1],id[,2],id[,3],sep="-"))
patdata=as.data.frame(cbind(as.character(patient),as.character(relative),
                            survival,
                            as.character(ethnicity),as.character(status),
                            as.character(dead),age,as.character(gender),
                            as.character(hist),as.character(metastasis),
                            as.character(focal),as.character(state),as.character(race),
                            as.character(stage)))
colnames(patdata)=c("pat","relative","survival","ethnicity","status","dead",
                    "age","gender","hist","metastasis","focal","state","race","stage")
mergedat=merge(patdata,snpdat,by.x="pat",by.y="patid2")
mergedat2=merge(clin,snpdat,by.x="patient.bcr_patient_barcode",by.y="patid2")
subset(mergedat,grepl("0/1", snp))

snpid <- substr(mergedat$snp, 1, 3)
mergedat$SNP <- as.factor(ifelse(snpid == "0/0",0,1))
mergedat$age2 <- as.factor(ifelse(as.matrix(mergedat$age) <= 60,0,1))

mergedat2$SNP <- as.factor(ifelse(snpid == "0/0",0,1))
mergedat$age2 <- as.factor(ifelse(as.matrix(mergedat$age) <= 60,0,1))

mergedat2$age2 <- as.factor(ifelse(mergedat2$patient.age_at_initial_pathologic_diagnosis <= 
                              60,0,1))
mergedat2$gt <- as.factor(ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA", 
                                     0,1))
mergedat2$gt2 <- as.factor(ifelse(mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="not identified" | 
                      mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="not identifed",
                      0,ifelse(mergedat2$patient.ras_family_gene_genotyping_outcome_lab_results_text=="NA","NA",1)))


mergedat2$gt <- as.factor(ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="not identified" |
                              mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="negative", 
                                 0,ifelse(mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA","NA",1)))

mergedat2$stage <- as.factor(ifelse(mergedat2$patient.stage_event.pathologic_stage=="stage i" |
                                   mergedat2$patient.stage_event.pathologic_stage=="stage ii", 
                                 "low",ifelse(mergedat2$patient.stage_event.pathologic_stage=="stage iii" |
                                    mergedat2$patient.stage_event.pathologic_stage=="stage iva" |
                                      mergedat2$patient.stage_event.pathologic_stage=="stage ivc",
                                   "high","NA")))



mergedat2$gt[mergedat2$patient.braf_gene_genotyping_outcome_lab_results_text=="NA"]="0" 

snprel=table(mergedat$relative,mergedat$SNP)
snpstage=table(mergedat$stage,mergedat$SNP)
metstage=table(mergedat$met,mergedat$SNP)
genstage=table(mergedat$gen,mergedat$SNP)
statstage=table(mergedat$status,mergedat$SNP)
racestage=table(mergedat$race,mergedat$SNP)
survstage=table(mergedat$surv,mergedat$SNP)
tab=table(mergedat2$gt,mergedat2$SNP)

fisher.test(as.matrix(racestage))
result=vector()

for (i in 2:ncol(mergedat2)) {
  if (length(levels(mergedat2[,i])) <=10 & length(levels(mergedat2[,i])) > 1){
  tab=table(mergedat2[,i],mergedat2$SNP)
  tab
  tabname=colnames(mergedat2)[i]
  tabname
  fishtest=fisher.test(tab)
  categ=paste(levels(mergedat2[,i]),",",collapse="")
  res=c(i,tabname,fishtest$p.value,categ)
  result=rbind(result,res)
}
}

result[result[,3]<=0.05,2:3]  
result_HABP=subset(mergedat2,grepl("0/1", snp))
result_ALL=mergedat2

write.table(result,"snp_assoc_clin.txt",row.names=FALSE,sep="\t")
write.table(result_HABP,"snp_assoc_clin_HABP.txt",row.names=FALSE,sep="\t")
write.table(result_ALL,"all_clin.txt",row.names=FALSE,sep="\t")

tab=table(mergedat2[,i],mergedat2$SNP)
fishtest=fisher.test(tab)

dim(subset(mergedat,grepl("0/0", snp)))
selected<-c("parent","child","sibling")
mergedat[mergedat$relative %in% selected,]
mergedat[mergedat$dead =="dead",]

getwd()
/gpfs/gsfs4/users/CCBR/projects/ccbr540/gdac.broadinstitute.org_THCA.Merge_Clinical.Level_1.2015040200.0.0

setwd("/Volumes/CCBR/projects/ccbr540/")
setwd("gdac.broadinstitute.org_THCA.Merge_Clinical.Level_1.2015040200.0.0")
load(".RData")

head(patdata)
