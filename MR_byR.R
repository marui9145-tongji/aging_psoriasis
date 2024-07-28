setwd("F:/R/mendel/telomore&pso")
getwd()
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github("MRCIEU/MRInstruments")
#install.packages("R.utils")
library(TwoSampleMR)
#remotes::install_github('MRCIEU/TwoSampleMR')
library(data.table)
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
library(tidyverse)
library(plyr)
library(ggplot2)
library(MendelianRandomization)

#telomere and psoriasis###############################################################################
#import GWAS summary of telomere length from UK biobank
telomere<-fread("F:/R/mendel/GWASSUMMARY/UKB_telomere_gwas_summarystats.tsv")

telomere<-data.frame(telomere)

head(telomere)

colnames(telomere)<-c("SNP","pval.exposure","chromosome","position.exposure","effect_allele.exposure",
                      "other_allele.exposure","eaf.exposure","beta.exposure","se.exposure")

telomere$id.exposure<-"Leukocyte telomere length"
telomere$exposure<-"Leukocyte telomere length"
telomere$samplesize.exposure<-464716

head(telomere)

#Select SNPs that are significantly or strongly correlated with GWAS exposure, n=38506
telomere<-subset(telomere,pval.exposure<5e-08)

#clump_data，n=153
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
telomere<-clump_data(telomere,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                     clump_p2 = 5e-08,pop = "EUR")

#import GWAS summary of psoriasis from Finngen Release 10
psoriasis<-fread("F:/R/mendel/GWASSUMMARY/finngen_R10_L12_PSORIASIS.tsv")
psoriasis<-data.frame(psoriasis)

head(psoriasis)

colnames(psoriasis)<-c("chromosome","position.outcome","other_allele.outcome",
                       "effect_allele.outcome","SNP","nearest_genes","pval.outcome","mlogp.outcome",
                       "beta.outcome","se.outcome","eaf.outcome","eafca.outcome","eafcon.outcome")

psoriasis$id.outcome<-"Psoriasis"
psoriasis$outcome<-"Psoriasis"
psoriasis$samplesize.outcome<-407876

#Searching for SNPs corresponding to exposure in the GWAS summary with 12 missing, n=141
total<-merge(telomere,psoriasis,by.x = "SNP",by.y = "SNP",all = F)

#Remove 1 SNP with GWAS significance related to the outcome, n=140
total<-subset(total,pval.outcome>5e-08)

#Extract exposure and outcome data separately
telomere<-total[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                   "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

psoriasis<-total[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                    "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(telomere)
head(psoriasis)

#Removed palindrome structure, regulated direction, deleted 6, n=134
#Due to the palindrome with moderate allele frequency, it has been deleted:rs2276182, rs2306646, rs4742448, rs56178008, rs6501181, rs670180
data_h<-harmonise_data(exposure_dat = telomere,outcome_dat = psoriasis,action = 2)
data_h<-data_h%>%subset(data_h$mr_keep==TRUE)

#steiger filter
data_h_steiger<-steiger_filtering(data_h)
data_h_steiger<-subset(data_h_steiger,steiger_dir==TRUE)
#Remove duplicates,n=134
data_h_steiger<-data_h_steiger[!duplicated(data_h_steiger$SNP),]

#Export the overall data to check for any confounding factors related to psoriasis, and identify 5 SNPs with confounding,
#rs2293579、rs66731853、rs429358、rs10805346、rs2230590
write.csv(data_h_steiger,"data_h_steiger.csv",quote=F,row.names=F)
#n=129
data_h_steiger<-fread("data_h_steiger.csv")

#F statistic
data_h_steiger$R2<-data_h_steiger$beta.exposure*data_h_steiger$beta.exposure*2*(data_h_steiger$eaf.exposure)*(1-data_h_steiger$eaf.exposure)
data_h_steiger$F<-(data_h_steiger$samplesize.exposure-2)*data_h_steiger$R2/(1-data_h_steiger$R2)

#Export F and R2 data without SNPs with F ≤ 10
write.csv(data_h_steiger[,c("SNP","R2","F")],"data_h_steiger_R2_F.csv",quote=F,row.names=F)

#MR
mr_outcome<-mr(data_h_steiger)
mr_outcome

#OR and 95%CI
OR <-generate_odds_ratios(mr_outcome)
OR

#export
write.csv(OR[,c("outcome","exposure","method","nsnp","or","or_lci95",
                "or_uci95","pval")],"data_MR.csv",quote=F,row.names=F)

#Scatter plot
p1<-mr_scatter_plot(mr_outcome,data_h_steiger)
p1[[1]]

#Sensitivity analysis - leave one out method
mr_outcome_loo<-mr_leaveoneout(data_h_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

#forest plot
mr_outcome_single<-mr_singlesnp(data_h_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

#funnel plot
p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

#merge plots
library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

#heterogeneity
#If Q_pval<0.05, it indicates the presence of outliers. Use the ivw_mre method for analysis
heterogeneity<-mr_heterogeneity(data_h_steiger)
heterogeneity
#export heterogeneity results
write.csv(heterogeneity[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity.csv",quote=F,row.names=F)

#ivw_mre
mr_outcome_mre<-mr(data_h_steiger,method_list = "mr_ivw_mre")
mr_outcome_mre
#OR and 95%CI
OR_mre <-generate_odds_ratios(mr_outcome_mre)
OR_mre#
#export 
write.csv(OR_mre[,c("outcome","exposure","method","nsnp","or","or_lci95",
                "or_uci95","pval")],"data_MR_mre.csv",quote=F,row.names=F)

#Multi effect test, combined with scatter plot, to determine whether the y-axis intercept is large or not. If it is large, it indicates that the exposure factor has not yet taken effect,
#The outcome has already been produced, which means there is significant confounding and the MR results are unreliable. So we need to ensure the level
pleiotropy <- mr_pleiotropy_test(data_h_steiger)
pleiotropy
#export 
write.csv(pleiotropy[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy.csv",quote=F,row.names=F)

#mrpresso
res_MRPRESSO<-mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        data=data_h_steiger, 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        NbDistribution = 10000,  
                        SignifThreshold = 0.05)
res_MRPRESSO
#export 
write.csv(res_MRPRESSO$`Main MR results`,"MRPRESSO_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_Outlier.csv",quote=F,row.names=F)

#Export Steiger data that needs to exclude outliers
write.csv(data_h_steiger,"data_h_steiger_exoutliers.csv",quote=F,row.names=F)

#Import and remove outlier SNPs（n=5,rs10774624、rs4695407、rs55747751、rs7772289、rs9398196）
data_h_steiger_exoutliers<-fread("data_h_steiger_exoutliers.csv")
head(data_h_steiger_exoutliers)

#MR
mr_outcome_exoutliers<-mr(data_h_steiger_exoutliers)
mr_outcome_exoutliers#

#OR and 95%CI
OR_exoutliers <-generate_odds_ratios(mr_outcome_exoutliers)
OR_exoutliers

#export
write.csv(OR_exoutliers[,c("outcome","exposure","method","nsnp","or","or_lci95",
                "or_uci95","pval")],"data_exoutliers_MR.csv",quote=F,row.names=F)

#scatter_plot
p1<-mr_scatter_plot(mr_outcome_exoutliers,data_h_steiger_exoutliers)
p1[[1]]

#leaveoneout
mr_outcome_loo<-mr_leaveoneout(data_h_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]
#如果横线与大多数结果不一致则表示该SNP存在较大异质性

#forest plot
mr_outcome_single<-mr_singlesnp(data_h_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

#funnel_plot
p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

#merge
library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

#heterogeneity
heterogeneity_exoutliers<-mr_heterogeneity(data_h_steiger_exoutliers)
heterogeneity_exoutliers

#export
write.csv(heterogeneity_exoutliers[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_exoutliers_heterogeneity.csv",quote=F,row.names=F)

#ivw_mre
mr_outcome_mre_exoutliers<-mr(data_h_steiger_exoutliers,method_list = "mr_ivw_mre")
mr_outcome_mre_exoutliers
#OR and 95%CI
OR_mre_exoutliers <-generate_odds_ratios(mr_outcome_mre_exoutliers)
OR_mre_exoutliers#查看OR结果
#export
write.csv(OR_mre_exoutliers[,c("outcome","exposure","method","nsnp","or","or_lci95",
                    "or_uci95","pval")],"data_exoutliers_MR_mre.csv",quote=F,row.names=F)

#pleiotropy
pleiotropy_exoutliers <- mr_pleiotropy_test(data_h_steiger_exoutliers)
pleiotropy_exoutliers

#export
write.csv(pleiotropy_exoutliers[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_exoutliers_pleiotropy.csv",quote=F,row.names=F)

#Bidirectional Two-sample MR. Subsequently, conduct reverse analysis
psoriasis<-fread("F:/R/mendel/GWASSUMMARY/finngen_R10_L12_PSORIASIS.tsv")
psoriasis<-data.frame(psoriasis)

head(psoriasis)

colnames(psoriasis)<-c("chromosome","position.exposure","other_allele.exposure",
                       "effect_allele.exposure","SNP","nearest_genes","pval.exposure","mlogp.exposure",
                       "beta.exposure","se.exposure","eaf.exposure","eafca.exposure","eafcon.exposure")

psoriasis$id.exposure<-"Psoriasis"
psoriasis$exposure<-"Psoriasis"
psoriasis$samplesize.exposure<-407876

head(psoriasis)

psoriasis<-subset(psoriasis,pval.exposure<5e-08)

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
psoriasis<-clump_data(psoriasis,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                     clump_p2 = 5e-08,pop = "EUR")#如果是亚洲人则为“EAS”

telomere<-fread("F:/R/mendel/GWASSUMMARY/UKB_telomere_gwas_summarystats.tsv")

telomere<-data.frame(telomere)

head(telomere)

colnames(telomere)<-c("SNP","pval.outcome","chromosome","position.outcome","effect_allele.outcome",
                      "other_allele.outcome","eaf.outcome","beta.outcome","se.outcome")

head(telomere)

telomere$id.outcome<-"Leukocyte telomere length"
telomere$outcome<-"Leukocyte telomere length"
telomere$samplesize.outcome<-464716

head(telomere)

total<-merge(psoriasis,telomere,by.x = "SNP",by.y = "SNP",all = F)

total<-subset(total,pval.outcome>5e-08)

psoriasis<-total[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                   "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

telomere<-total[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                    "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(psoriasis)
head(telomere)

data_h_f<-harmonise_data(exposure_dat = psoriasis,outcome_dat = telomere,action = 2)
data_h_f<-data_h_f%>%subset(data_h_f$mr_keep==TRUE)

data_h_f_steiger<-steiger_filtering(data_h_f)
data_h_f_steiger<-subset(data_h_f_steiger,steiger_dir==TRUE)

data_h_f_steiger<-data_h_f_steiger[!duplicated(data_h_f_steiger$SNP),]


data_h_f_steiger$R2<-data_h_f_steiger$beta.exposure*data_h_f_steiger$beta.exposure*2*(data_h_f_steiger$eaf.exposure)*(1-data_h_f_steiger$eaf.exposure)
data_h_f_steiger$F<-(data_h_f_steiger$samplesize.exposure-2)*data_h_f_steiger$R2/(1-data_h_f_steiger$R2)

write.csv(data_h_f_steiger[,c("SNP","R2","F")],"data_h_f_steiger_R2_F.csv",quote=F,row.names=F)
write.csv(data_h_f_steiger,"data_h_f_steiger.csv",quote=F,row.names=F)

mr_outcome_f<-mr(data_h_f_steiger)
mr_outcome_f

mr_outcome_f$b_uci95<-mr_outcome_f$b+1.96*mr_outcome_f$se
mr_outcome_f$b_lci95<-mr_outcome_f$b-1.96*mr_outcome_f$se
mr_outcome_f
write.csv(mr_outcome_f[,c("outcome","exposure","method","nsnp","b","se","b_uci95","b_uci95","pval")],
          "data_f_MR.csv",quote=F,row.names=F)


p1<-mr_scatter_plot(mr_outcome_f,data_h_f_steiger)
p1[[1]]


mr_outcome_loo<-mr_leaveoneout(data_h_f_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_h_f_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]


p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]


ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列


heterogeneity_f<-mr_heterogeneity(data_h_f_steiger)
heterogeneity_f

write.csv(heterogeneity_f[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity_f.csv",quote=F,row.names=F)

mr_outcome_f_mre<-mr(data_h_f_steiger,method_list = "mr_ivw_mre")
mr_outcome_f_mre

mr_outcome_f_mre$b_uci95<-mr_outcome_f_mre$b+1.96*mr_outcome_f_mre$se
mr_outcome_f_mre$b_lci95<-mr_outcome_f_mre$b-1.96*mr_outcome_f_mre$se

mr_outcome_f_mre

write.csv(mr_outcome_f_mre[,c("outcome","exposure","method","nsnp","b","b_lci95",
                               "b_uci95","pval")],"data_f_MR_mre.csv",quote=F,row.names=F)

pleiotropy_f <- mr_pleiotropy_test(data_h_f_steiger)
pleiotropy_f

write.csv(pleiotropy_f[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_f.csv",quote=F,row.names=F)

res_MRPRESSO_f<-mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        data=data_h_f_steiger, 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        NbDistribution = 10000,  
                        SignifThreshold = 0.05)
res_MRPRESSO_f

write.csv(res_MRPRESSO_f$`Main MR results`,"MRPRESSO_f_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO_f$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_f_Outlier.csv",quote=F,row.names=F)

write.csv(data_h_f_steiger,"data_h_f_steiger_exoutliers.csv",quote=F,row.names=F)

data_h_f_steiger_exoutliers<-fread("data_h_f_steiger_exoutliers.csv")
head(data_h_f_steiger_exoutliers)

mr_outcome_f_exoutliers<-mr(data_h_f_steiger_exoutliers)
mr_outcome_f_exoutliers

mr_outcome_f_exoutliers$b_uci95<-mr_outcome_f_exoutliers$b+1.96*mr_outcome_f_exoutliers$se
mr_outcome_f_exoutliers$b_lci95<-mr_outcome_f_exoutliers$b-1.96*mr_outcome_f_exoutliers$se

mr_outcome_f_exoutliers

write.csv(mr_outcome_f_exoutliers[,c("outcome","exposure","method","nsnp","b","b_uci95","b_lci95","se","pval")],
          "data_f_exoutliers_MR.csv",quote=F,row.names=F)

p1<-mr_scatter_plot(mr_outcome_f_exoutliers,data_h_f_steiger_exoutliers)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_h_f_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_h_f_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_f_exoutliers<-mr_heterogeneity(data_h_f_steiger_exoutliers)
heterogeneity_f_exoutliers

write.csv(heterogeneity_f_exoutliers[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_f_exoutliers_heterogeneity.csv",quote=F,row.names=F)

pleiotropy_f_exoutliers <- mr_pleiotropy_test(data_h_f_steiger_exoutliers)
pleiotropy_f_exoutliers
write.csv(pleiotropy_f_exoutliers[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_f_exoutliers_pleiotropy.csv",quote=F,row.names=F)

#Same univariable Mendelian Randomization process for GrimAge, HannumAge, intrinsic epigenetic age acceleratio PheoAge, 

#LTL and psoriatic comorbidities##############################################################

#1.LTL and hypertension#######################################################################
#import GWAS summary of telomere length from UKB
telomere<-fread("UKB_telomere_gwas_summarystats.tsv")

telomere<-data.frame(telomere)

head(telomere)

colnames(telomere)<-c("SNP","pval.exposure","chromosome","position.exposure","effect_allele.exposure",
                      "other_allele.exposure","eaf.exposure","beta.exposure","se.exposure")

head(telomere)

telomere$id.exposure<-"Leukocyte telomere length"
telomere$exposure<-"Leukocyte telomere length"
telomere$samplesize.exposure<-464716

head(telomere)

telomere<-subset(telomere,pval.exposure<5e-08)

telomere<-clump_data(telomere,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                     clump_p2 = 5e-08,pop = "EUR")#如果是亚洲人则为“EAS”

hypertension<-fread("finngen_R10_I9_HYPTENS.tsv")
hypertension<-data.frame(hypertension)

head(hypertension)

colnames(hypertension)<-c("chromosome","position.outcome","other_allele.outcome",
                       "effect_allele.outcome","SNP","nearest_genes","pval.outcome","mlogp.outcome",
                       "beta.outcome","se.outcome","eaf.outcome","eafca.outcome","eafcon.outcome")

hypertension$id.outcome<-"Hypertension"
hypertension$outcome<-"Hypertension"
hypertension$samplesize.outcome<-412113

head(hypertension)

total<-merge(telomere,hypertension,by.x = "SNP",by.y = "SNP",all = F)

total<-subset(total,pval.outcome>5e-08)

telomere<-total[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                   "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

hypertension<-total[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                    "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(telomere)
head(hypertension)

#for being palindromic with intermediate allele frequencies:rs2276182, rs2306646, rs4742448, rs56178008, rs6501181, rs670180
data_h<-harmonise_data(exposure_dat = telomere,outcome_dat = hypertension,action = 2)
data_h<-data_h%>%subset(data_h$mr_keep==TRUE)

#Export data-h to check for confounding factors and remove 5（rs10805346、rs2230590、rs2293579、rs429358、rs66731853），n=129
write.csv(data_h,"data_h.csv",quote = F,row.names = F)
#Import and remove mixed SNPs from data_h
data_h<-fread("data_h.csv")

#steiger
data_h_steiger<-steiger_filtering(data_h)
data_h_steiger<-subset(data_h_steiger,steiger_dir==TRUE)

data_h_steiger<-data_h_steiger[!duplicated(data_h_steiger$SNP),]

#F
data_h_steiger$R2<-data_h_steiger$beta.exposure*data_h_steiger$beta.exposure*2*(data_h_steiger$eaf.exposure)*(1-data_h_steiger$eaf.exposure)
data_h_steiger$F<-(data_h_steiger$samplesize.exposure-2)*data_h_steiger$R2/(1-data_h_steiger$R2)

write.csv(data_h_steiger[,c("SNP","R2","F")],"data_steiger_R2_F.csv",quote=F,row.names=F)

#MR
mr_outcome<-mr(data_h_steiger)
mr_outcome

OR <-generate_odds_ratios(mr_outcome)
OR

write.csv(OR[,c("outcome","exposure","method","nsnp","or","or_lci95",
                "or_uci95","pval")],"data_steiger_MR.csv",quote=F,row.names=F)

p1<-mr_scatter_plot(mr_outcome,data_h_steiger)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_h_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_h_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity<-mr_heterogeneity(data_h_steiger)
heterogeneity
write.csv(heterogeneity[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity_steiger.csv",quote=F,row.names=F)

mr_outcome_mre<-mr(data_h_steiger,method_list = "mr_ivw_mre")
mr_outcome_mre

OR_mre <-generate_odds_ratios(mr_outcome_mre)
OR_mre

write.csv(OR_mre[,c("outcome","exposure","method","nsnp","or","or_lci95",
                    "or_uci95","pval")],"data_MR_mre.csv",quote=F,row.names=F)

pleiotropy <- mr_pleiotropy_test(data_h_steiger)
pleiotropy

write.csv(pleiotropy[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_steiger.csv",quote=F,row.names=F)

res_MRPRESSO<-mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        data=data_h_steiger, 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        NbDistribution = 10000,  
                        SignifThreshold = 0.05)
res_MRPRESSO

write.csv(res_MRPRESSO$`Main MR results`,"MRPRESSO_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_Outlier.csv",quote=F,row.names=F)

#Export Steiger data that needs to exclude outliers
write.csv(data_h_steiger,"data_h_steiger_exoutliers.csv",quote=F,row.names=F)

#Import and remove outlier SNPs（n=14,rs113525195、rs115610405、rs13062095、rs17464525、rs1879100、rs1957937、rs2977608、
#rs3865523、rs38664、rs4530278、rs4695407、rs5742915、rs61748181、rs8105767）
#n=115
data_h_steiger_exoutliers<-fread("data_h_steiger_exoutliers.csv")
head(data_h_steiger_exoutliers)

#MR
mr_outcome_exoutliers<-mr(data_h_steiger_exoutliers)
mr_outcome_exoutliers

OR_exoutliers <-generate_odds_ratios(mr_outcome_exoutliers)
OR_exoutliers

write.csv(OR_exoutliers[,c("outcome","exposure","method","nsnp","or","or_lci95",
                "or_uci95","pval")],"data_exoutliers_MR.csv",quote=F,row.names=F)

p1<-mr_scatter_plot(mr_outcome_exoutliers,data_h_steiger_exoutliers)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_h_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_h_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_exoutliers<-mr_heterogeneity(data_h_steiger_exoutliers)
heterogeneity_exoutliers

write.csv(heterogeneity_exoutliers[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_exoutliers_heterogeneity.csv",quote=F,row.names=F)

mr_outcome_mre_exoutliers<-mr(data_h_steiger_exoutliers,method_list = "mr_ivw_mre")
mr_outcome_mre_exoutliers
OR_mre_exoutliers <-generate_odds_ratios(mr_outcome_mre_exoutliers)
OR_mre_exoutliers

write.csv(OR_mre_exoutliers[,c("outcome","exposure","method","nsnp","or","or_lci95",
                    "or_uci95","pval")],"data_MR_mre_exoutliers.csv",quote=F,row.names=F)

pleiotropy_exoutliers <- mr_pleiotropy_test(data_h_steiger_exoutliers)
pleiotropy_exoutliers

write.csv(pleiotropy_exoutliers[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_exoutliers_pleiotropy.csv",quote=F,row.names=F)

#Reverse analysis
hypertension<-fread("finngen_R10_I9_HYPTENS.tsv")
hypertension<-data.frame(hypertension)

head(hypertension)

colnames(hypertension)<-c("chromosome","position.exposure","other_allele.exposure",
                       "effect_allele.exposure","SNP","nearest_genes","pval.exposure","mlogp.exposure",
                       "beta.exposure","se.exposure","eaf.exposure","eafca.exposure","eafcon.exposure")

hypertension$id.exposure<-"Hypertension"
hypertension$exposure<-"Hypertension"
hypertension$samplesize.exposure<-412113

head(hypertension)

hypertension<-subset(hypertension,pval.exposure<5e-08)

hypertension<-clump_data(hypertension,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                      clump_p2 = 5e-08,pop = "EUR")#如果是亚洲人则为“EAS”

telomere<-fread("UKB_telomere_gwas_summarystats.tsv")

telomere<-data.frame(telomere)

head(telomere)

colnames(telomere)<-c("SNP","pval.outcome","chromosome","position.outcome","effect_allele.outcome",
                      "other_allele.outcome","eaf.outcome","beta.outcome","se.outcome")

head(telomere)

telomere$id.outcome<-"Leukocyte telomere length"
telomere$outcome<-"Leukocyte telomere length"
telomere$samplesize.outcome<-464716

head(telomere)

total<-merge(hypertension,telomere,by.x = "SNP",by.y = "SNP",all = F)

total<-subset(total,pval.outcome>5e-08)

hypertension<-total[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                    "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

telomere<-total[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                   "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(hypertension)
head(telomere)

#for being palindromic with intermediate allele frequencies:
#rs10832778, rs1558901, rs1745414, rs2274224, rs2302661, rs28455998, rs4977575, rs7088877, rs9851392
data_h_f<-harmonise_data(exposure_dat = hypertension,outcome_dat = telomere,action = 2)
data_h_f<-data_h%>%subset(data_h$mr_keep==TRUE)

#steiger
data_h_f_steiger<-steiger_filtering(data_h_f)
data_h_f_steiger<-subset(data_h_f_steiger,steiger_dir==TRUE)

data_h_f_steiger<-data_h_f_steiger[!duplicated(data_h_f_steiger$SNP),]

data_h_f_steiger$R2<-data_h_f_steiger$beta.exposure*data_h_f_steiger$beta.exposure*2*(data_h_f_steiger$eaf.exposure)*(1-data_h_f_steiger$eaf.exposure)
data_h_f_steiger$F<-(data_h_f_steiger$samplesize.exposure-2)*data_h_f_steiger$R2/(1-data_h_f_steiger$R2)

write.csv(data_h_f_steiger[,c("SNP","R2","F")],"data_f_steiger_R2_F.csv",quote=F,row.names=F)

#MR
mr_outcome_f<-mr(data_h_f_steiger)
mr_outcome_f

mr_outcome_f$b_uci95<-mr_outcome_f$b+1.96*mr_outcome_f$se
mr_outcome_f$b_lci95<-mr_outcome_f$b-1.96*mr_outcome_f$se
mr_outcome_f
write.csv(mr_outcome_f[,c("outcome","exposure","method","nsnp","b","b_lci95","b_uci95","se","pval")],
          "data_f_steiger_MR.csv",quote=F,row.names=F)

p1<-mr_scatter_plot(mr_outcome,data_h_f_steiger)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_h_f_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_h_f_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_f<-mr_heterogeneity(data_h_f_steiger)
heterogeneity_f

write.csv(heterogeneity_f[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_f_heterogeneity.csv",quote=F,row.names=F)

mr_outcome_mre_f<-mr(data_h_f_steiger,method_list = "mr_ivw_mre")
mr_outcome_mre_f

mr_outcome_mre_f$b_uci95<-mr_outcome_mre_f$b+1.96*mr_outcome_mre_f$se
mr_outcome_mre_f$b_lci95<-mr_outcome_mre_f$b-1.96*mr_outcome_mre_f$se
mr_outcome_mre_f
write.csv(mr_outcome_mre_f[,c("outcome","exposure","method","nsnp","b","b_lci95",
                    "b_uci95","pval")],"data_f_MR_mre.csv",quote=F,row.names=F)

pleiotropy_f <- mr_pleiotropy_test(data_h_f_steiger)
pleiotropy_f
write.csv(pleiotropy_f[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_f.csv",quote=F,row.names=F)

res_MRPRESSO_f<-mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        data=data_h_f_steiger, 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        NbDistribution = 10000,  
                        SignifThreshold = 0.05)
res_MRPRESSO_f

write.csv(res_MRPRESSO_f$`Main MR results`,"MRPRESSO_f_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO_f$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_f_Outlier.csv",quote=F,row.names=F)
write.csv(data_h_f_steiger,"data_f_h_steiger_exoutliers.csv",quote=F,row.names=F)

data_f_h_steiger_exoutliers<-fread("data_f_h_steiger_exoutliers.csv")
head(data_f_h_steiger_exoutliers)

mr_outcome_f<-mr(data_f_h_steiger_exoutliers)
mr_outcome_f

mr_outcome_f$b_uci95<-mr_outcome_f$b+1.96*mr_outcome_f$se
mr_outcome_f$b_lci95<-mr_outcome_f$b-1.96*mr_outcome_f$se
mr_outcome_f
write.csv(mr_outcome_f[,c("outcome","exposure","method","nsnp","b","b_lci95","b_uci95","se","pval")],
          "data_f_outliers_MR.csv",quote=F,row.names=F)

p1<-mr_scatter_plot(mr_outcome_f,data_f_h_steiger_exoutliers)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_f_h_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_f_h_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_f<-mr_heterogeneity(data_f_h_steiger_exoutliers)
heterogeneity_f

write.csv(heterogeneity_f[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_f_outliers_heterogeneity.csv",quote=F,row.names=F)

pleiotropy_f <- mr_pleiotropy_test(data_f_h_steiger_exoutliers)
pleiotropy_f

write.csv(pleiotropy_f[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_f_outliers_pleiotropy.csv",quote=F,row.names=F)

#multivariable MR#########################################################################
#remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(TwoSampleMR)

library(data.table)
#devtools::install_github("shaoming1995/sofmMVIV")
#library(sofmMVIV)
#MVMR：hypertension##########################
psoriasis<-fread("finngen_R10_L12_PSORIASIS.tsv")
psoriasis<-data.frame(psoriasis)

head(psoriasis)

colnames(psoriasis)<-c("chromosome","position.exposure","other_allele.exposure",
                       "effect_allele.exposure","SNP","nearest_genes","pval.exposure","mlogp.exposure",
                       "beta.exposure","se.exposure","eaf.exposure","eafca.exposure","eafcon.exposure")

psoriasis$id.exposure<-"Psoriasis"
psoriasis$exposure<-"Psoriasis"
psoriasis$samplesize.exposure<-407876

head(psoriasis)

EXP1<-subset(psoriasis,pval.exposure<5e-08)

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
EXP1<-clump_data(EXP1,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                      clump_p2 = 5e-08,pop = "EUR")

telomere<-fread("UKB_telomere_gwas_summarystats.tsv")

telomere<-data.frame(telomere)

head(telomere)

colnames(telomere)<-c("SNP","pval.outcome","chromosome","position.outcome","effect_allele.outcome",
                      "other_allele.outcome","eaf.outcome","beta.outcome","se.outcome")

head(telomere)

telomere$id.outcome<-"Leukocyte telomere length"
telomere$outcome<-"Leukocyte telomere length"
telomere$samplesize.outcome<-464716

head(telomere)

total1<-merge(EXP1,telomere,by.x = "SNP",by.y = "SNP",all = F)

total1<-subset(total1,pval.outcome>5e-08)

EXP1<-total1[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                    "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

OUT1<-total1[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                   "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(EXP1)
head(OUT1)

#Removing the following SNPs for being palindromic with intermediate allele frequencies:rs2664280
data_h_f1<-harmonise_data(exposure_dat = EXP1,outcome_dat = OUT1,action = 2)
data_h_f1<-data_h_f1%>%subset(data_h_f1$mr_keep==TRUE)

#steiger
data_h_f_steiger1<-steiger_filtering(data_h_f1)
data_h_f_steiger1<-subset(data_h_f_steiger1,steiger_dir==TRUE)

data_h_f_steiger1<-data_h_f_steiger1[!duplicated(data_h_f_steiger1$SNP),]

hypertension<-fread("finngen_R10_I9_HYPTENS.tsv")
hypertension<-data.frame(hypertension)

head(hypertension)

colnames(hypertension)<-c("chromosome","position.exposure","other_allele.exposure",
                       "effect_allele.exposure","SNP","nearest_genes","pval.exposure","mlogp.exposure",
                       "beta.exposure","se.exposure","eaf.exposure","eafca.exposure","eafcon.exposure")

hypertension$id.exposure<-"Hypertension"
hypertension$exposure<-"Hypertension"
hypertension$samplesize.exposure<-412113

head(hypertension)

EXP2<-subset(hypertension,pval.exposure<5e-08)

EXP2<-clump_data(EXP2,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                      clump_p2 = 5e-08,pop = "EUR")#如果是亚洲人则为“EAS”

total2<-merge(EXP2,telomere,by.x = "SNP",by.y = "SNP",all = F)

total2<-subset(total2,pval.outcome>5e-08)

EXP2<-total2[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
               "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]

OUT2<-total2[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
              "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(EXP2)
head(OUT2)

#Removing the following SNPs for being palindromic with intermediate allele frequencies:
#rs10832778, rs1558901, rs1745414, rs2274224, rs2302661, rs28455998, rs4977575, rs7088877, rs9851392
data_h_f2<-harmonise_data(exposure_dat = EXP2,outcome_dat = OUT2,action = 2)
data_h_f2<-data_h_f2%>%subset(data_h_f2$mr_keep==TRUE)

#steiger过滤
data_h_f_steiger2<-steiger_filtering(data_h_f2)
data_h_f_steiger2<-subset(data_h_f_steiger2,steiger_dir==TRUE)

data_h_f_steiger2<-data_h_f_steiger2[!duplicated(data_h_f_steiger2$SNP),]

exp_name<-c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
            "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")
out_name<-c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
            "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")

#To find the instrumental variables of GWAS SUMMAY for exposure 1 in exposure 2
data_h_steiger_mv1_2<-merge(data_h_f_steiger1[,c("SNP","outcome")],hypertension,by="SNP",all = F)

#To find the instrumental variables of GWAS SUMMAY for exposure 2 in exposure 1
data_h_steiger_mv2_1<-merge(data_h_f_steiger2[,c("SNP","outcome")],psoriasis,by="SNP",all = F)

#merge
MVIV1_1<-data_h_f_steiger1[,exp_name]
MVIV1_2<-data_h_steiger_mv2_1[,exp_name]
MVIV2_2<-data_h_f_steiger2[,exp_name]
MVIV2_1<-data_h_steiger_mv1_2[,exp_name]
exposure_dat<-rbind(MVIV1_1,MVIV1_2,MVIV2_2,MVIV2_1)
write.csv(exposure_dat,"exposure_dat.csv",quote = F,row.names = F)

exposure_dat_rev<-read.csv("exposure_dat.csv",header = T)

head(exposure_dat_rev)

telomere_data<-merge(exposure_dat_rev,telomere,by="SNP",all = F)

telomere_data<-telomere_data[!duplicated(telomere_data$SNP),]

telomere_data<-telomere_data[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                                "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

mvdat<-mv_harmonise_data(exposure_dat_rev,telomere_data)
res<-mv_multiple(mvdat)
write.csv(res,"银屑病与高血压-多变量MR.csv",quote = F,row.names = F)

mv_lasso_res<-mv_subset(mvdat,features = mv_lasso_feature_selection(mvdat),
                        intercept = FALSE,
                        instrument_specific = FALSE,
                        pval_threshold = 5e-08,
                        plots = FALSE)

mv_lasso_res
write.csv(mv_lasso_res,"避免银屑病与高血压共线性-多变量MR.csv",quote = F,row.names = F)

library(MendelianRandomization)
str(mvdat)
MRMVInput<-mr_mvinput(bx=mvdat$exposure_beta,
                    bxse = mvdat$exposure_se,
                    by=mvdat$outcome_beta,
                    byse = mvdat$outcome_se,
                    correlation=matrix())

mvegger<-mr_mvegger(MRMVInput)
mvegger

mvivw<-mr_mvivw(MRMVInput)
mvivw

mv_median<-mr_mvmedian(MRMVInput)
mv_median

mv_lasso<-mr_mvlasso(MRMVInput)
mv_lasso

mv_IVs<-MVMR::strength_mvmr(MRMVInput)
mv_IVs

mv_pleiotropy<-MVMR::pleiotropy_mvmr(MRMVInput)
mv_pleiotropy

#same multivariable MR process for other diseases

#mediation MR-CRP##########################################################################
#C-Reactive protein level Dataset: ieu-b-35（PMID	30388399，Sample size	204,402，Author	Ligthart, S）
#1.X to M, i.e. exposure to psoriasis, outcome to CRP############################################
psoriasis<-fread("finngen_R10_L12_PSORIASIS.tsv")
psoriasis<-data.frame(psoriasis)

head(psoriasis)

colnames(psoriasis)<-c("chromosome","position.exposure","other_allele.exposure",
                       "effect_allele.exposure","SNP","nearest_genes","pval.exposure","mlogp.exposure",
                       "beta.exposure","se.exposure","eaf.exposure","eafca.exposure","eafcon.exposure")

psoriasis$id.exposure<-"Psoriasis"
psoriasis$exposure<-"Psoriasis"
psoriasis$samplesize.exposure<-407876

head(psoriasis)

X<-subset(psoriasis,pval.exposure<5e-08)

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
X<-clump_data(X,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-08,
                         clump_p2 = 5e-08,pop = "EUR")#如果是亚洲人则为“EAS”

head(X)

library(ieugwasr)
#devtools::install_github("MRCIEU/ieugwasr@v0.2.1")
M<-extract_outcome_data(snps = X$SNP,
                          outcomes = "ieu-b-35",
                          proxies = T,
                          maf_threshold = 0.01,
                          access_token = NULL)

M$id.outcome<-"CRP"
M$outcome<-"CRP"
M$samplesize.outcome<-204402
head(M)

M<-M[,c(1:12)]
head(M)

total_XM<-merge(X,M,by.x = "SNP",by.y = "SNP",all = F)
total_XM<-subset(total_XM,pval.outcome>5e-08)

X<-total_XM[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                       "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]
M<-total_XM[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                  "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(X)
head(M)

#Removing the following SNPs for being palindromic with intermediate allele frequencies:rs2664280
data_XM_h<-harmonise_data(exposure_dat = X,outcome_dat = M,action = 2)
data_XM_h<-data_XM_h%>%subset(data_XM_h$mr_keep==TRUE)

write.csv(data_XM_h,"data_XM_h.csv",quote = F,row.names = F)

data_XM_h_steiger<-steiger_filtering(data_XM_h)
data_XM_h_steiger<-subset(data_XM_h_steiger,steiger_dir==TRUE)

data_XM_h_steiger<-data_XM_h_steiger[!duplicated(data_XM_h_steiger$SNP),]

data_XM_h_steiger$R2<-data_XM_h_steiger$beta.exposure*data_XM_h_steiger$beta.exposure*2*(data_XM_h_steiger$eaf.exposure)*(1-data_XM_h_steiger$eaf.exposure)
data_XM_h_steiger$F<-(data_XM_h_steiger$samplesize.exposure-2)*data_XM_h_steiger$R2/(1-data_XM_h_steiger$R2)

write.csv(data_XM_h_steiger[,c("SNP","R2","F")],"data_XM_h_steiger_R2_F.csv",quote=F,row.names=F)
write.csv(data_XM_h_steiger,"data_XM_h_steiger.csv",quote=F,row.names=F)

PS_CRP_mr<-mr(data_XM_h_steiger)
PS_CRP_mr
write.csv(PS_CRP_mr,"PS_CRP_mr.csv",quote = F,row.names = F)

p1<-mr_scatter_plot(PS_CRP_mr,data_XM_h_steiger)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_XM_h_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_XM_h_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity<-mr_heterogeneity(data_XM_h_steiger)
heterogeneity
write.csv(heterogeneity[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity.csv",quote=F,row.names=F)

mr_outcome_mre<-mr(data_XM_h_steiger,method_list = "mr_ivw_mre")
mr_outcome_mre
write.csv(mr_outcome_mre,"data_MR_mre.csv",quote=F,row.names=F)

pleiotropy <- mr_pleiotropy_test(data_XM_h_steiger)
pleiotropy
write.csv(pleiotropy[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy.csv",quote=F,row.names=F)

res_MRPRESSO<-mr_presso(BetaOutcome = "beta.outcome", 
                          BetaExposure = "beta.exposure", 
                          SdOutcome = "se.outcome", 
                          SdExposure = "se.exposure", 
                          data=data_XM_h_steiger, 
                          OUTLIERtest = TRUE, 
                          DISTORTIONtest = TRUE, 
                          NbDistribution = 10000,  
                          SignifThreshold = 0.05)
res_MRPRESSO
write.csv(res_MRPRESSO$`Main MR results`,"MRPRESSO_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_Outlier.csv",quote=F,row.names=F)

write.csv(data_XM_h_steiger,"data_XM_h_steiger_exoutliers.csv",quote=F,row.names=F)

data_XM_h_steiger_exoutliers<-fread("data_XM_h_steiger_exoutliers.csv")
head(data_XM_h_steiger_exoutliers)

PS_CRP_mr_exoutliers<-mr(data_XM_h_steiger_exoutliers)
PS_CRP_mr_exoutliers
write.csv(PS_CRP_mr_exoutliers,"PS_CRP_mr_exoutliers.csv",quote = F,row.names = F)

p1<-mr_scatter_plot(PS_CRP_mr_exoutliers,data_XM_h_steiger_exoutliers)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_XM_h_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_XM_h_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_exoutliers<-mr_heterogeneity(data_XM_h_steiger_exoutliers)
heterogeneity_exoutliers
write.csv(heterogeneity_exoutliers[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity_exoutliers.csv",quote=F,row.names=F)

pleiotropy_exoutliers <- mr_pleiotropy_test(data_XM_h_steiger_exoutliers)
pleiotropy_exoutliers
write.csv(pleiotropy_exoutliers[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_exoutliers.csv",quote=F,row.names=F)

#Reverse-CRP and psoriasis,n=57
M_f<-extract_instruments(outcomes = "ieu-b-35",
                           p1=5e-08,
                           clump = T,
                           r2=0.001,
                           kb=10000,
                           p2=5e-08,
                           access_token = NULL)
M_f$id.exposure<-"CRP"
M_f$exposure<-"CRP"
M_f$samplesize.exposure<-204402

head(M_f)

M_f<-M_f[,c(1:12)]

head(M_f)

psoriasis<-fread("finngen_R10_L12_PSORIASIS.tsv")
psoriasis<-data.frame(psoriasis)

head(psoriasis)

colnames(psoriasis)<-c("chromosome","position.outcome","other_allele.outcome",
                       "effect_allele.outcome","SNP","nearest_genes","pval.outcome","mlogp.outcome",
                       "beta.outcome","se.outcome","eaf.outcome","eafca.outcome","eafcon.outcome")

psoriasis$id.outcome<-"Psoriasis"
psoriasis$outcome<-"Psoriasis"
psoriasis$samplesize.outcome<-407876

head(psoriasis)

total_MX<-merge(M_f,psoriasis,by.x = "SNP",by.y = "SNP",all=F)

total_MX<-subset(total_MX,pval.outcome>5e-08)

M_f<-total_MX[,c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                  "se.exposure","pval.exposure","id.exposure","exposure","samplesize.exposure")]
X_f<-total_MX[,c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome",
                      "se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

head(M_f)
head(X_f)

#for incompatible alleles:rs10778215, rs10778215
#Removing the following SNPs for being palindromic with intermediate allele frequencies:
#rs10778215, rs10778215, rs10778215, rs10778215, rs11108056
data_f_MX_h<-harmonise_data(exposure_dat = M_f,outcome_dat = X_f,action = 2)
data_f_MX_h<-data_f_MX_h%>%subset(data_f_MX_h$mr_keep==TRUE)

write.csv(data_f_MX_h,"data_f_MX_h.csv",quote = F,row.names = F)

data_f_MX_h_steiger<-steiger_filtering(data_f_MX_h)
data_f_MX_h_steiger<-subset(data_f_MX_h_steiger,steiger_dir==TRUE)

data_f_MX_h_steiger<-data_f_MX_h_steiger[!duplicated(data_f_MX_h_steiger$SNP),]

data_f_MX_h_steiger$R2<-data_f_MX_h_steiger$beta.exposure*data_f_MX_h_steiger$beta.exposure*2*(data_f_MX_h_steiger$eaf.exposure)*(1-data_f_MX_h_steiger$eaf.exposure)
data_f_MX_h_steiger$F<-(data_f_MX_h_steiger$samplesize.exposure-2)*data_f_MX_h_steiger$R2/(1-data_f_MX_h_steiger$R2)

write.csv(data_f_MX_h_steiger[,c("SNP","R2","F")],"data_f_MX_h_steiger_R2_F.csv",quote=F,row.names=F)
write.csv(data_f_MX_h_steiger,"data_f_MX_h_steiger.csv",quote=F,row.names=F)

CRP_PS_mr<-mr(data_f_MX_h_steiger)
CRP_PS_mr
write.csv(CRP_PS_mr,"CRP_PS_mr.csv",quote = F,row.names = F)

p1<-mr_scatter_plot(CRP_PS_mr,data_f_MX_h_steiger)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_f_MX_h_steiger)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_f_MX_h_steiger)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

library(ggpubr)
ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_f<-mr_heterogeneity(data_f_MX_h_steiger)
heterogeneity_f
write.csv(heterogeneity_f[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity_f.csv",quote=F,row.names=F)

mr_outcome_mre_f<-mr(data_f_MX_h_steiger,method_list = "mr_ivw_mre")
mr_outcome_mre_f

write.csv(mr_outcome_mre_f[,c("outcome","exposure","method","nsnp","b","se","pval")],
          "data_MR_mre_f.csv",quote=F,row.names=F)

pleiotropy_f <- mr_pleiotropy_test(data_f_MX_h_steiger)
pleiotropy_f
write.csv(pleiotropy_f[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_f.csv",quote=F,row.names=F)

res_MRPRESSO_f<-mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        data=data_f_MX_h_steiger, 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        NbDistribution = 10000,  
                        SignifThreshold = 0.05)
res_MRPRESSO_f
write.csv(res_MRPRESSO_f$`Main MR results`,"MRPRESSO_f_main results.csv",quote=F,row.names=F)
write.csv(res_MRPRESSO_f$`MR-PRESSO results`$`Outlier Test`,"MRPRESSO_f_Outlier.csv",quote=F,row.names=F)
write.csv(data_f_MX_h_steiger,"data_f_MX_h_steiger_exoutliers.csv",quote=F,row.names=F)

data_f_MX_h_steiger_exoutliers<-fread("data_f_MX_h_steiger_exoutliers.csv")
head(data_f_MX_h_steiger_exoutliers)

CRP_PS_mr_f_exoutliers<-mr(data_f_MX_h_steiger_exoutliers)
CRP_PS_mr_f_exoutliers
write.csv(CRP_PS_mr_f_exoutliers,"CRP_PS_mr_exoutliers.csv",quote = F,row.names = F)

p1<-mr_scatter_plot(CRP_PS_mr_f_exoutliers,data_f_MX_h_steiger_exoutliers)
p1[[1]]

mr_outcome_loo<-mr_leaveoneout(data_f_MX_h_steiger_exoutliers)
p3<-mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single<-mr_singlesnp(data_f_MX_h_steiger_exoutliers)
p2<-mr_forest_plot(mr_outcome_single)
p2[[1]]

p4<-mr_funnel_plot(mr_outcome_single)
p4[[1]]

ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],labels = c("A","B","C","D"),nrow = 2,ncol = 2)#2行2列

heterogeneity_f_exoutliers<-mr_heterogeneity(data_f_MX_h_steiger_exoutliers)
heterogeneity_f_exoutliers
write.csv(heterogeneity_f_exoutliers[,c("outcome","exposure","method","Q","Q_df","Q_pval")],
          "data_heterogeneity_f_exoutliers.csv",quote=F,row.names=F)

pleiotropy_f_exoutliers <- mr_pleiotropy_test(data_f_MX_h_steiger_exoutliers)
pleiotropy_f_exoutliers
write.csv(pleiotropy_f_exoutliers[,c("outcome","exposure","egger_intercept","se","pval")],
          "data_pleiotropy_f_exoutliers.csv",quote=F,row.names=F)


#2.M to Y,exposure is CRP, outcome is LTL###############################################
#Reverse：LTL and CRP

#3.MVMR:psoriasis+CRP to CRP###########################################################

#Calculate 95CI% of mediation effect
remotes::install_github("quantPsych/rmediation")
#library(RMediation)
resmed<-medci(mu.x=0.0268764049345607,#a
              mu.y=-0.0218918024941885,#b
              se.x=0.00841134910157432,#a-se
              se.y=0.0106816060507102,#b-se
              rho=0,
              alpha=0.05,
              type="prodclin")
resmed

#Same mediation MR analysis for other mediators