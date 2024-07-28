rm(list = ls())
getwd()
setwd("E:/NHANES/BAPA&pso")

#########################################statistical analysis for nhanes database
#import xlsx
library("readxl")
totalforr <- read_excel("totalforr.xlsx",sheet = "Sheet1")
#Kruskal-Wallis H test##########################################################
#study design
library("survey")
nhanes<-svydesign(data=totalforr,id=~SDMVPSU,strata = ~SDMVSTRA,weights = ~fiweight,nest = TRUE)
#svyranktest
svyranktest(kdm~psonew,nhanes,test = "KruskalWallis")
svyranktest(kdm_advance~psonew,nhanes,test = "KruskalWallis")
svyranktest(phenoage~psonew,nhanes,test = "KruskalWallis")
svyranktest(phenoage_advance~psonew,nhanes,test = "KruskalWallis")

#colour
#https://github.com/BlakeRMills/MetBrewer
#install.packages("MetBrewer")
library("MetBrewer")
# All Palettes
display_all(sequential = FALSE, colorblind_only = FALSE)
isfahan <- MetBrewer::met.brewer("Isfahan1")
vangogh <-MetBrewer::met.brewer("VanGogh3")
benedictus <-MetBrewer::met.brewer("Benedictus")
tsimshian <- MetBrewer::met.brewer("Tsimshian")
length(vangogh)
isfahan[1]
vangogh[3]
tsimshian[4]

#set font
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))

#Forest plot######################################################################
#install.packages("forestploter")
library("grid")
library("forestploter")
library("readxl")
forplot <- read_excel("Forest.xlsx",sheet = "Sheet1")

#NA to blank or NA will be transformed to character
forplot$`Biological age` <- ifelse(is.na(forplot$`Biological age`), "", forplot$`Biological age`)
forplot$`P-trend` <- ifelse(is.na(forplot$`P-trend`), "", forplot$`P-trend`)
forplot$`P-trend2` <- ifelse(is.na(forplot$`P-trend2`), "", forplot$`P-trend2`)
forplot$`P-trend3` <- ifelse(is.na(forplot$`P-trend3`), "", forplot$`P-trend3`)
forplot$se1 <- (log(forplot$hi) - log(forplot$OR))/1.96
forplot$se2 <- (log(forplot$hi2) - log(forplot$OR2))/1.96
forplot$se3 <- (log(forplot$hi3) - log(forplot$OR3))/1.96

# Add blank column for the forest plot to display CI.
# Adjust the column width with space. 
forplot$` ` <- paste(rep(" ", 20), collapse = " ")

# Create confidence interval column to display
forplot$`OR (95% CI)` <- ifelse(is.na(forplot$se1), "",
                                sprintf("%.2f (%.2f, %.2f)",
                                        forplot$OR, forplot$low, forplot$hi))
forplot$`OR (95% CI)2` <- ifelse(is.na(forplot$se2), "",
                                 sprintf("%.2f (%.2f, %.2f)",
                                         forplot$OR2, forplot$low2, forplot$hi2))
forplot$`OR (95% CI)3` <- ifelse(is.na(forplot$se3), "",
                                 sprintf("%.2f (%.2f, %.2f)",
                                         forplot$OR3, forplot$low3, forplot$hi3))

attributes(forplot)
p <- forest(forplot[,c(1:2,18,6)],
            est = forplot$OR,
            lower = forplot$low, 
            upper = forplot$hi,
            ci_column = 3,
            xlim = c(0, 4.5),
            ticks_at = c(0.5, 2, 3, 4),
            ref_line = 1)

# Print plot
plot(p)

#set font
windowsFonts(myFont1 = windowsFont("Times New Roman"))

tm <- forest_theme(base_size = 9,
                   base_family = "myFont1",
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = benedictus[3],
                   ci_fill = benedictus[3],
                   ci_alpha = 0.8,  #
                   ci_lty = 1,# 可信区间的线型
                   ci_lwd = 2,# 可信区间的线宽
                   ci_Theight = 0.2, # Set an T end at the end of CI
                   # 设置T字在可信区间末端的高度，默认是NULL 
                   # Reference line width/type/color  # 设置参考线的外观
                   refline_lwd = 2,# 参考线的线宽
                   refline_lty = "dashed",# 参考线的线型
                   refline_col = "grey20",# 参考线的颜色
                   # Vertical line width/type/color # 设置垂直线的外观
                   vertline_lwd = 1,# 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lty = "dashed",# 垂直线的线型
                   vertline_col = "grey20")#垂直线的颜色
                   # 设置脚注的字体大小、字体样式和颜色
                   #footnote_cex = 0.6,            # 脚注字体大小
                   #footnote_fontface = "italic",  # 脚注字体样式
                   #footnote_col = "red4"          # 脚注文本的颜色）

pt <- forest(forplot[,c(1:2,18,6)],
             est = forplot$OR,
             lower = forplot$low, 
             upper = forplot$hi,
             is_summary = c(rep(FALSE, nrow(forplot)-1), FALSE),
             ci_column = 3,
             ref_line = 1,
             xlim = c(0, 4.5),
             colgap = unit(1,"mm"),
             ticks_at = c(0.5, 1, 2, 3, 4),
             theme = tm)
plot(pt)

# Add blank column for the second CI column
forplot$`   ` <- paste(rep(" ", 20), collapse = " ")

attributes(forplot)

p <- forest(forplot[,c(1:2,18,6,22,10)],
            est = list(forplot$OR,
                       forplot$OR2),
            lower = list(forplot$low,
                         forplot$low2), 
            upper = list(forplot$hi,
                         forplot$hi2),
            ci_column = c(3, 5),
            ref_line = c(1,1),
            nudge_y = 0.2,
            x_trans = c("none", "none"),
            xlim = list(c(0, 4), c(0, 4)),
            graphwidth = unit(80,"mm"),
            colgap = unit(0.5,"mm"),
            ticks_at = list(c(0.5, 1, 2, 3, 4), c(0.5, 1, 2, 3, 4)),
            theme = tm)

plot(p)

attributes(forplot)

# Add blank column for the second CI column
forplot$`     ` <- paste(rep(" ", 20), collapse = " ")

p <- forest(forplot[,c(1:2,18,6,22,10,23,14)],
            est = list(forplot$OR,
                       forplot$OR2,
                       forplot$OR3),
            lower = list(forplot$low,
                         forplot$low2,
                         forplot$low3), 
            upper = list(forplot$hi,
                         forplot$hi2,
                         forplot$hi3),
            ci_column = c(3, 5, 7),
            ref_line = c(1,1,1),
            nudge_y = 0.2,
            x_trans = c("none", "none","none"),
            xlim = list(c(0, 4), c(0, 4),c(0, 3.5)),
            graphwidth = unit(120,"mm"),
            colgap = unit(0.5,"mm"),
            ticks_at = list(c(0.5, 1, 2, 3, 4), c(0.5, 1, 2, 3, 4), c(0.5, 1, 2, 3, 4)),
            theme = tm)

plot(p)

# 把第十行变成红色
p <- edit_plot(p, row = 19, gp = gpar(col = "red4"))
p <- edit_plot(p, row = 24, gp = gpar(col = "red4"))

# 把全部的文本变成粗体
p <- edit_plot(p,
               row = 1:25,
               gp = gpar(fontface = "bold"))
plot(p)

#add border
p<-add_border(p,
               row = 1,
               col = NULL,
               part = c("header"),
               where = c("bottom"),
               gp = gpar(lwd = 2))

p<-add_border(p,
               row = 1,
               col = NULL,
               part = c("header"),
               where = c("top"),
               gp = gpar(lwd = 2))

p<-add_border(p,
               row = 25,
               col = NULL,
               part = c("header"),
               where = c("bottom", "left", "top", "right"),
               gp = gpar(lwd = 2))

# Print plot
plot(p)

#import xlsx
library("readxl")
totalforcor <- read_excel("totalforcor.xlsx",sheet = "Sheet1")

#install.packages("GGally")
#install.packages("rstatix")

library(GGally)
library(rstatix)

vig_ggally("ggpairs")
?ggpairs

library(ggplot2)
library(dplyr)
library(tidyr)

#set font
windowsFonts(myFont1 = windowsFont("Times New Roman"))

plot_ba = function(totalforcor, agevar, label) {
  
  cor_label = totalforcor %>%
    pivot_longer(all_of(agevar), names_to = "method", values_to = "measure") %>%
    mutate(method = factor(method, levels = agevar, labels = label)) %>%
    group_by(method) %>%
    summarise(cor(measure, age, use="complete.obs"))
  
  colnames(cor_label) = c("method","r")
  cor_label$r = round(cor_label$r, 3)
  cor_label$r = paste0("r = ", cor_label$r)
  
  plot = totalforcor %>%
    pivot_longer(all_of(agevar), names_to = "method", values_to = "measure") %>%
    mutate(method = factor(method, levels = agevar, labels = label)) %>%
    ggplot(., aes(x = age,y = measure, group = method)) +
    geom_point(shape = 1,colour=benedictus[3]) +
    geom_smooth(method = lm, se=TRUE, color = "black", linetype = 2, fullrange=T, size = 1) +
    facet_wrap(~ method, scales = "free") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    scale_x_continuous(name ="Chronological Age") +
    scale_y_continuous(name ="Biological Age") +
    geom_label(data = cor_label, aes(label = r),
               x = Inf, y = -Inf, hjust=1, vjust=0,
               label.size=0, inherit.aes = FALSE, size = 6,alpha = 1, family="myFont1")+
    theme_bw() +
    theme(axis.text = element_text(face = "bold", family = "myFont1",size = 14,colour = "black"), 
          axis.text.x = element_text(face = "bold", family = "myFont1",size = 14,colour = "black"), 
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold",family="myFont1"), 
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size=14, face = "bold",family="myFont1"), 
          legend.text = element_text(size = 13,face = "bold", family = "myFont1",),
          strip.background = element_rect(fill = "white"), 
          strip.text = element_text(size=14,face = "bold",family="myFont1"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(plot)
  
  
}
#select biological age variables
agevar = c("KDM-BA","PhenoAge")
#prepare labels
label = c("KDM-BA",
          "PhenoAge")
#plot age vs bioage
plot_ba(totalforcor, agevar, label)

library(dplyr)
advance <- totalforr %>% select(kdm_advance,phenoage_advance)

#heatmap########################################################################
library("ggplot2")
library("ggpubr")
p1<-ggplot(advance, aes(kdm_advance, phenoage_advance)) + 
  scale_x_continuous(name ="KDM-BA acceleration") +
  scale_y_continuous(name ="PhenoAge acceleration") +
  geom_point(shape = 1,colour=benedictus[3])+ 
  geom_smooth(method="lm",se=TRUE, color = "black", linetype = 2, fullrange=T, linewidth = 1,formula = y ~ x) + 
  theme_bw()+
  theme(axis.text = element_text(face = "bold", family = "myFont1",size = 14,colour = "black"), 
        axis.text.x = element_text(face = "bold", family = "myFont1",size = 14,colour = "black"), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold",family="myFont1"), 
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=14, face = "bold",family="myFont1"), 
        legend.text = element_text(size = 13,face = "bold", family = "myFont1",),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=14,face = "bold",family="myFont1"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()+
        stat_cor(method = 'spearman', aes(x = kdm_advance, y = phenoage_advance))+
        annotate('text',label='r = 0.227',x = 110,y = 20,
                 size=6,color='black'))
p1

#correlation
r <- cor.test(advance$kdm_advance,advance$phenoage_advance,method = 'spearman')
r$estimate
r$p.value

#import xlsx
library("readxl")
totalmi<- read_excel("totalmi_zscore.xlsx",sheet = "Sheet1")

str(totalmi)
attributes(totalmi)

totalmi$hyper<-as.factor(totalmi$hyper)
totalmi$CHD<-as.factor(totalmi$CHD)
totalmi$Hcho<-as.factor(totalmi$Hcho)
totalmi$PAD700<-as.factor(totalmi$PAD700)
totalmi$DM<-as.factor(totalmi$DM)
totalmi$CKD<-as.factor(totalmi$CKD)
totalmi$psonew<-as.factor(totalmi$psonew)
totalmi$cancer<-as.factor(totalmi$cancer)
totalmi$agec<-as.factor(totalmi$agec)
totalmi$bmicnew<-as.factor(totalmi$bmicnew)
totalmi$educnew<-as.factor(totalmi$educnew)
totalmi$PIRcnew<-as.factor(totalmi$PIRcnew)
totalmi$smknew<-as.factor(totalmi$smknew)
totalmi$alqnew<-as.factor(totalmi$alqnew)
totalmi$RIAGENDR<-as.factor(totalmi$RIAGENDR)
totalmi$RIDRETH1<-as.factor(totalmi$RIDRETH1)

#RCS，knot##################################################################
library("splines")
library("rms")
library("dplyr")

imp1<-subset(totalmi, Imputation =="10")

#Modeling, set df=7, 8, 9 respectively, construct 3 models and compare them, with node number=df-4
md1 = lrm(psonew ~ rcs(phenoage_advance , 3), data = imp1)
md2 = lrm(psonew ~ rcs(phenoage_advance , 4), data = imp1)
md3 = lrm(psonew ~ rcs(phenoage_advance , 5), data = imp1)
#AIC
AIC(md1)
AIC(md2)
AIC(md3)
#BIC
BIC(md1)
BIC(md2)
BIC(md3)

#FDR#########################################################################
#P-trend adjusted for FDR
#crude model
crude <- c(0.0008,0.1229,0.0001,1.44045E-05)
round(p.adjust(crude, "BH"), 4)

#model1
model1 <- c(0.1255,0.2768,0.0314,0.0007)
round(p.adjust(model1, "BH"), 4)

#model2
model2 <- c(0.3369,0.3398,0.0993,0.0008)
round(p.adjust(model2, "BH"), 4)

#RCS FDR adjusted P-value
Poverall<-c(0.62,0.24,0.39,0.0011)
round(p.adjust(Poverall, "BH"), 4)

Pnonlinear<-c(0.76,0.23,0.91,0.62)
round(p.adjust(Pnonlinear, "BH"), 4)

#MR analysis P-value adjusted for FDR
IVW <- c(0.030270454,0.674513162609503,0.336246236488667,0.348411948910852,0.522061322012991)
round(p.adjust(IVW, "BH"), 5)

Weightedmode <- c(0.017954716,0.814409637102852,0.558330541042247,0.484174319231243,0.794189626386805)
round(p.adjust(Weightedmode, "BH"), 5)

Weightedmedian<- c(0.006331004,0.68026821704495,0.572535533933234,0.348411948910852,0.961559353790386)
round(p.adjust(Weightedmedian, "BH"), 5)

egger<- c(0.40367795,0.981798546939564,0.265835260638123,0.520527414104841,0.850820829294493)
round(p.adjust(egger, "BH"), 5)

PRESSO<- c(0.038075393552883,0.6400623,0.3253643,0.4064977,0.4668151)
round(p.adjust(PRESSO, "BH"), 5)

continous<-c(0.3324,0.3375,0.0926,0.0004)
round(p.adjust(continous, "BH"), 4)

IVW <- c(0.033,0.610,0.106,0.001,0.902)
round(p.adjust(IVW, "BH"), 5)

IVW_mre <- c(0.033,0.106,0.902)
round(p.adjust(IVW_mre, "BH"), 3)

IVW_exoutliers <- c(0.060,0.002,0.612)
round(p.adjust(IVW_exoutliers, "BH"), 3)

Weightedmode <-c(0.078,0.850,0.143,0.179,0.433)
round(p.adjust(Weightedmode, "BH"), 5)

Weightedmedian <-c(0.143,0.944,0.042,0.006,0.463)
round(p.adjust(Weightedmedian, "BH"), 5)

egger<- c(0.601,0.312,0.366,0.737,0.064)
round(p.adjust(egger, "BH"), 5)

PRESSO<- c(0.035,0.615,0.145,0.904)
round(p.adjust(PRESSO, "BH"), 4)

#weighted RCS##########################################################################
#import xlsx
library("readxl")
totalforrcs <- read_excel("totalforrrcs.xlsx",sheet = "Sheet1")
str(totalforrcs)
totalforrcs$RIAGENDR<-as.factor(totalforrcs$RIAGENDR)
totalforrcs$RIDRETH1<-as.factor(totalforrcs$RIDRETH1)
totalforrcs$RIDRETH1<-as.factor(totalforrcs$RIDRETH1)
totalforrcs$agec<-as.factor(totalforrcs$agec)
totalforrcs$bmicnew<-as.factor(totalforrcs$bmicnew)
totalforrcs$educnew<-as.factor(totalforrcs$educnew)
totalforrcs$smknew<-as.factor(totalforrcs$smknew)
totalforrcs$alqnew<-as.factor(totalforrcs$alqnew)
totalforrcs$PIRcnew<-as.factor(totalforrcs$PIRcnew)
totalforrcs$hyper<-as.factor(totalforrcs$hyper)
totalforrcs$DM<-as.factor(totalforrcs$DM)
totalforrcs$CHD<-as.factor(totalforrcs$CHD)
totalforrcs$CKD<-as.factor(totalforrcs$CKD)
totalforrcs$Hcho<-as.factor(totalforrcs$Hcho)
totalforrcs$PAD700<-as.factor(totalforrcs$PAD700)
totalforrcs$yearcycle<-as.factor(totalforrcs$yearcycle)
totalforrcs$psonew<-as.factor(totalforrcs$psonew)
summary(totalforrcs)

#first imputated database
totalforrcs_1<-subset(totalforrcs,totalforrcs$`_Imputation_`==1)
str(totalforrcs_1)
#exposure:kdm, kdm_advance, phenoage, phenoage_advance
#outcome:psonew
#covariates:
#model1:RIAGENDR, agec, RIDRETH1, educnew, PIRcnew, bmicnew, PAD700, alqnew, smknew.
#model2:mode1+hyper, DM, CKD, Hcho, CHD.

##analysis_or_1
analysis_OR_1<-data.frame(
  kdm=totalforrcs_1$kdm,
  RIAGENDR=1,
  agec=1,
  RIDRETH1=3,
  educnew=2,
  PIRcnew=3,
  bmicnew=1,
  PAD700=2,
  alqnew=1,
  smknew=1,
  hyper=0,
  DM=0,
  CKD=0,
  Hcho=0,
  CHD=0)

##as.factor
str(analysis_OR_1)
analysis_OR_1$RIAGENDR <- as.factor(analysis_OR_1$RIAGENDR)#1
analysis_OR_1$agec <- as.factor(analysis_OR_1$agec)#2
analysis_OR_1$RIDRETH1 <- as.factor(analysis_OR_1$RIDRETH1)#3
analysis_OR_1$educnew <- as.factor(analysis_OR_1$educnew)#4
analysis_OR_1$PIRcnew <- as.factor(analysis_OR_1$PIRcnew)#5
analysis_OR_1$bmicnew <- as.factor(analysis_OR_1$bmicnew)#6
analysis_OR_1$PAD700 <- as.factor(analysis_OR_1$PAD700)#7
analysis_OR_1$alqnew <- as.factor(analysis_OR_1$alqnew)#8
analysis_OR_1$smknew <- as.factor(analysis_OR_1$smknew)#9
analysis_OR_1$hyper <- as.factor(analysis_OR_1$hyper)#10
analysis_OR_1$DM <- as.factor(analysis_OR_1$DM)#11
analysis_OR_1$CKD <- as.factor(analysis_OR_1$CKD)#12
analysis_OR_1$Hcho <- as.factor(analysis_OR_1$Hcho)#13
analysis_OR_1$CHD <- as.factor(analysis_OR_1$CHD)#14
str(analysis_OR_1)

#design
library("survey")
NHANES_survey<-svydesign(data=totalforrcs_1,id=~SDMVPSU,strata = ~SDMVSTRA,weights = ~fiweight,nest = TRUE)
#frequency of kdm-ba

library('rms')
#knots=3，P10，P50, P90
#knots=4,P5,P35,P65,P95
#knots=5,P5,P27.5,P50,P72.5,P95

#calculation for weighted frequency
weighted_quantile90 <- quantile(totalforrcs_1$kdm, probs = 0.9, type = 6, weights = totalforrcs_1$fiweight)
weighted_quantile90

weighted_quantile50 <- quantile(totalforrcs_1$kdm, probs = 0.5, type = 6, weights = totalforrcs_1$fiweight)
weighted_quantile50

weighted_quantile10 <- quantile(totalforrcs_1$kdm, probs = 0.1, type = 6, weights = totalforrcs_1$fiweight)
weighted_quantile10



svy.fit <- svyglm(psonew~rcs(kdm,c(23.03507,40.73952,68.08127))+RIAGENDR+agec+RIDRETH1+
                  educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD,
                  design=NHANES_survey,data=totalforrcs_1,family= quasibinomial)

##Calculate the linear prediction value of the reference value, where the variable lnr is called link, not fit
lnr<- as.data.frame(predict(svy.fit,se= TRUE,newdata=analysis_OR_1,ref.zero = TRUE)) 
lnr_ref<- lnr$link-lnr$link[(which(analysis_OR_1$kdm==median(analysis_OR_1$kdm)))[1]] 
OR_ref<-exp(lnr_ref)
analysis_OR_1$OR_ref<-OR_ref

library('ggplot2')
##Use the calculated redrawing to set reference values
p1<-ggplot()+theme_classic()+
  geom_line(data=analysis_OR_1, aes(kdm,OR_ref),linetype="solid",linewidth=1,alpha = 1)+
  scale_x_continuous(limits=c(-4,16),breaks=seq(-4,16,4))+   
  scale_y_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.3))+  
  geom_hline(yintercept=1,linetype=2,linewidth=0.25)
p1

library(boot)
##bootstrap
bs <- function(formula, data, indices) {
  d <- data[indices,]
  survey <- svydesign(data=d, id=~SDMVPSU, strata=~SDMVSTRA, weights=~fiweight, nest=TRUE)
  svy.fit <- svyglm(formula,design=survey,data=d,family= quasibinomial)
  lnr<- predict(svy.fit,newdata=analysis_OR_1,ref.zero = TRUE)
  lnr_ref<- lnr-lnr[(which(analysis_OR_1$kdm==median(analysis_OR_1$kdm)))[1]] #利用which找到中位数对应的位置
  return(lnr_ref)
}

set.seed(1234)
results <- boot(data=totalforrcs_1, statistic=bs,
                R=1000, formula=psonew~rcs(kdm,c(23.03507,40.73952,68.08127))+RIAGENDR+agec+RIDRETH1+
                educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD)
print(results)


##Bootstrap uses percentiles for function calculations
OR_ci<- function(Index){
  ORR<- boot.ci(results, type="perc",index = Index)
  out<- data.frame(
    low=as.numeric(ORR$perc)[4],
    high=as.numeric(ORR$perc)[5]
  )
  return(out)
}
library("purrr")
##Calculate the corresponding CI for each observation, obtain the corresponding number of observations through nrow, and correspond to the number of observations in the original dataset
ci<- map_dfr(1:nrow(analysis_OR_1),OR_ci)
ci

##Add the estimated value and confidence interval to analysis_oR_1
analysis_OR_1$OR_est<-exp(lnr_ref)
analysis_OR_1$lowlimit<-exp(ci$low)
analysis_OR_1$uplimit<-exp(ci$high)

##plot
p1<-ggplot()+theme_classic()+
  geom_line(data=analysis_OR_1, aes(kdm,OR_est),linetype="solid",linewidth=1,alpha = 0.9,colour="black")+
  geom_ribbon(data=analysis_OR_1, aes(kdm,ymin = lowlimit, ymax = uplimit),alpha = 0.2,fill="red")+
  scale_x_continuous(limits=c(0,155),breaks=seq(0,155,30))+   
  scale_y_continuous(limits=c(0,13),breaks=seq(0,13,3),name = "OR (95%CI) of psoriasis")+
  geom_hline(yintercept=1,linetype=2,linewidth=0.6)+
  xlab("KDM-BA (years)")+
  theme(axis.title = element_text(size = 14, color = "black", face = "bold",family = "A"),
        axis.text = element_text(size=14, color = "black", face = "bold",family = "A"),
        panel.grid = element_blank())
p1

#P-overall,P-nonlinear#######################################################

#exposure:kdm, kdm_advance, phenoage, phenoage_advance
#outcome:psonew
#covariates:
#model1:RIAGENDR, agec, RIDRETH1, educnew, PIRcnew, bmicnew, PAD700, alqnew, smknew.
#model2:mode1+hyper, DM, CKD, Hcho, CHD.

#summary
svy.fit <- svyglm(psonew~rcs(kdm_advance,3)+RIAGENDR+agec+RIDRETH1+
                    educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD,
                  design=NHANES_survey,data=totalforrcs_1,family= quasibinomial)
summary(svy.fit)

svy.fit1 <- svyglm(psonew~rcs(kdm_advance,4)+RIAGENDR+agec+RIDRETH1+
                     educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD,
                   design=NHANES_survey,data=totalforrcs_1,family= quasibinomial)
summary(svy.fit1)

svy.fit2 <- svyglm(psonew~rcs(kdm_advance,5)+RIAGENDR+agec+RIDRETH1+
                     educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD,
                   design=NHANES_survey,data=totalforrcs_1,family= quasibinomial)
summary(svy.fit2)

#svy.fit
svy.fit1 <- svyglm(psonew~rcs(phenoage_advance,c(-1.08121,-0.1216926,1.197989))+RIAGENDR+agec+RIDRETH1+
                   educnew+PIRcnew+bmicnew+PAD700+alqnew+smknew+hyper+DM+CKD+Hcho+CHD,
                   design=NHANES_survey,data=totalforrcs_1,family= quasibinomial)

summary(svy.fit1)

#install.packages("aod")
library('aod')
#total
wald.test(Sigma = vcov(svy.fit1), b = coef(svy.fit1), Terms = 2:3)
#nonlinear
wald.test(Sigma = vcov(svy.fit1), b = coef(svy.fit1), Terms = 3)









