libname BAPA 'E:\NHANES\BAPA&pso';

/**************************************import 2003-2004 cycles data********************************************/
proc import out=demo_2003/*general information*/
datafile="E:\NHANES\dataset\2003\demo_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data demo_2003;
set demo_2003;
keep SEQN RIAGENDR RIDAGEYR RIDRETH1 DMDEDUC2 INDFMPIR WTMEC2YR SDMVPSU SDMVSTRA yearcycle;
run;

proc sort data=demo_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=bmi_2003/*body mass index*/
datafile="E:\NHANES\dataset\2003\bmi_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data bmi_2003;
set bmi_2003;
keep SEQN BMXBMI;
run;

proc sort data=bmi_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=cotinine_2003/*cotinine*/
datafile="E:\NHANES\dataset\2003\cotinine_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data cotinine_2003;
set cotinine_2003;
keep SEQN LBXCOT;
run;

proc sort data=cotinine_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=alq_2003/*alcoh using*/
datafile="E:\NHANES\dataset\2003\alq_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data alq_2003;
set alq_2003;
keep SEQN ALQ180;
run;

proc sort data=alq_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=deq_2003/*psoriasis*/
datafile="E:\NHANES\dataset\2003\deq_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data deq_2003;
set deq_2003;
keep SEQN DEQ053;
run;

proc sort data=deq_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=BAfem_2003/*KDM-BA for females*/
datafile="E:\NHANES\dataset\2003\BA&PA_2003\results\kdm_fem2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data BAfem_2003;
set BAfem_2003;
keep SEQN kdm kdm_advance;
run;

proc sort data=BAfem_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=BAmale_2003/*KDM-BA for males*/
datafile="E:\NHANES\dataset\2003\BA&PA_2003\results\kdm_male2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data BAmale_2003;
set BAmale_2003;
keep SEQN kdm kdm_advance;
run;

proc sort data=BAmale_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=PA_2003/*PhenoAge*/
datafile="E:\NHANES\dataset\2003\BA&PA_2003\results\PA_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data PA_2003;
set PA_2003;
keep SEQN phenoage phenoage_advance;
run;

proc sort data=PA_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=bp_2003/*hypertension*/
datafile="E:\NHANES\dataset\2003\bp_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data bp_2003;
set bp_2003;
keep SEQN hyper;
run;

proc sort data=bp_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=CHD_2003/*coronary heart disease*/
datafile="E:\NHANES\dataset\2003\CHD_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data CHD_2003;
set CHD_2003;
keep SEQN CHD;
run;

proc sort data=CHD_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=DM_2003/*type 2 diabetes*/
datafile="E:\NHANES\dataset\2003\dm2_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data DM_2003;
set DM_2003;
keep SEQN DM;
run;

proc sort data=DM_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=CKD_2003/*chronic kidney disease*/
datafile="E:\NHANES\dataset\2003\ckd_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data CKD_2003;
set CKD_2003;
keep SEQN CKD;
run;

proc sort data=CKD_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=lipid_2003/*hypercholesterolemia*/
datafile="E:\NHANES\dataset\2003\lipid_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data lipid_2003;
set lipid_2003;
keep SEQN Hcho;
run;

proc sort data=lipid_2003;/*sorting by SEQN*/
by SEQN;
run;

proc import out=paq_2003/*physical activity*/
datafile="E:\NHANES\dataset\2003\paq_2003.xlsx"
dbms=xlsx replace;
sheet="sheet1";
getnames=yes;
run;

data paq_2003;
set paq_2003;
keep SEQN PAD700;
run;

proc sort data=paq_2003;/*sorting by SEQN*/
by SEQN;
run;

/*merge database*/
data total_2003;
merge demo_2003 bmi_2003 cotinine_2003 alq_2003 deq_2003 BAfem_2003 BAmale_2003 PA_2003 bp_2003 CHD_2003 DM_2003 CKD_2003 lipid_2003
paq_2003 mcq_2003;
by SEQN;
run;

data totalduplicate;/*check duplication*/
set total_2003;
run;

proc sort nouniquekey out=totalduplicate;
by SEQN;
proc print;
run;
/*no duplication，N=10122*/

/*rename variables*/
data total_2003;
set total_2003;
rename ALQ180=alcoh;
rename DEQ053=pso;
rename MCQ220=cancer;
run;

/*recalculate sample weight,fiweight=1/3 of sample weight*/
data total_2003;
set total_2003;
fiweight=WTMEC2YR/3;
run;

/**************************************same importion processing for 2005-2006, and 2009-2010 cycles********************************************/

/*merge*/
data BAPA.total;/*n=31007*/
merge total_2003 total_2005 total_2009;
by SEQN;
run;

/*exclude psoriasis missing data,n=13418*/
data total;
set BAPA.total;
if pso^=. and pso^=7 and pso^=9;
run;

/*exclude subjects younger than 20 years old,n=13418*/
data total;
set total;
if RIDAGEYR>=20;
run;

/*exclude biological aging biomarkers missing data,n=10912*/
data total;
set total;
if kdm^=. and phenoage^=.;
run;

/*exclude psoriatic comorbidities missing data，n=10883*/
data total;
set total;
if hyper^=. and DM^=. and Hcho^=. and CKD^=. and CHD^=.;
run;

proc freq data=total;
tables RIAGENDR RIDAGEYR RIDRETH1 DMDEDUC2 INDFMPIR SDMVPSU SDMVSTRA fiweight yearcycle pso
alcoh LBXCOT BMXBMI kdm kdm_advance phenoage phenoage_advance;
run;

/*变量nmiss*//*ntotal=10883*/
/*RIAGENDR RIDAGEYR RIDRETH1 DMDEDUC2(7or9)=13 INDFMPIR=720 SDMVPSU SDMVSTRA fiweight yearcycle*/
/*alcoh=812*/
/*LBXCOT=3*/
/*BMXBMI=84*/

proc freq data=total;
tables hyper DM CHD Hcho CKD PAD700;
run;

data total;
set total;
if pso=1 then psonew=1;
if pso=2 then psonew=0;
run;
 
/*set new variables 
agec:20-60，>60；bmic：normal、overweight、obesity；educ:lower than high school、 high school、higher than high school，
cotinine，<0.011 ng/mL=none，0.011-10=slight，≥10 ng/mL为severe，
INDFMPIR represented in low (PIR≤1.3), middle (1.3-3.5), high (PIR≥3.5)*/
data total;
set total;
if RIDAGEYR^=. then do;
if RIDAGEYR<60 and RIDAGEYR>=20 then agec=1;
if RIDAGEYR>=60 then agec=2;
end; 

if BMXBMI^=. then do;
if BMXBMI<25 and BMXBMI>=0 then bmic=1;
if BMXBMI>=25 and BMXBMI<30 then bmic=2;
if BMXBMI>=30 then bmic=3;
end;

if DMDEDUC2^=. then do;
if DMDEDUC2<=2 then educ=1;
if DMDEDUC2=3 then educ=2;
if DMDEDUC2>3 and DMDEDUC2<7 then educ=3;
if DMDEDUC2>=7 then educ=.;
end;

if LBXCOT^=. then do;
if LBXCOT=0.011 then smk=1;
if LBXCOT>0.011 and LBXCOT<10 then smk=2;
if LBXCOT>=10 then smk=3;
end;

if INDFMPIR^=. then do;
if INDFMPIR<1.3 then PIRc=1;
if INDFMPIR>=1.3 and INDFMPIR<3.5 then PIRc=2;
if INDFMPIR>=3.5 then PIRc=3;
end; 

run;

proc freq data=total;
tables agec bmic educ smk PIRc hyper CHD Hcho PAD700;
run;

/*correlation between KDM-BA, PhenoAge and chronological age*/
proc corr data=total pearson spearman plots=matrix(histogram) plots(maxpoints=100000);
var RIDAGEYR kdm phenoage;
run;

proc surveyfreq data=total;
tables RIAGENDR agec RIDRETH1 bmic educ PIRc smk alcoh hyper DM CKD Hcho CHD PAD700;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
where psonew=0;
run;

proc univariate data=total normal plot;
var kdm kdm_advance phenoage phenoage_advance;
run;

proc surveymeans data=total mean stderr;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
var kdm kdm_advance phenoage phenoage_advance;
where psonew=0;
run;

proc surveyreg data=total;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
model phenoage_advance=psonew/vadjust=none anova solution clparm;
run;

/*X2 test*/
/*variables：RIAGENDR agec RIDRETH1 bmic educ PIRc smk alcoh hyper DM CKD Hcho CHD PAD700*/
proc surveyfreq data=total;
tables CHD*psonew/row nowt chisq wchisq;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
run;

/*Normalization for kdm_advance and phenoage_advance，mean=0，SD=1*/
proc standard data=total out=total_zscore mean=0 std=1;
var kdm_advance phenoage_advance;
run;
proc print data=total_zscore(obs=6);
run;

data total_zscore;/*tertiles*/
set total_zscore;

if kdm^=. then do;
if kdm<33.772901549 then kdmT=1;
if kdm>=33.772901549 and kdm<47.250336148 then kdmT=2;
if kdm>=47.250336148 then kdmT=3;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.437533583 then kdm_advanceT=1;
if kdm_advance>=-0.437533583 and kdm_advance<0.240498593 then kdm_advanceT=2;
if kdm_advance>=0.240498593 then kdm_advanceT=3;
end;

if phenoage^=. then do;
if phenoage<31.918910852 then phenoageT=1;
if phenoage>=31.918910852 and phenoage<46.355834752 then phenoageT=2;
if phenoage>=46.355834752 then phenoageT=3;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.509356109 then phenoage_advanceT=1;
if phenoage_advance>=-0.509356109 and phenoage_advance<0.1039726113 then phenoage_advanceT=2;
if phenoage_advance>=0.1039726113 then phenoage_advanceT=3;
end;

run;

data total_zscore;/*quantiles*/
set total_zscore;

if kdm^=. then do;
if kdm<30.730842925 then kdmQ=1;
if kdm>=30.730842925 and kdm<40.085588332 then kdmQ=2;
if kdm>=40.085588332 and kdm<51.344047093 then kdmQ=3;
if kdm>=51.344047093 then kdmQ=4;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.617517634 then kdm_advanceQ=1;
if kdm_advance>=-0.617517634 and kdm_advance<-0.10130017 then kdm_advanceQ=2;
if kdm_advance>=-0.10130017 and kdm_advance<0.4244121044 then kdm_advanceQ=3;
if kdm_advance>=0.4244121044 then kdm_advanceQ=4;
end;

if phenoage^=. then do;
if phenoage<28.253669438 then phenoageQ=1;
if phenoage>=28.253669438 and phenoage<39.328905807 then phenoageQ=2;
if phenoage>=39.328905807 and phenoage<49.767785729 then phenoageQ=3;
if phenoage>=49.767785729 then phenoageQ=4;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.670596269 then phenoage_advanceQ=1;
if phenoage_advance>=-0.670596269 and phenoage_advance<-0.218391258 then phenoage_advanceQ=2;
if phenoage_advance>=-0.218391258 and phenoage_advance<0.308272293 then phenoage_advanceQ=3;
if phenoage_advance>=0.308272293 then phenoage_advanceQ=4;
end;

run;

data total_zscore;/*quintiles*/
set total_zscore;

if kdm^=. then do;
if kdm<28.672551733 then kdmP=1;
if kdm>=28.672551733 and kdm<36.245791951 then kdmP=2;
if kdm>=36.245791951 and kdm<44.415169363 then kdmP=3;
if kdm>=44.415169363 and kdm<54.22722268 then kdmP=4;
if kdm>=54.22722268 then kdmP=5;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.735331992 then kdm_advanceP=1;
if kdm_advance>=-0.735331992 and kdm_advance<-0.298550237 then kdm_advanceP=2;
if kdm_advance>=-0.298550237 and kdm_advance<0.1038409515 then kdm_advanceP=3;
if kdm_advance>=0.1038409515 and kdm_advance<0.5608222684 then kdm_advanceP=4;
if kdm_advance>=0.5608222684 then kdm_advanceP=5;
end;

if phenoage^=. then do;
if phenoage<25.846359992 then phenoageP=1;
if phenoage>=25.846359992 and phenoage<35.118369546 then phenoageP=2;
if phenoage>=35.118369546 and phenoage<43.699261628 then phenoageP=3;
if phenoage>=43.699261628 and phenoage<52.024624478 then phenoageP=4;
if phenoage>=52.024624478 then phenoageP=5;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.781024878 then phenoage_advanceP=1;
if phenoage_advance>=-0.781024878 and phenoage_advance<-0.3961624 then phenoage_advanceP=2;
if phenoage_advance>=-0.3961624 and phenoage_advance<-0.032953332 then phenoage_advanceP=3;
if phenoage_advance>=-0.032953332 and phenoage_advance<0.4666604974 then phenoage_advanceP=4;
if phenoage_advance>=0.4666604974 then phenoage_advanceP=5;
end;

run;

proc surveyfreq data=total_zscore;
tables kdmT kdm_advanceT phenoageT phenoage_advanceT
kdmQ kdm_advanceQ phenoageQ phenoage_advanceQ
kdmP kdm_advanceP phenoageP phenoage_advanceP;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
run;

/*crude model*/
proc surveylogistic data=total_zscore;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
class phenoage_advanceP(param=reference ref="1");
/*kdm kdm_advance phenoage phenoage_advance*/
model psonew(ref="0")=phenoage_advanceP/covb;
run;

/*crude model-P-trend*/
proc surveylogistic data=total_zscore;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
/*kdm kdm_advance phenoage phenoage_advance*/
model psonew(ref="0")=kdmP/covb;
run;

/*MICE*/
/*variables nmiss*//*ntotal=10872*/
/*RIAGENDR RIDAGEYR RIDRETH1 DMDEDUC2(7or9)=13 INDFMPIR=720 SDMVPSU SDMVSTRA fiweight yearcycle*/
/*alcoh=812*/
/*LBXCOT=3*/
/*BMXBMI=84*/
/****************************多重填补**********************************************/
data total_missing;
set total;
alqnew=alcoh;
smknew=smk;
bmicnew=bmic;
educnew=educ;
PIRcnew=PIRc;
run;

data total_missing;/*check the frequence of missing data*/
set total_missing;
if alqnew=. then alqnew_flag =1; else alqnew_flag=0;/*1076*/
if smknew=. then smknew_flag=1; else smknew_flag=0;/*3*/
if bmicnew=. then bmicnew_flag=1; else bmicnew_flag=0;/*118*/
if educnew=. then educnew_flag=1;else educnew_flag=0;/*15*/
if PIRcnew=. then PIRcnew_flag=1;else PIRcnew_flag=0;/*784*/
run;

proc freq data=total_missing;
tables alqnew_flag smknew_flag bmicnew_flag educnew_flag PIRcnew_flag;
run;

/*MICE*/
proc mi data=total_missing seed=1000 out=BAPA.totalmi nimpute=10 /*MINIMUM=(.,.,.,.,.)*/;
class alqnew smknew bmicnew educnew PIRcnew; 
var alqnew smknew bmicnew educnew PIRcnew fiweight SDMVPSU SDMVSTRA; 
*fcs nbiter=100 reg(cotnew bminew prinew/details) logistic(marcnew alcohnew smknew/details likelihood=augment) ; 
*fcs nbiter=100 regpmm(cotnew bminew pirnew/details) logistic(marcnew alcohnew smknew/details likelihood=augment); 
*fcs nbiter=100 reg(cotnew bminew pirnew/details) discrim(marcnew alcohnew smknew/classeffects=include details); 
fcs logistic(alqnew smknew bmicnew educnew PIRcnew/details likelihood=augment) nbiter=100; 
run; 

data totalmi;
set BAPA.totalmi;
run;

proc freq data=totalmi;
tables alqnew smknew bmicnew educnew PIRcnew;
run;

/*adjusted for covariables*/
%let factors=RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew PIRcnew PAD700 hyper CHD Hcho DM CKD; 
proc surveylogistic data=totalmi_zscore;
by _imputation_; 
class RIAGENDR(param=reference ref="1");
class agec(param=reference ref="1");
class RIDRETH1(param=reference ref="1");
class bmicnew(param=reference ref="1");
class educnew(param=reference ref="3");
class smknew(param=reference ref="1");
class alqnew(param=reference ref="0");
class PIRcnew(param=reference ref="3");
class PAD700(param=reference ref="2");
class hyper(param=reference ref="0");
class CHD(param=reference ref="0");
class DM(param=reference ref="0");
class Hcho(param=reference ref="0");
class CKD(param=reference ref="0");
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
/*kdm kdm_advance phenoage phenoage_advance*/
model psonew(ref="0")=phenoage_advance &factors/covb;
ods output ParameterEstimates=lgsparmtcsgdm;
run;

proc mianalyze parms(classvar=classval)=lgsparmtcsgdm;
modeleffects Intercept phenoage_advance /*RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew CKD PIRcnew*/;
ods output ParameterEstimates=mi_parms;
run;

data mi_parms;
set mi_parms;
where parm ne 'intercept';
OR=exp(estimate);
LCL_OR=exp(LCLMean);
UCL_OR=exp(UCLMean);
run;

proc print data=mi_parms noobs;
var parm OR LCL_OR UCL_OR;
title 'Combined odds ratio estimates and confidence limits';
run;

proc surveyfreq data=totalmi_zscore;
tables psonew;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
run;

/*frequency*/
proc freq data=totalmi_zscore;
WEIGHT fiweight;
tables kdm kdm_advance phenoage phenoage_advance;
run;

data totalmi_zscore;
set totalmi_zscore;

if kdm^=. then do;
if kdm<33.772901549 then kdmT=1;
if kdm>=33.772901549 and kdm<47.250336148 then kdmT=2;
if kdm>=47.250336148 then kdmT=3;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.437533583 then kdm_advanceT=1;
if kdm_advance>=-0.437533583 and kdm_advance<0.240498593 then kdm_advanceT=2;
if kdm_advance>=0.240498593 then kdm_advanceT=3;
end;

if phenoage^=. then do;
if phenoage<31.918910852 then phenoageT=1;
if phenoage>=31.918910852 and phenoage<46.355834752 then phenoageT=2;
if phenoage>=46.355834752 then phenoageT=3;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.509356109 then phenoage_advanceT=1;
if phenoage_advance>=-0.509356109 and phenoage_advance<0.1039726113 then phenoage_advanceT=2;
if phenoage_advance>=0.1039726113 then phenoage_advanceT=3;
end;

run;

data totalmi_zscore;
set totalmi_zscore;

if kdm^=. then do;
if kdm<30.730842925 then kdmQ=1;
if kdm>=30.730842925 and kdm<40.085588332 then kdmQ=2;
if kdm>=40.085588332 and kdm<51.344047093 then kdmQ=3;
if kdm>=51.344047093 then kdmQ=4;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.617517634 then kdm_advanceQ=1;
if kdm_advance>=-0.617517634 and kdm_advance<-0.10130017 then kdm_advanceQ=2;
if kdm_advance>=-0.10130017 and kdm_advance<0.4244121044 then kdm_advanceQ=3;
if kdm_advance>=0.4244121044 then kdm_advanceQ=4;
end;

if phenoage^=. then do;
if phenoage<28.253669438 then phenoageQ=1;
if phenoage>=28.253669438 and phenoage<39.328905807 then phenoageQ=2;
if phenoage>=39.328905807 and phenoage<49.767785729 then phenoageQ=3;
if phenoage>=49.767785729 then phenoageQ=4;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.670596269 then phenoage_advanceQ=1;
if phenoage_advance>=-0.670596269 and phenoage_advance<-0.218391258 then phenoage_advanceQ=2;
if phenoage_advance>=-0.218391258 and phenoage_advance<0.308272293 then phenoage_advanceQ=3;
if phenoage_advance>=0.308272293 then phenoage_advanceQ=4;
end;

run;

data totalmi_zscore;
set totalmi_zscore;

if kdm^=. then do;
if kdm<28.672551733 then kdmP=1;
if kdm>=28.672551733 and kdm<36.245791951 then kdmP=2;
if kdm>=36.245791951 and kdm<44.415169363 then kdmP=3;
if kdm>=44.415169363 and kdm<54.22722268 then kdmP=4;
if kdm>=54.22722268 then kdmP=5;
end;

if kdm_advance^=. then do;
if kdm_advance<-0.735331992 then kdm_advanceP=1;
if kdm_advance>=-0.735331992 and kdm_advance<-0.298550237 then kdm_advanceP=2;
if kdm_advance>=-0.298550237 and kdm_advance<0.1038409515 then kdm_advanceP=3;
if kdm_advance>=0.1038409515 and kdm_advance<0.5608222684 then kdm_advanceP=4;
if kdm_advance>=0.5608222684 then kdm_advanceP=5;
end;

if phenoage^=. then do;
if phenoage<25.846359992 then phenoageP=1;
if phenoage>=25.846359992 and phenoage<35.118369546 then phenoageP=2;
if phenoage>=35.118369546 and phenoage<43.699261628 then phenoageP=3;
if phenoage>=43.699261628 and phenoage<52.024624478 then phenoageP=4;
if phenoage>=52.024624478 then phenoageP=5;
end;

if phenoage_advance^=. then do;
if phenoage_advance<-0.781024878 then phenoage_advanceP=1;
if phenoage_advance>=-0.781024878 and phenoage_advance<-0.3961624 then phenoage_advanceP=2;
if phenoage_advance>=-0.3961624 and phenoage_advance<-0.032953332 then phenoage_advanceP=3;
if phenoage_advance>=-0.032953332 and phenoage_advance<0.4666604974 then phenoage_advanceP=4;
if phenoage_advance>=0.4666604974 then phenoage_advanceP=5;
end;

run;

proc surveyfreq data=totalmi_zscore;
tables kdmT kdm_advanceT phenoageT phenoage_advanceT
kdmQ kdm_advanceQ phenoageQ phenoage_advanceQ
kdmP kdm_advanceP phenoageP phenoage_advanceP;
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
run;

/*multivariables logistical analysis*/
%let factors=RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew PIRcnew PAD700 hyper CHD Hcho DM CKD; 
proc surveylogistic data=totalmi_zscore;
by _imputation_; 
class RIAGENDR(param=reference ref="1");
class agec(param=reference ref="1");
class RIDRETH1(param=reference ref="1");
class bmicnew(param=reference ref="1");
class educnew(param=reference ref="3");
class smknew(param=reference ref="1");
class alqnew(param=reference ref="0");
class PIRcnew(param=reference ref="3");
class PAD700(param=reference ref="2");
class hyper(param=reference ref="0");
class CHD(param=reference ref="0");
class DM(param=reference ref="0");
class Hcho(param=reference ref="0");
class CKD(param=reference ref="0");
class phenoage_advanceP(param=reference ref="1");
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
/*kdm kdm_advance phenoage phenoage_advance*/
model psonew(ref="0")=phenoage_advanceP &factors/covb;
ods output ParameterEstimates=lgsparmtcsgdm;
run;

proc mianalyze parms(classvar=classval)=lgsparmtcsgdm;
class phenoage_advanceP;
modeleffects Intercept phenoage_advanceP/*RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew CKD PIRcnew*/;
ods output ParameterEstimates=mi_parms;
run;

data mi_parms;
set mi_parms;
where parm ne 'intercept';
OR=exp(estimate);
LCL_OR=exp(LCLMean);
UCL_OR=exp(UCLMean);
run;

proc print data=mi_parms noobs;
var parm phenoage_advanceP OR LCL_OR UCL_OR;
title 'Combined odds ratio estimates and confidence limits';
run;

/*分位后逻辑回归-P-trend*/
/*校正协变量*/
%let factors=RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew PIRcnew PAD700 hyper CHD Hcho DM CKD; 
proc surveylogistic data=totalmi_zscore;
by _imputation_; 
class RIAGENDR(param=reference ref="1");
class agec(param=reference ref="1");
class RIDRETH1(param=reference ref="1");
class bmicnew(param=reference ref="1");
class educnew(param=reference ref="3");
class smknew(param=reference ref="1");
class alqnew(param=reference ref="0");
class PIRcnew(param=reference ref="3");
class PAD700(param=reference ref="2");
class hyper(param=reference ref="0");
class CHD(param=reference ref="0");
class DM(param=reference ref="0");
class Hcho(param=reference ref="0");
class CKD(param=reference ref="0");
strata sdmvstra;
cluster sdmvpsu;
WEIGHT fiweight;
/*kdm kdm_advance phenoage phenoage_advance*/
model psonew(ref="0")=phenoage_advanceP &factors/covb;
ods output ParameterEstimates=lgsparmtcsgdm;
run;

proc mianalyze parms(classvar=classval)=lgsparmtcsgdm;
modeleffects Intercept phenoage_advanceP/*RIAGENDR agec RIDRETH1 bmicnew educnew smknew alqnew CKD PIRcnew*/;
ods output ParameterEstimates=mi_parms;
run;

data mi_parms;
set mi_parms;
where parm ne 'intercept';
OR=exp(estimate);
LCL_OR=exp(LCLMean);
UCL_OR=exp(UCLMean);
run;

proc print data=mi_parms noobs;
var parm OR LCL_OR UCL_OR;
title 'Combined odds ratio estimates and confidence limits';
run;
