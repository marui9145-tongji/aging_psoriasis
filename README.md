# aging_psoriasis
Statistical analysis based on NHANES database to explore the association between psoriasis and aging by SAS and R
Univariable, multivariable and mediation Mendelian randomization to investigate the independent effects of psoriasis on aging by R

These codes are for the paper "Association between aging and psoriasis: evidence from observational study and bidirectional two-sample Mendelian randomization analysis" submitted to Nature Communications. These codes are contributed by Dr. Rui Ma.

Cross_sectional_study_NHANES_byR.R is a Rscript used for performing Kruskal-Wallis H test, weighted restricted cubic spline analysis, calculating FDR adjusted P-value based on NHANES database, and drawing forest plot and heatmap by R software (version 4.3.3).

Cross_sectional_study_NHANES_bySAS.sas is a shell script used for statistical analysis based on NHANES database by SAS software (version 9.4 TS1M7).

MR_byR.R is a Rscript used for conducting univariable, multivariable, and mediation Mendelian Randomization analysis by R software (version 4.3.3).

In addition: Computation of KDM-BA, PhenoAge was conducted using the R package “BioAge” (https://github.com/dayoonkwon/BioAge); Two-sample MR analyses were conducted using the TwoSampleMR package (version 0.5.5)(https://github.com/MRCIEU/TwoSampleMR) in R . MR-PRESSO was performed using MRPRESSO (version 1.0) ((https://github.com/rondolab/MR-PRESSO)) in R; Multivariable MR analysis was conducted using both the MVMR (version 0.2.0) (https://github.com/WSpiller/RMVMR) and MendelianRandomization (version 0.5.0) ((https://cran.r-project.org/web/packages/MendelianRandomization) packages in R. 
