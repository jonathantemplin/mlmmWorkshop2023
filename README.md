# Multilevel Measurement Models Workshop 2023

## Instructors: [Lesa Hoffman](https://lesahoffman.com) and [Jonathan Templin](https://jonathantemplin.com)

This repository contains the materials for the Multilevel Measurement Models Workshop 2023, given 27-30 June 2023 at the Universität Mannheim. This workshop was held in the Research Training Group "Statistical Modeling in Psychology", funded by the German Research Foundation (DFG). Details about the event can be found at [https://www.uni-mannheim.de/smip-summerschool/](https://www.uni-mannheim.de/smip-summerschool/). 

Workshop videos are available via YouTube at [https://www.youtube.com/playlist?list=PLSmMs4UgmSMgay00hLnxQNsZHNXN8_ZIt](https://www.youtube.com/playlist?list=PLSmMs4UgmSMgay00hLnxQNsZHNXN8_ZIt).

### Repository Information

The folder structure of the repository is set to correspond with the analysis files found in the root folder. All other paths for necessary files are relative.

Please note, Stan output files (the empty model folder) are not included in the repository. These files will be created when the models are run.

## Workshop Description

This workshop will focus on the use of latent variable measurement models in multilevel sampling designs (e.g., persons within clusters, occasions within persons, stimuli crossed with persons). Course time will be allocated to traditional lectures, guided practice building models, and opportunities for individual data analysis (or further independent practice through instructor-provided data analysis activities). Day 1 will focus on latent variable measurement models for normal, binary, and ordinal responses (all in slope–intercept form) and introduce Stan for MCMC estimation. Day 2 will present concepts of multilevel models using observed outcomes and transition into three-level models for item responses nested in persons nested in clusters. Day 3 will extend multilevel models to include latent variable measurement models with level-specific discrimination parameters. Finally, Day 4 will make connections to models in which item parameters are treated as random effects instead of fixed effects (i.e., for predicting sources of item difficulty and discrimination, as in explanatory item response models). All instructional sessions will be recorded for future participant use.

Prerequisite knowledge and skills include: (1) familiarity with R software for data analysis, (2) some familiarity with Markov Chain Monte Carlo (MCMC) estimation (i.e., have estimated models with MCMC before), (3) some prior knowledge of latent variable measurement models (i.e., confirmatory factor analysis for continuous responses; item response theory for binary and ordinal responses), and (4) some prior knowledge of multilevel models (i.e., hierarchical linear models, mixed-effects models). The course will use Stan software as run through R (using CMDStanR), but no prior experience with Stan is assumed. Participants who wish to use their own devices during the workshop should install Stan ahead of time. No readings will be required ahead of time.


# Guide to Workshop Scripts


| Unit | Syntax | Results | Description |
| --- | --- | --- | --- |
| 01 | 01_Introduction_Bayes_Stan_IRT.docx |  | Quarto document with Word-formatted slides for Introduction to Bayesian Estimation with Stan |
|  | 01_unidimensionalIRT.R | | R script for IRT models (Marginal ML and Bayesian) |
|    | 01_model01.stan | 01_model01.RData | Unidimenionsal 1PL Model |
|    | 01_model02.stan | 01_model02.RData | Unidimenionsal 2PL Model (Standardized Factor) |
|    | 01_model03.stan | 01_model03.RData | Unidimenionsal 2PL Model (Maker Item Discrimination) |
| 02 | 02_Multilevel_Observed.pdf |  | PDF of slides for Multilevel (Observed) Models |
|  | 02_mlm_Bayes.R |  | R script for Multilevel (Non-Measurement) Models with Bayesian Estimation |
|    | 02_mlm_ML.docx |  | R markdown output for Multilevel (Non-Measurement) Models with Maximum Likelihood Estimation |
|    | 02_mlm_ML.R |  | R script for Multilevel (Non-Measurement) Models with Maximum Likelihood Estimation |
|    | 02_model01.stan | 02_model01.RData | Empty (Non-Multilevel) Linear Model for Sum Score |
|    | 02_model02.stan | 02_model02a.RData | Random Intercept Model for Sum Score |
|    |  | 02_model02b.RData | Multilevel Model of Sum Score with Free/Reduced Lunch Smushed |
|    |  | 02_model02c.RData | Multilevel Model of Sum Score with Free/Reduced Lunch at L1 and L2 |
|    |  | 02_model02d.RData | Multilevel Model of Sum Score with Cluster Centered Free/Reduced Lunch at L1 and L2 |
|    | 02_model03.stan | 02_model03.RData | Empty (Non-Multilevel) Generalized Linear Model for Free/Reduced Lunch Status |
|    | 02_model04.stan | 02_model04.RData | Random Intercept Model for Free/Reduced Lunch Status |
|    | 02_model05.stan | 02_model05a.RData | Random Linear Slope Model for Free/Reduced Lunch Status |  
|    |  | 02_model05b.RData | Random Linear Slope Model for Free/Reduced Lunch Status with Cross-Level Interaction |  
| 03 | 03_Multilevel_Measurement.pdf | | PDF of slides for Multilevel Measurement Models |
| | 03_Multilevel_Measurement_Equations.docx | Multilevel measurement model equations and notation |
|  | 03_mlmm.R | | R script for Multilevel Measurement Models with Stan |
|      | 03_model01.stan | 03_model01.RData | Empty (Non-Measurement) Two-Level Model with Correlated Random Item Intercepts |
|      | 03_model02.stan | 03_model02.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Within-School Discriminations Fixed=1 | 
|      | 03_model03.stan | 03_model03.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Estimated Within-School Discriminations using Standardized Theta | 
|      | 03_model04.stan | 03_model04.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Estimated Within-School Discriminations using Item1=Marker |
|      | 03_model05.stan | 03_model05.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Estimated Level-Specific WS (Item1=Marker) and BS (Item10=Marker) Discriminations |
|      | 03_model06.stan | 03_model06.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
|      | 03_model07.stan | 03_model07.RData | Within-School and Between-School Measurement Model without Random Item Intercepts and with Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
|      | 03_model08.stan | 03_model08.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Free/Reduced Lunch MLM Predictor and Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
| Other | Other_Materials_from_Lesa_Hoffman.docx | | Word document with additional materials from Lesa Hoffman |
| | | | | 
