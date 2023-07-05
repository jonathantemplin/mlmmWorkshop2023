
# Guide to Workshop Scripts


| Unit | Syntax | Results | Description |
| --- | --- | --- | --- |
| 01 | 01_unidimensionalIRT.R | | R script for IRT models (Marginal ML and Bayesian) |
|    | 01_model01.stan | 01_model01.RData | Unidimenionsal 1PL Model |
|    | 01_model02.stan | 01_model02.RData | Unidimenionsal 2PL Model (Standardized Factor) |
|    | 01_model03.stan | 01_model03.RData | Unidimenionsal 2PL Model (Maker Item Discrimination) |
| 02 | 02_mlm_Bayes.R |  | R script for Multilevel (Non-Measurement) Models with Bayesian Estimation |
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
| 03 | 03_mlmm.R | | R script for Multilevel Measurement Models with Stan |
|      | 03_model01.stan | 03_model01.RData | Empty (Non-Measurement) Two-Level Model with Correlated Random Item Intercepts |
|      | 03_model02.stan | 03_model02.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Within-School Discriminations Fixed=1 | 
|      | 03_model03.stan | 03_model03.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Estimated Within-School Discriminations using Standardized Theta | 
|      | 03_model04.stan | 03_model04.RData | Within-School Measurement Model with Correlated Random Item Intercepts and Estimated Within-School Discriminations using Item1=Marker |
|      | 03_model05.stan | 03_model05.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Estimated Level-Specific WS (Item1=Marker) and BS (Item10=Marker) Discriminations |
|      | 03_model06.stan | 03_model06.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
|      | 03_model07.stan | 03_model07.RData | Within-School and Between-School Measurement Model without Random Item Intercepts and with Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
|      | 03_model08.stan | 03_model08.RData | Within-School and Between-School Measurement Model with Uncorrelated Random Item Intercepts and Free/Reduced Lunch MLM Predictor and Estimated Level-Constrained WS (Item1=Marker) and BS (Item1=Marker) Discriminations |
| | | | | 
