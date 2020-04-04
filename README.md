# [Hierarchical spatial capture-recapture models: modelling population density in stratified populations.](https://doi.org/10.1111/2041-210X.12135)

### Royle J.A. and S.J. Converse

### Methods in Ecology and Evolution 

## Abstract:
1. Capture–recapture studies are often conducted on populations that are stratified by space, time or other factors. In this paper, we develop a Bayesian spatial capture–recapture (SCR) modelling framework for stratified populations – when sampling occurs within multiple distinct spatial and temporal strata.
2. We describe a hierarchical model that integrates distinct models for both the spatial encounter history data from capture–recapture sampling, and also for modelling variation in density among strata. We use an implementation of data augmentation to parameterize the model in terms of a latent categorical stratum or group membership variable, which provides a convenient implementation in popular BUGS software packages.
3. We provide an example application to an experimental study involving small‐mammal sampling on multiple trapping grids over multiple years, where the main interest is in modelling a treatment effect on population density among the trapping grids.
4. Many capture–recapture studies involve some aspect of spatial or temporal replication that requires some attention to modelling variation among groups or strata. We propose a hierarchical model that allows explicit modelling of group or strata effects. Because the model is formulated for individual encounter histories and is easily implemented in the BUGS language and other free software, it also provides a general framework for modelling individual effects, such as are present in SCR models.

## Code 
1. [IPM_analysis](./IPM_analysis/): This folder contains the code to load data and run the Integrated Population Model and goodness-of-fit for Chukchi Sea polar bears. It also contains the code to generate the JAGS file for the Bayesian implementation.

2. [recruitment_analysis](./recruitment_analysis/): This folder contains the code to load data and run the recruitment analysis (yearlings [c1s] per adult female). It also contains the code to generate the JAGS file for the Bayesian implementation.

3. [density_extrapolation_analysis](./density_extrapolation_analysis/): This folder contains the code to run the density extrapolation.


## Data

Code (JAGS and NIMBLE) and example data for implementing hierarchical spatial capture-recapture models

The original analysis was completed with the JAGS code. As of early 2020, this code is not compiling, apparently due to changes in JAGS. I haven't had to time work out how to fix it.  

The NIMBLE version of the code (written by Nathan Hostetter) does compile and run efficiently. 
