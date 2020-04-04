# [Hierarchical spatial capture-recapture models: modelling population density in stratified populations.](https://doi.org/10.1111/2041-210X.12135)

### Royle J.A. and S.J. Converse

### Methods in Ecology and Evolution 

## Abstract:
1. Capture–recapture studies are often conducted on populations that are stratified by space, time or other factors. In this paper, we develop a Bayesian spatial capture–recapture (SCR) modelling framework for stratified populations – when sampling occurs within multiple distinct spatial and temporal strata.
2. We describe a hierarchical model that integrates distinct models for both the spatial encounter history data from capture–recapture sampling, and also for modelling variation in density among strata. We use an implementation of data augmentation to parameterize the model in terms of a latent categorical stratum or group membership variable, which provides a convenient implementation in popular BUGS software packages.
3. We provide an example application to an experimental study involving small‐mammal sampling on multiple trapping grids over multiple years, where the main interest is in modelling a treatment effect on population density among the trapping grids.
4. Many capture–recapture studies involve some aspect of spatial or temporal replication that requires some attention to modelling variation among groups or strata. We propose a hierarchical model that allows explicit modelling of group or strata effects. Because the model is formulated for individual encounter histories and is easily implemented in the BUGS language and other free software, it also provides a general framework for modelling individual effects, such as are present in SCR models.

## Code 
1. Analysis_MEE.r: This folder contains the original JAGS code used to run (an example data set) for the paper. This currently does not compile, due apparently to changes in JAGS. I have not updated it.  

2. Analysis_MEE_updated.r: This folder contains a NIMBLE implementation of the model, courtesy of Nathan Hostetter. 

## Data

Example_MEE.RData: R workspace containing the data objects necessary to run the model. 
