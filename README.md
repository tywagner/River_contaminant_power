# River_contaminant_power
Wagner, T., McLaughlin, P., Smalling, K., Breitmeyer, S., Gordon, S., and Noe, G.B. The statistical power to detect regional temporal trends in riverine contaminants in the Chesapeake Bay Watershed, USA. Science of the Total Environment.

R scripts and JAGS model code for fitting the censored hierarchical components of variance model described in Wagner et al. and for performing power simulations. 

Files:
vc.txt: model written in “BUGS” language for fitting using JAGS

roar_power_sims_contaminants3.R: R code for performing power simulations

roar_power_sims_contaminants3_skipYear2.R: R code for performing power simulations across different sampling frequencies

All code and simulations ran for Wagner et al. were done on the Roar Supercomputer at Penn State University Institute for Computational and Data Sciences. Some changes to the R code may be necessary to run on different platforms. 
