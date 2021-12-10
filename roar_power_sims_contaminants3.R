#get job ID from cluster
jobid = as.integer(Sys.getenv("INPUT")) 
set.seed(jobid)

#load libraries 
library(dplyr)
library(jagsUI)
library(lme4)

# Read in mcmc output
#out <- readRDS(file='estro_vc.rds')
out <- readRDS(file='estro_vc.rds')
# Number of simulation (and thus MCMC samples used)
num.sim <- 1
# Chain length from analysis
chainLength <- out$mcmc.info$n.samples
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
#ID = seq( 1 , chainLength , floor(chainLength/out$mcmc.info$n.thin) )
ID <- seq(1,chainLength,3)

#set levels for n.sites, max.year and beta1
siteLevels <- c(10,20,30,50,80,100)
maxYearLevels <- c(5,10,20)
beta1Levels <- c(-0.05,-.07,-0.1,-0.15,-0.2)
#store length of each level
nSiteLevels <- length(siteLevels)
nMaxYearLevels <- length(maxYearLevels)
nBeta1Levels <- length(beta1Levels)
#number of samples for each cell
nRep <- 100


# Settings
# Number of sample sites
n.sites <- rep(rep(siteLevels,nBeta1Levels*nMaxYearLevels),nRep)[jobid]

# Years
max.year <- rep(rep(rep(maxYearLevels,each=nSiteLevels),nBeta1Levels),nRep)[jobid]
year <- 1:max.year
# Proportion of data to censor
cenProp <- 0.10

# Number of years
n.years <- length(year)
# Number of observations per site/year
obs <- rep(1:3, n.sites*n.years)
# length(obs)

# Set effect (trend) magnitudes to evaluate
# The code currently looks over beta1 to allow for evaluating multiple slopes during one simulation
# but we might not want to do this for the server??
beta1 <- rep(rep(beta1Levels,each=nMaxYearLevels*nSiteLevels),nRep)[jobid]

# create empty results object to bind jags results into. This is a list, with each element 
# of the list corresponding to a different effect size, i.e., corresponding to a different value in beta1, above
results <- list()
for(m in 1:length(beta1)){ # length(beta1) * n.years
  results[[m]] <- list(NULL)
}

# Create container to to hold posterior estimates of beta1 from each simulation run
bPosterior <- list()

start.time = Sys.time()  #start timer
for(j in 1:length(beta1)){ # loop over beta1's
  
  
  # parameters to save during simulations
  ### matrix to save result
  pars.to.save <- c("beta0", "beta1.mean","sigma.e","sigma.slope",
                    "sigma.site", "sigma.year", "sigma.site.year")
  resultCols<-c('parameter','Mean','SD','2.5%','97.5%','rHat','iterNumber','effect_size')
  res = array(NA,dim=c(length(pars.to.save),length(resultCols))) %>%
    data.frame()
  names(res)<-resultCols
  res$parameter<-pars.to.save
  
  # Loop over the number of simulations
  for(sim in 1:num.sim){
    
    
    # Read in parameters, ID[sim] is grabbing a different set of parameters from the joint posterior
    # for every simulation to account for parameter uncertainty
    # Fixed intercept
    beta0 <- mean(out$sims.list$beta0[ID])
    # beta1 <- out$sims.list$beta1.mean[ID]
    # Random effects
    sigma.e <- mean(out$sims.list$sigma.e[ID]) # residual sd
    sigma.slope <- mean(out$sims.list$sigma.slope[ID]) # slope sd
    sigma.site <- mean(out$sims.list$sigma.site[ID]) # among site sd
    sigma.year <- mean(out$sims.list$sigma.year[ID]) # among year sd
    sigma.site.year <- mean(out$sims.list$sigma.site.year[ID]) # site-by-year interaction sd
    
    # Generate random effects for data simulation
    # site effect
    site.true <- rnorm(n.sites, mean = beta0, sd = sigma.site)
    # slope (time effect)
    trend.true <- rnorm(n.sites, mean = beta1[j], sd = sigma.slope)
    # year effect
    year.true <- rnorm(n.years, mean = 0, sd = sigma.year)
    # site-by-year effect
    site.year.true <- matrix(rnorm(n.sites*n.years, mean = 0, sd = sigma.site.year), nrow=n.sites, ncol=n.years)
    
    
    # Create some identifiers for fitting the model
    # Site identifier
    sites <- rep(1:n.sites, each=n.years)
    # Year identifier
    years <- rep(1:n.years, length=n.sites * n.years)
    # Time series for each site
    time <- rep(year, n.sites)
    # Standardize time for use as predictor variable
    z_year <- year-1
    
    # Simulate dummy 'experiment' for lmer model below to get design matrices
    y <- rnorm(length(obs), site.true[sites] + trend.true[sites] * z_year + 
                 year.true[years] + site.year.true[sites,years], sigma.e)
    
    # Add y to data frame
    dat <- data.frame(y, sites, time, z_year)
    
    # Fit lmer model to obtain design matrices for simulating actual data
    skip_to_next <- FALSE
    tryCatch(m1 <- lmer(y ~ 1 + z_year +  (z_year||sites) + (1|time) + (1|sites:time), data=dat),
             error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    # Get design matrix Z and X for simulating data below
    design_matrix_Z <- getME(m1, "Z") %>% as.matrix()
    design_matrix_X <- getME(m1, "X") %>% as.matrix()
    # design_matrix_Z %>% head(1)
    
    # Bundle true random effects
    random_effects_intercepts <- c(site.true, trend.true, year.true, site.year.true)
    
    # Bundle fixed effects
    fixed.effects <- c(beta0, beta1[j])
    
    
    # create observations 
    # multiply the design matrix Z by the vector of random effects and the fixed 
    # effects by the design matrix X and 
    # add in a residual error. 
    # Simualte log-normally distributed values 
    dat2 <- dat %>%
      mutate(measurement = exp(rnorm(n = length(obs), mean = design_matrix_X %*% fixed.effects + design_matrix_Z %*% random_effects_intercepts, 
                                     sd = sigma.e)
      ))
    
    # Log-transform simulated response measurements
    dat2 <- dat2 %>% 
      mutate(measurement = log(measurement))
    
    
    # Censor some proportion of the data
    dat2$censorLimitVec <- as.numeric(rep(quantile(dat2$measurement,cenProp ),length = dim(dat2)[1]))
    dat2$isNotCensored <- ( dat2$measurement > dat2$censorLimitVec ) # Must tell JAGS which observations are ABOVE Censoring limit
    dat2$measurement[!dat2$isNotCensored] <- NA
    
    
    
    # Load data
    data <- list(y = as.numeric(dat2$measurement), n = dim(dat2)[1], site = dat2$sites, year = dat2$time, yearX = dat2$z_year,
                 J = length(unique(dat2$sites)), Y = length(unique(dat2$time)),
                 censorLimitVec = dat2$censorLimitVec,
                 isNotCensored = as.numeric(dat2$isNotCensored) )
    
    
    # Initial values
    yInit <- rep( NA , nrow(dat2) )
    yInit[!dat2$isNotCensored] <- dat2$censorLimitVec[!dat2$isNotCensored]*0.5
    yInit <- yInit
    
    inits <- function (){
      list (beta0=rnorm(1), beta1.mean=rnorm(1), sigma.e=runif(1), sigma.slope = runif(1), 
            sigma.site = runif(1), sigma.year = runif(1), sigma.site.year = runif(1), y=yInit )
    }
    
    # Parameters monitored
    parameters <- c("beta0", "beta1.mean","sigma.e","sigma.slope",
                    "sigma.site", "sigma.year", "sigma.site.year")
    
    
    # MCMC settings
    ni <- 50000
    nt <- 1
    nb <- 25000
    nc <- 3
    
    # Jags could crash due to random number generation, so using try to avoid nullifying the whole simulation
    jagsWorked<-FALSE
    try(expr={
      out <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, 
                    model.file = "vc.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
                    n.burnin = nb, parallel = T)
      jagsWorked<-TRUE
    })
    # print(out)
    
    if(!jagsWorked){
      #if jags throws an error, then just save NAs as the result to indicate which models failed
      res[,c('Mean','SD','2.5%','97.5%','rHat')] <- as.numeric(NA)
      res$iterNumber<- sim
      res$effect_size <- beta1[j]
      results[[j]]<-rbind(results[[j]],res)
      next
    }
    
    #save the results summary and simulation info
    res[,c('Mean','SD','2.5%','97.5%','rHat')] <- out$summary[c(c("beta0", "beta1.mean","sigma.e","sigma.slope",
                                                                  "sigma.site", "sigma.year", "sigma.site.year")),c(1,2,3,7,8)]
    res$iterNumber <- sim
    res$effect_size <- beta1[j]
    
    #bind with previous simulations
    results[[j]]<-rbind(results[[j]],res)
    
    #Retain posterior samples for b - effect of temp
    name <- paste("b:",j,sim,sep='')
    bPosterior[[name]] <- out$sims.list$beta1.mean
    
    #save 'true' concentrations
    conc_start <- exp(beta0)
    conc_end <- exp(beta0 + max.year*beta1)
    results[[1]]$conc_start <- conc_start
    results[[1]]$conc_end <- conc_end
    
    #save n.sites, max.year and beta1
    results[[1]]$beta1 <- beta1
    results[[1]]$n.sites <- n.sites
    results[[1]]$max.year <-  max.year
    results[[1]]$jobid <-  jobid
    
  } # end num.sim loop
}  # end scenario loop
end.time = Sys.time() #end timer
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
results


#save results
write.csv(results, file=paste0('results',jobid,'.csv'))
