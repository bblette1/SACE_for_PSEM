# Set working directory to "No Early Treatment Effects Scenario" folder
rm(list=ls())

nsims <- 2000
set.seed(1234)

library(rootSolve)
library(numDeriv)
library(plyr)

source("simulate_EER.R")
source("analyze_EER.R")
source("low_eefun_EER.R")
source("up_eefun_EER.R")
source("../M Estimation Helper Functions/compute.R")
source("../M Estimation Helper Functions/utilities.R")

# Define helper function for Imbens-Manski interval computation
f_cstar <- function(c, low, up, maxsig) {
  pnorm(c + (up - low) / maxsig) - pnorm(-c) - 0.95
}

# Values corresponding to CEP(1,0) - CEP(0,0) estimand (plots' x-axis)
# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 0, S^tau(0) = 0) = a
varyingprobs1 <- seq(0.7, 0.1, by = -0.1)
# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 1, S^tau(0) = 0) = b
varyingprobs2 <- seq(0.1, 0.7, by = 0.1)

# Sim across proportion of non-cases with S^tau observed
# 1 corresponds to full cohort
for (casecohort in c(1, 0.25, 0.1)) {

# Results array: save six results per sim
# 1 CEP_l, 2 sigma^2_l, 3 CEP_u, 4 sigma^2_u
# 5 EUI lower limit, 6 EUI upper limit
results <- array(NA, dim = c(3, 4, length(varyingprobs1), nsims, 6))

# Sim across sensitivity parameter ranges
for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-1, 1)
  if (h == 3) beta0range <- c(-2.5, 2.5)
  
  # Sim across sample sizes
  for (i in 1:4) {

    if (i == 1) n <- 200
    # Skip low sample size and small proportion sampled
    if (i == 1 & casecohort == 0.1) next
    if (i == 2) n <- 400
    if (i == 3) n <- 800
    if (i == 4) n <- 1600
   
    # Sim across different values of estimand: CEP(1,0) - CEP(0,0)
    for (j in 1:length(varyingprobs1)) {
      
      # P(Y^tau(1) = a, Y^tau(0) = b) for (a,b) in {(0,0),(0,1),(1,0),(1,1)}
      Y_tau_probvector <- c(.8, 0, 0, .2)
      # P(S^tau(1) = a, S^tau(0) = b | Y^tau(0) = 1)
      S_star_probvector_ytau0eq1 <- c(.2, 0, .3, .5)
      # P(S^tau(1) = a, S^tau(0) = b | Y^tau(0) = 0)
      S_star_probvector_ytau0eq0 <- c(.4, 0, .6, 0) 
      # P(Y(1) = 1 | S^tau(1) = a, S^tau(0) = b, Y^tau(0) = 0, Y^tau(0) = 0)
      Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      # P(Y(0) = 1 | S^tau(1) = a, S^tau(0) = b, Y^tau(0) = 0, Y^tau(0) = 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
            
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      TrueDiff <- TrueCEP_10 - TrueCEP_00      
      
      for (k in 1:nsims) {

        # Simulate full observed trial data given parameters for each sim
        observed_data <- simulate_data(Y_tau_probvector,
                                       S_star_probvector_ytau0eq1,
                                       S_star_probvector_ytau0eq0,
                                       Y_1_probvector, Y_0_probvector, n,
                                       ccprob = casecohort)
        
        # Get point estimates if full cohort and prepare data for
        # variance estimation
	      if (casecohort == 1) {
	        ad <- analyze_data(data = observed_data, brange=beta0range)
		      thetahat_low_B <- as.numeric(c(ad[1:4], ad[5], ad[7], ad[9]))
	  	    thetahat_up_B <-  as.numeric(c(ad[1:4], ad[6], ad[8], ad[10]))
      		low_fun <- low_eefun
      		up_fun <- up_eefun
      		observed_data_wide <- count(observed_data,
      		                            c('Y', 'Z', 'Y_tau', 'S_star'))
      		observed_data_wide$Group <- 1:nrow(observed_data_wide)
	      }
        
        # Get point estimates if not full cohort and prepare data for
        # variance estimation
    	  if (casecohort < 1) {
    	    # R: indicator S^tau == S_star is measured
    	  	observed_data$R <- 1*!is.na(observed_data$S_star)
    	  	# set S_star to zero instead of NA
    	  	observed_data$S_star[observed_data$R == 0] <- 0
    		  ad <- analyze_data_cc(data=observed_data, brange=beta0range)
    	  	thetahat_low_B <- as.numeric(c(ad[1:4], ad[5], ad[7], ad[9],
    	  	                               ad[11]))
    	  	thetahat_up_B <-  as.numeric(c(ad[1:4], ad[6], ad[8], ad[10],
    	  	                               ad[11]))
      		low_fun <- low_eefun_cc
      		up_fun <- up_eefun_cc
      		observed_data_wide <- count(observed_data,
      		                            c('Y', 'Z', 'Y_tau', 'S_star', 'R'))
      		observed_data_wide$Group <- 1:nrow(observed_data_wide)
    		}

        # Calculate variance estimate for the lower bound of the ignorance
        # intervals using 'geex' package source code
    	  mats <-
    	    compute_matrices(list(eeFUN = low_fun,
    	                          splitdt = split(observed_data_wide,
    	                                       f = observed_data_wide$Group)),
    	                     numDeriv_options = list(method = 'simple'),
    	                     theta = thetahat_low_B,
    	                     beta0range = beta0range)
    
    	  A <- apply(simplify2array(Map(`*`, mats$A_i,
    	                                observed_data_wide$freq)),
    	             1:2, sum)
    	  B <- apply(simplify2array(Map(`*`, mats$B_i,
    	                                observed_data_wide$freq)),
    	             1:2, sum)
    
    	  Sigma <- solve(A) %*% B %*% t(solve(A))
        
    	  # Extract variance for estimand of interest from sandwich matrix
      	var_Ignorance_low_B_M <- Sigma[7,7]
    
      	# Same process for upper bound
    	  if (min(beta0range) == max(beta0range)) {
    	  	var_Ignorance_up_B_M <- var_Ignorance_low_B_M
    	  } else {
    		  mats <- 
    		    compute_matrices(list(eeFUN = up_fun,
    		                          splitdt = split(observed_data_wide,
    		                                    f = observed_data_wide$Group)),
    		                     numDeriv_options = list(method = 'simple'),
    		                     theta = thetahat_up_B, beta0range = beta0range)
    		  
    	  	A <- apply(simplify2array(Map(`*`, mats$A_i,
    	  	                              observed_data_wide$freq)),
    	  	           1:2, sum)
    	  	B <- apply(simplify2array(Map(`*`, mats$B_i,
    	  	                              observed_data_wide$freq)),
    	  	           1:2, sum)
    	  	Sigma <- solve(A) %*% B %*% t(solve(A))
      	  var_Ignorance_up_B_M <- Sigma[7,7]
    	  } # end of 'else'
    		
    	  results[h, i, j, k, 1] <- as.numeric(ad$II_low)
    	  results[h, i, j, k, 2] <- var_Ignorance_low_B_M
    	  results[h, i, j, k, 3] <- as.numeric(ad$II_up)
    	  results[h, i, j, k, 4] <- var_Ignorance_up_B_M
    
      	# compute Imbens-Manski interval
      	maxsig <- max(sqrt(results[h, i, j, k, 2]),
      	              sqrt(results[h, i, j, k, 4]))
        cB <- uniroot(f_cstar, c(1.64, 1.96), low = results[h, i, j, k, 1],
                      up = results[h, i, j, k, 3], maxsig = maxsig)$root       
      	results[h, i, j, k, 5] <- results[h, i, j, k, 1] - 
      	                          cB*sqrt(results[h, i, j, k, 2])
      	results[h, i, j, k, 6] <- results[h, i, j, k, 3] +
      	                          cB*sqrt(results[h, i, j, k, 4])
      
    	} # end of 'k' loop 

    } # end of 'j' loop
  
    print(paste(casecohort, h, i))

  } # end of 'i' loop
  
} # end of 'h' loop

outfile <- paste("results", nsims, "sims.", casecohort, "cc.Rdata", sep="")
save(results, file = outfile)

} # end of 'casecohort' loop