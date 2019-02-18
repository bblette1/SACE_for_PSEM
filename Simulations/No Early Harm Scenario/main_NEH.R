# Set working directory to "No Early Harm Scenario" folder
rm(list=ls())

nsims <- 20
set.seed(1234)

library(rootSolve)
library(numDeriv)
library(plyr)

source("simulate_NEH.R")
source("analyze_NEH.R")
source("low_eefun_NEH.R")
source("up_eefun_NEH.R")
source("../../M Estimation Helper Functions/compute.R")
source("../../M Estimation Helper Functions/utilities.R")

# Define helper function for Imbens-Manski interval computation
f_cstar <- function(c, low, up, maxsig){
  pnorm(c + (up - low) / maxsig) - pnorm(-c) - 0.95
}

# Values corresponding to CEP(1,0) - CEP(0,0) estimand (plots' x-axis)
# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 0, S^tau(0) = 0) = a
varyingprobs1 <- seq(0.7, 0.1, by = -0.05)
# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 1, S^tau(0) = 0) = b
varyingprobs2 <- seq(0.1, 0.7, by = 0.05)

# Sim across proportion of non-cases with S^tau observed
# 1 corresponds to full cohort
for (casecohort in c(0.1, 0.25, 1)){	

# Results array: save six results per sim
# 1 CEP_l, 2 sigma^2_l, 3 CEP_u, 4 sigma^2_u
# 5 EUI lower limit, 6 EUI upper limit
results <- array(NA,dim=c(3,3,length(varyingprobs1),nsims,6))

# Sim across sensitivity parameter ranges
for (h in 1:3) {

  if (h == 1) beta0range <- beta1range <- beta2range <- beta3range <- c(0,0)
  if (h == 2) beta0range <- beta1range <- beta2range <- beta3range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- beta1range <- beta2range <- beta3range <- c(-1, 1)
  
  # Sim across sample sizes
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Sim across different values of estimand: CEP(1,0) - CEP(0,0)
    for (j in 1:length(varyingprobs1)) {
      
      # P(Y^tau(1) = a, Y^tau(0) = b) for (a,b) in {(0,0),(0,1),(1,0),(1,1)}
      Y_tau_probvector <- c(.7, .2, 0, .1)
      # P(S^tau(1) = a, S^tau(0) = b | Y^tau(1) = 0, Y^tau(0) = 0)
      S_star_probvector_ytaueq00 <- c(.4, 0, .6, 0)
      # P(S^tau(1) = a, S^tau(0) = NA | Y^tau(1) = 0, Y^tau(0) = 1)
      S_star_probvector_ytaueq01 <- c(.4, .6)
      # P(Y(1) = 1 | S^tau(1) = a, S^tau(0) = b, Y^tau(0) = 0, Y^tau(0) = 0)
      Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      # P(Y(0) = 1 | S^tau(1) = a, S^tau(0) = b, Y^tau(0) = 0, Y^tau(0) = 0)
      Y_0_probvector <- c(.5, 0, .5, 0)
            
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      TrueDiff <- TrueCEP_10 - TrueCEP_00      
      
      for (k in 1:nsims) {

        # Simulate full observed trial data given parameters for each sim
        observed_data <- simulate_data(Y_tau_probvector,
                                       S_star_probvector_ytaueq00,
                                       S_star_probvector_ytaueq01,
                                       Y_1_probvector, Y_0_probvector,
                                       n, ccprob = casecohort)	  
        
        # Get point estimates if full cohort and prepare data for
        # variance estimation
    	  if (casecohort == 1) {
      	 	ad <- analyze_data(data = observed_data, brange0 = beta0range,
      	 	                   brange1 = beta1range, brange2 = beta2range,
      	 	                   brange3 = beta3range)
      	 	# Specify positions of output vector from analyze_data()
      	 	# that correspond to items in estimating equation vector
      		thetahat_low_C <- as.numeric(c(ad[1:5], ad[8:9], ad[12:15],
      		                               ad[18:24], ad[27]))
      		thetahat_up_C <- as.numeric(c(ad[1:5], ad[10:17], ad[20:22],
      		                              ad[25:26], ad[28]))
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
      	 	ad <- analyze_data_cc(data = observed_data,
      	 	                        brange0 = beta0range,
      	 	                        brange1 = beta1range,
      	 	                        brange2 = beta2range,
      	 	                        brange3 = beta3range)
      		thetahat_low_C <- as.numeric(c(ad[1:5], ad[8:9], ad[12:15],
      		                               ad[18:24], ad[27], ad[29]))
      		thetahat_up_C <- as.numeric(c(ad[1:5], ad[10:17], ad[20:22],
      		                              ad[25:26], ad[28:29]))		
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
    	                     theta = thetahat_low_C, beta0range = beta0range,
    	                     beta1range = beta1range, beta2range = beta2range,
    	                     beta3range = beta3range)
    
    	  A <- apply(simplify2array(Map(`*`, mats$A_i,
    	                                observed_data_wide$freq)), 1:2, sum)
    	  B <- apply(simplify2array(Map(`*`, mats$B_i,
    	                                observed_data_wide$freq)), 1:2, sum)
    
    	  Sigma <- solve(A) %*% B %*% t(solve(A))
    
    	  # Extract variance for estimand of interest from sandwich matrix
      	var_Ignorance_low_C_M <- Sigma[19, 19]

      	# Same process for upper bound
    	  if (min(beta0range) == max(beta0range)) {
    	  	var_Ignorance_up_C_M <-  var_Ignorance_low_C_M
    	  } else {
    		  mats <- 
    		    compute_matrices(list(eeFUN = up_fun,
    		                          splitdt = split(observed_data_wide,
    		                                    f = observed_data_wide$Group)),
    		                     numDeriv_options = list(method = 'simple'),
    		                     theta = thetahat_up_C, beta0range = beta0range,
    		                     beta1range = beta1range,
    		                     beta2range = beta2range,
    		                     beta3range = beta3range)
    		  
    	  	A <- apply(simplify2array(Map(`*`, mats$A_i,
    	  	                              observed_data_wide$freq)), 1:2, sum)
    	  	B <- apply(simplify2array(Map(`*`, mats$B_i,
    	  	                              observed_data_wide$freq)), 1:2, sum)
    	  	
    	  	Sigma <- solve(A) %*% B %*% t(solve(A))
    	  	var_Ignorance_up_C_M <- Sigma[19, 19]
    	  } # end of 'else'

		
	      results[h, i, j, k, 1] <- as.numeric(ad$II_low)
	      results[h, i, j, k, 2] <- var_Ignorance_low_C_M
	      results[h, i, j, k, 3] <- as.numeric(ad$II_up)
	      results[h, i, j, k, 4] <- var_Ignorance_up_C_M

      	# compute Imbens-Manski interval
      	maxsig <- max(sqrt(results[h, i, j, k, 2]),
      	              sqrt(results[h, i, j, k, 4]))
        cC <- uniroot(f_cstar, c(1.64, 1.96), low = results[h, i, j, k, 1],
                      up = results[h, i, j, k, 3], maxsig = maxsig)$root       
      	results[h, i, j, k, 5] <- results[h, i, j, k, 1] - 
      	                          cC*sqrt(results[h, i, j, k, 2])
      	results[h, i, j, k, 6] <- results[h, i, j, k, 3] + 
      	                          cC*sqrt(results[h, i, j, k, 4])
      
      } # end of 'k' loop 

    } # end of 'j' loop

  } # end of 'i' loop

} # end of 'h' loop

outfile <- paste("results", nsims, "sims.", casecohort, "cc.Rdata", sep="")
save(results, file=outfile)

} # end of 'casecohort' loop