# Simulate trial data under NETE scenario
simulate_data <- function(Y_tau_probvector, S_star_probvector_ytau0eq1,
                          S_star_probvector_ytau0eq0, Y_1_probvector,
                          Y_0_probvector, n, ccprob) {

	multi <- rmultinom(1, n, Y_tau_probvector)
	Y_tau_1 <- c(rep(0, multi[1, 1] + multi[2, 1]),
	             rep(1, multi[3, 1] + multi[4, 1]))
	Y_tau_0 <- c(rep(0, multi[1, 1] + multi[3, 1]),
	             rep(1, multi[2, 1] + multi[4, 1]))
      
	multi2 <- rmultinom(1, length(Y_tau_0[Y_tau_0 == 1]),
	                    S_star_probvector_ytau0eq1)
	S_star_1 <- rep(0, n)
	S_star_1[Y_tau_0 == 1] <- c(rep(0, multi2[1, 1] + multi2[2, 1]),
	                            rep(1, multi2[3, 1] + multi2[4, 1]))
	S_star_0 <- rep(0,n)
	S_star_0[Y_tau_0 == 1] <- c(rep(0, multi2[1, 1] + multi2[3, 1]),
	                            rep(1, multi2[2, 1] + multi2[4, 1]))
       
	multi3 <- rmultinom(1, length(Y_tau_0[Y_tau_0 == 0]),
	                    S_star_probvector_ytau0eq0)
	S_star_1[Y_tau_0 == 0] <- c(rep(0, multi3[1, 1] + multi3[2, 1]),
	                            rep(1, multi3[3, 1] + multi3[4, 1]))
	S_star_0[Y_tau_0 == 0] <- c(rep(0, multi3[1, 1] + multi3[3, 1]),
	                            rep(1, multi3[2, 1] + multi3[4, 1]))
        
	Y_1 <- Y_tau_1
	Y_1[Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 0 & S_star_0 == 0] <- 
	  rbinom(length(Y_1[Y_1 == 0 & S_star_1 == 0 & S_star_0 == 0]), 1,
	         Y_1_probvector[1])
	Y_1[Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 1 & S_star_0 == 0] <- 
	  rbinom(length(Y_1[Y_1 == 0 & S_star_1 == 1 & S_star_0 == 0]), 1,
	         Y_1_probvector[3])

	Y_0 <- Y_tau_0
	Y_0[Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 0 & S_star_0 == 0] <- 
	  rbinom(length(Y_0[Y_0 == 0 & S_star_1 == 0 & S_star_0 == 0]), 1,
	         Y_0_probvector[1])
	Y_0[Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 1 & S_star_0 == 0] <- 
	  rbinom(length(Y_0[Y_0 == 0 & S_star_1 == 1 & S_star_0 == 0]), 1,
	         Y_0_probvector[3])
        
	# Observed data (S_star == S^tau in manuscript)
	Z <- rbinom(n, 1, 0.5)
	Y <- Y_1*Z + Y_0*(1 - Z)
	Y_tau <- Y_tau_1*Z + Y_tau_0*(1 - Z)
	S_star <- S_star_1*Z + S_star_0*(1 - Z)
       
	# Case-cohort sampling
	R <- rbinom(n, 1, ccprob)
	S_star[R == 0 & Y == 0] <- NA

	# Output observed data frame
	data.frame(Z, Y, Y_tau, S_star, ID = 1:n)
}