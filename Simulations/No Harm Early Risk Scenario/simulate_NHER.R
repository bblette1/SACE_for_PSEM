###########################################################################
# simulate_data() 
# Description: Function that simulates trial data under No Early Harm
# assumption set

simulate_data <- function(Y_tau_probvector, S_star_probvector_ytaueq00,
                          S_star_probvector_ytaueq01, Y_1_probvector,
                          Y_0_probvector, n, ccprob) {
  #########################################################################
  # Arguments                  :
  #  Y_tau_probvector          : (4 x 1) vector containing 
  #                              P(Y^tau(1) = a, Y^tau(0) = b) for (a, b) in 
  #                              {(0, 0), (0, 1), (1, 0), (1, 1)}
  #  S_star_probvector_ytaueq00: same input format as above for
  #                              P(S^tau(1) = a, S^tau(0) = b | 
  #                                              Y^tau(1) = 0, Y^tau(0) = 0)
  #  S_star_probvector_ytaueq01: same input format as above for
  #                              P(S^tau(1) = a, S^tau(0) = b | 
  #                                              Y^tau(1) = 0, Y^tau(0) = 1)
  #  Y_1_probvector            : same input format as above for
  #                              P(Y(1) = 1 | S^tau(1) = a, S^tau(0) = b,
  #                                           Y^tau(0) = 0, Y^tau(0) = 0)
  #  Y_0_probvector            : same input format as above for
  #                              P(Y(0) = 1 | S^tau(1) = a, S^tau(0) = b,
  #                                           Y^tau(0) = 0, Y^tau(0) = 0)
  #  n                         : sample size
  #  ccprob                    : case-cohort sampling probability:
  #                              P(R = 1 | Y = 0)
  
  # Generate Y^tau(1) and Y^tau(0)
	multi <- rmultinom(1, n, Y_tau_probvector)
	Y_tau_1 <- c(rep(0, multi[1]), rep(0, multi[2]),
	             rep(1, multi[3]), rep(1, multi[4]))
	Y_tau_0 <- c(rep(0, multi[1]), rep(1, multi[2]),
	             rep(0, multi[3]), rep(1, multi[4]))
      
	# Generate S^tau(1) and S^tau(0)
	S_star_1 <- rep(NA, n)
	S_star_0 <- rep(NA, n)
	multi3 <- rmultinom(1, length(Y_tau_0[Y_tau_0 == 0 & Y_tau_1 == 0]),
	                    S_star_probvector_ytaueq00)
	S_star_1[Y_tau_0 == 0 & Y_tau_1 == 0] <- c(rep(0, multi3[1]),
	                                           rep(0, multi3[2]),
	                                           rep(1, multi3[3]),
	                                           rep(1, multi3[4]))
	S_star_1[Y_tau_0 == 1 & Y_tau_1 == 0] <-
	  rbinom(length(S_star_1[Y_tau_0 == 1 & Y_tau_1 == 0]), 1,
	         S_star_probvector_ytaueq01[2])
	S_star_0[Y_tau_0 == 0 & Y_tau_1 == 0] <- c(rep(0, multi3[1]),
	                                           rep(1, multi3[2]),
	                                           rep(0, multi3[3]),
	                                           rep(1, multi3[4]))
        
	# Generate Y(1) and Y(0)
	Y_1 <- Y_tau_1
	Y_0 <- Y_tau_0
	index <- Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 0 & S_star_0 == 0
	Y_1[index] <- rbinom(length(Y_1[index]), 1, Y_1_probvector[1])
	Y_0[index] <- rbinom(length(Y_0[index]), 1, Y_0_probvector[1])

	index <- Y_tau_1 == 0 & Y_tau_0 == 0 & S_star_1 == 1 & S_star_0 == 0
	Y_1[index] <- rbinom(length(Y_1[index]), 1, Y_1_probvector[3])
	Y_0[index] <- rbinom(length(Y_0[index]), 1, Y_0_probvector[3])

	index <- (Y_tau_1 == 0 & Y_tau_0 == 1 & S_star_1 == 0 &
	          is.na(S_star_1) == F & is.na(S_star_0) == T)
	Y_1[index] <- rbinom(length(Y_1[index]), 1, Y_1_probvector[1])

	index <- (Y_tau_1 == 0 & Y_tau_0 == 1 & S_star_1 == 1 &
	          is.na(S_star_1) == F & is.na(S_star_0) == T)
	Y_1[index] <- rbinom(length(Y_1[index]), 1, Y_1_probvector[3])

	# Observed data (S_star == S^tau in manuscript)
	Z <- rbinom(n, 1, 0.5)
	Y <- Y_1*Z + Y_0*(1 - Z)
	Y_tau <- Y_tau_1*Z + Y_tau_0*(1 - Z)
	S_star <- S_star_1
	S_star[Z == 0] <- S_star_0[Z == 0]
       
	# case cohort sampling
	R <- rbinom(n, 1, ccprob)
	S_star[R == 0 & Y == 0] <- NA

	# Output observed data frame
	data.frame(Z, Y, Y_tau, S_star, ID=1:n)
}