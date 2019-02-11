# Estimating equation vector for ignorance interval upper bound under NEH
# Used for M estimation variance calculation

# Full-cohort estimating equations
up_eefun <- function(data, beta0range, beta1range, beta2range, beta3range) {

  function(theta) {
    with(data, c( (1-Y_tau)*(1-Z)*(Y-theta[1]),												                    # risk_0 = theta[1]
              		theta[2]*(1-Y_tau)*Z - (1-Y_tau)*(1-Z),											            # p_ytau0eq0_giv_ytau1eq0 = v = theta[2]
              		(1-Y_tau)*Z*((!is.na(S_star) & S_star==1)-theta[3]),									  # p_1 = theta[3]
              		theta[4]*theta[2] + theta[5]*(1-theta[2]) - theta[3],									  # p_10_first = theta[4]
              		exp(min(beta3range))*theta[5]/(1-theta[5]) - theta[4]/(1-theta[4]),			# P(S_star(1) = 1 | Y_tau(1) = 0, Y_tau(0) = 1)_first = theta[5]
              
              		theta[1] - theta[6]*(1-theta[4]) - theta[7]*theta[4],									  # risk_0_00 = theta[6] 
              		exp(max(beta0range))*theta[7]/(1-theta[7]) - theta[6]/(1-theta[6]),			# risk_0_10 = theta[7]
              
              		Z*(1-Y_tau)*(!is.na(S_star) & S_star==0)*(Y-theta[8]),								  # row1 = theta[8]
              		(1-Z)*(1-Y_tau-theta[9]),													                      # P(Y_tau(0) = 0) = theta[9]
              		(1-Y_tau)*(!is.na(S_star) & S_star==0)*Z - Z*theta[10],								  # P(Y_tau(1) = 0, S_star(1) = 0) = theta[10]
              		(1-theta[4])*theta[9] - theta[11]*theta[10],										        # P(Y_tau(0) = 0, S_star(0) = 0 | Y_tau(1) = 0, S*(1) = 0) = theta[11]
              		theta[8] - theta[12]*theta[11] - theta[13]*(1-theta[11]),								# risk_1_00 = theta[12]
              		exp(min(beta1range))*theta[13]/(1-theta[13]) - theta[12]/(1-theta[12]),	# risk_1_0star = theta[13]
              
              		Z*(1-Y_tau)*(!is.na(S_star) & S_star==1)*(Y-theta[14]),							  	# row2 = theta[14]
              		(1-Y_tau)*(!is.na(S_star) & S_star==1)*Z - Z*theta[15],								  # P(Y_tau(1) = 0, S_star(1) = 1) = theta[15]
              		theta[4]*theta[11] - theta[16]*theta[15],										            # P(Y_tau(0) = 0, S_star(0) = 0 | Y_tau(1) = 0, S*(1) = 1) = theta[16]
              		theta[14] - theta[17]*theta[16] - theta[18]*(1-theta[16]),							# risk_1_10 = theta[17]
              		exp(max(beta2range))*theta[18]/(1-theta[18]) - theta[17]/(1-theta[17]),	# risk_1_1star = theta[18]
              
              		theta[19] - theta[17] + theta[7] + (theta[12] - theta[6])								# CEP(1,0) - CEP(0,0) = theta[19]

		))
	}

}


# Case-cohort sampling estimating equations
up_eefun_cc <- function(data, beta0range, beta1range, beta2range, beta3range) {

  function(theta) {
    with(data, c( (1-Y_tau)*(1-Z)*(Y-theta[1]),												                    # risk_0 = theta[1]
              		theta[2]*(1-Y_tau)*Z - (1-Y_tau)*(1-Z),											            # p_ytau0eq0_giv_ytau1eq0 = v = theta[2]
              		(1-Y_tau)*Z*(S_star-theta[3])*(1/theta[20]*(1-Y)*R+Y),								  # p_1 = theta[3]
              		theta[4]*theta[2] + theta[5]*(1-theta[2]) - theta[3],									  # p_10_first = theta[4]
              		exp(min(beta3range))*theta[5]/(1-theta[5]) - theta[4]/(1-theta[4]),			# P(S_star(1) = 1 | Y_tau(1) = 0, Y_tau(0) = 1)_first = theta[5]
              
              		theta[1] - theta[6]*(1-theta[4]) - theta[7]*theta[4],									  # risk_0_00 = theta[6] 
              		exp(max(beta0range))*theta[7]/(1-theta[7]) - theta[6]/(1-theta[6]),			# risk_0_10 = theta[7]
              
              		Z*(1-Y_tau)*(1-S_star)*(1/theta[20]*(1-Y)*R+Y)*(Y-theta[8]),						# row1 = theta[8]
              		(1-Z)*(1-Y_tau-theta[9]),													                      # P(Y_tau(0) = 0) = theta[9]
              		Z*((1-S_star)*(1-Y_tau)*(1/theta[20]*(1-Y)*R+Y) - theta[10]),						# P(Y_tau(1) = 0, S_star(1) = 0) = theta[10]
              		(1-theta[4])*theta[9] - theta[11]*theta[10],										        # P(Y_tau(0) = 0, S_star(0) = 0 | Y_tau(1) = 0, S*(1) = 0) = theta[11]
              		theta[8] - theta[12]*theta[11] - theta[13]*(1-theta[11]),								# risk_1_00 = theta[12]
              		exp(min(beta1range))*theta[13]/(1-theta[13]) - theta[12]/(1-theta[12]),	# risk_1_0star = theta[13]
              
              		Z*(1-Y_tau)*S_star*(1/theta[20]*(1-Y)*R+Y)*(Y-theta[14]),								# row2 = theta[14]
              		Z*(S_star*(1-Y_tau)*(1/theta[20]*(1-Y)*R+Y) - theta[15]),								# P(Y_tau(1) = 0, S_star(1) = 1) = theta[15]
              		theta[4]*theta[11] - theta[16]*theta[15],											          # P(Y_tau(0) = 0, S_star(0) = 0 | Y_tau(1) = 0, S*(1) = 1) = theta[16]
              		theta[14] - theta[17]*theta[16] - theta[18]*(1-theta[16]),							# risk_1_10 = theta[17]
              		exp(max(beta2range))*theta[18]/(1-theta[18]) - theta[17]/(1-theta[17]),	# risk_1_1star = theta[18]
              
              		theta[19] - theta[17] + theta[7] + (theta[12] - theta[6]),							# CEP(1,0) - CEP(0,0) = theta[19]
              
              		(1-Y)*(R-theta[20])														                          # pi_hat = theta[20]		

		))
	}

}