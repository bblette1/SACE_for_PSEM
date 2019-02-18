ScenarioB_low_eefun_min <- function(data, beta0range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*(1-Z)*(Y-theta[1]),						# risk_0 = theta[1]
                  (1-Ytau)*S_star*Z*(Y-theta[2]),							# risk_1_10 = theta[2]
                  (1-Ytau)*Z*(1-S_star-theta[3]),						# p_00 = theta[3]
                  Z*(1-Ytau)*(1-S_star)*(Y-theta[4]), 					# risk_1_00 = theta[4]
                  theta[1] - theta[5]*theta[3] - theta[6]*(1-theta[3]),			# risk_0_00 = theta[5] 
                  exp(min(beta0range))*theta[6]/(1-theta[6]) - theta[5]/(1-theta[5]), # risk_0_10 = theta[6]
                  ( theta[7] - (1- theta[4]/theta[5]) )*(contrast == "VE") +
                  ( theta[7] - (theta[4] - theta[5]) )*(contrast == "Difference") +
                  ( theta[7] - (log(theta[4]) - log(theta[5])) )*(contrast == "logRR")
    ))
  }
  
}

ScenarioB_low_eefun_max <- function(data, beta0range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*(1-Z)*(Y-theta[1]),						# risk_0 = theta[1]
                  (1-Ytau)*S_star*Z*(Y-theta[2]),							# risk_1_10 = theta[2]
                  (1-Ytau)*Z*(1-S_star-theta[3]),						# p_00 = theta[3]
                  Z*(1-Ytau)*(1-S_star)*(Y-theta[4]), 					# risk_1_00 = theta[4]
                  theta[1] - theta[5]*theta[3] - theta[6]*(1-theta[3]),			# risk_0_00 = theta[5] 
                  exp(max(beta0range))*theta[6]/(1-theta[6]) - theta[5]/(1-theta[5]), # risk_0_10 = theta[6]
                  ( theta[7] - (1- theta[4]/theta[5]) )*(contrast == "VE") +
                  ( theta[7] - (theta[4] - theta[5]) )*(contrast == "Difference") +
                  ( theta[7] - (log(theta[4]) - log(theta[5])) )*(contrast == "logRR")
    ))
  }
  
}


ScenarioB_low_eefun_cc_min <- function(data, beta0range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*(1-Z)*(Y-theta[1]),						# risk_0 = theta[1]
                  (1-Ytau)*S_star*Z*(1/theta[8]*(1-Y)*R+Y)*(Y-theta[2]),							# risk_1_10 = theta[2]
                  (1-Ytau)*Z*(1/theta[8]*(1-Y)*R+Y)*(1-S_star-theta[3]),		# p_00 = theta[3]
                  Z*(1-Ytau)*(1-S_star)*(1/theta[8]*(1-Y)*R+Y)*(Y-theta[4]), 		# risk_1_00 = theta[4]
                  theta[1] - theta[5]*theta[3] - theta[6]*(1-theta[3]),			# risk_0_00 = theta[5] 
                  exp(min(beta0range))*theta[6]/(1-theta[6]) - theta[5]/(1-theta[5]), # risk_0_10 = theta[6]
                  ( theta[7] - (1- theta[4]/theta[5]) )*(contrast == "VE") +
                  ( theta[7] - (theta[4] - theta[5]) )*(contrast == "Difference") +
                  ( theta[7] - (log(theta[4]) - log(theta[5])) )*(contrast == "logRR"),
                  (1-Y)*(R-theta[8])
    ))
  }
  
}

ScenarioB_low_eefun_cc_max <- function(data, beta0range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*(1-Z)*(Y-theta[1]),						# risk_0 = theta[1]
                  (1-Ytau)*S_star*Z*(1/theta[8]*(1-Y)*R+Y)*(Y-theta[2]),							# risk_1_10 = theta[2]
                  (1-Ytau)*Z*(1/theta[8]*(1-Y)*R+Y)*(1-S_star-theta[3]),		# p_00 = theta[3]
                  Z*(1-Ytau)*(1-S_star)*(1/theta[8]*(1-Y)*R+Y)*(Y-theta[4]), 		# risk_1_00 = theta[4]
                  theta[1] - theta[5]*theta[3] - theta[6]*(1-theta[3]),			# risk_0_00 = theta[5] 
                  exp(max(beta0range))*theta[6]/(1-theta[6]) - theta[5]/(1-theta[5]), # risk_0_10 = theta[6]
                  ( theta[7] - (1- theta[4]/theta[5]) )*(contrast == "VE") +
                  ( theta[7] - (theta[4] - theta[5]) )*(contrast == "Difference") +
                  ( theta[7] - (log(theta[4]) - log(theta[5])) )*(contrast == "logRR"),
                  (1-Y)*(R-theta[8])
    ))
  }
  
}
