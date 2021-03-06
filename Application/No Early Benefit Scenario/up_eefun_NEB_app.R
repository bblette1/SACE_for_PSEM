ScenarioC_up_eefun_minmin <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1]),							# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2]),								# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(min(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9]),								# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(min(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR")  # CEP(1,0) - CEP(0,0) = theta[12]                  
    ))
  }
  
}

ScenarioC_up_eefun_minmax <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1]),							# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2]),								# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(min(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9]),								# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(max(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR")  # CEP(1,0) - CEP(0,0) = theta[12]                  
    ))
  }
  
}

ScenarioC_up_eefun_maxmin <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1]),							# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2]),								# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(max(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9]),								# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(min(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR")  # CEP(1,0) - CEP(0,0) = theta[12]                  
    ))
  }
  
}

ScenarioC_up_eefun_maxmax <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1]),							# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2]),								# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(max(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9]),								# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(max(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR")  # CEP(1,0) - CEP(0,0) = theta[12]                  
    ))
  }
  
}







ScenarioC_up_eefun_cc_minmin <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1])*(1/theta[13]*(1-Y)*R+Y),			# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2])*(1/theta[13]*(1-Y)*R+Y),				# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(min(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9])*(1/theta[13]*(1-Y)*R+Y),				# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(min(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR"),  # CEP(1,0) - CEP(0,0) = theta[12]                  
                  
                  (1-Y)*(R-theta[13])														# pi_hat = theta[13]
                  
    ))
  }
  
}

ScenarioC_up_eefun_cc_minmax <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1])*(1/theta[13]*(1-Y)*R+Y),			# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2])*(1/theta[13]*(1-Y)*R+Y),				# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(min(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9])*(1/theta[13]*(1-Y)*R+Y),				# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(max(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR"),  # CEP(1,0) - CEP(0,0) = theta[12]                  
                  
                  (1-Y)*(R-theta[13])														# pi_hat = theta[13]
                  
    ))
  }
  
}

ScenarioC_up_eefun_cc_maxmin <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1])*(1/theta[13]*(1-Y)*R+Y),			# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2])*(1/theta[13]*(1-Y)*R+Y),				# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(max(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9])*(1/theta[13]*(1-Y)*R+Y),				# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(min(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR"),  # CEP(1,0) - CEP(0,0) = theta[12]                  
                  
                  (1-Y)*(R-theta[13])														# pi_hat = theta[13]
                  
    ))
  }
  
}

ScenarioC_up_eefun_cc_maxmax <- function(data, beta0range, beta1range, contrast) {
  
  function(theta) {
    with(data, c( (1-Ytau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1])*(1/theta[13]*(1-Y)*R+Y),			# risk_1_00 = theta[1]
                  (1-Ytau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2])*(1/theta[13]*(1-Y)*R+Y),				# risk_1_10 = theta[2]
                  
                  Z*(1-Ytau-theta[3]),														# P(Ytau(1) = 0) = theta[3]
                  (1-Z)*(1-Ytau-theta[4]),													# P(Ytau(0) = 0) = theta[4]
                  theta[5] - theta[3]/theta[4],													# P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
                  (1-Z)*(1-Ytau)*(Y-theta[6]),													# P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
                  theta[6] - theta[7]*theta[5] - theta[8]*(1-theta[5]),									# risk_0 = theta[7]
                  exp(max(beta0range))*theta[8]/(1-theta[8]) - theta[7]/(1-theta[7]),						# risk_new = theta[8]
                  
                  (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9])*(1/theta[13]*(1-Y)*R+Y),				# p_10 = theta[9] 
                  theta[7] - theta[10]*(1-theta[9]) - theta[11]*theta[9],								# risk_0_00 = theta[10]
                  exp(max(beta1range))*theta[11]/(1-theta[11]) - theta[10]/(1-theta[10]),						# risk_0_10 = theta[11]
                  
                  ( theta[12] - (1- theta[2]/theta[11]) + (1- theta[1]/theta[10]) )*(contrast == "VE") +
                  ( theta[12] - (theta[2] - theta[11]) + (theta[1] - theta[10]) )*(contrast == "Difference") +
                  ( theta[12] - (log(theta[2]) - log(theta[11])) + (log(theta[1]) - log(theta[10])) )*(contrast == "logRR"),  # CEP(1,0) - CEP(0,0) = theta[12]                  
                  
                  (1-Y)*(R-theta[13])														# pi_hat = theta[13]
                  
    ))
  }
  
}