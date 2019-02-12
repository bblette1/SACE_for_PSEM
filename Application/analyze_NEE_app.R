##########################################################################
# analyze_data_NEE()
# Description: Function to analyze data with at least one biomarker of
# interest under No Early Effect assumption set and three contrast options

analyze_data_NEE <- function(data, S_star_choice, brange,
                             contrast = "logRR") {
  ########################################################################
  # Arguments:
  #            data: data frame containing the following variables
  #                     Z     : indicator of treatment
  #                     Y     : indicator of outcome
  #                     Y_tau : indicator of early outcome
  #                     Other : intermediate biomarker values
  #   S_star_choice: column number of `data` corresponding to the
  #                  biomarker of interest
  #          brange: (2 x 1) vector containing the specified lower and upper 
  #                  bounds of the range for sensitivity parameter beta_0
  #        contrast: string of value "logRR", "Difference", or "VE"
  #                  specifying scale of estimand of interest

  # Specify estimand contrast
  h <- function(x, y) {
    ans <- log(x) - log(y)
    if (contrast=="Difference") { ans <- x - y }
    if (contrast=="VE") { ans <- 1 - (x / y) }
    return(ans)
  }

  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Y_tau
  S_star <- data[ , S_star_choice]
  R <- data$R
  S_star[R == 0] <- 0
  brange <- c(-beta, beta)
  
  # estimate probability control has S_star observed: P(R = 1 | Y = 0)
  f0 <- function(x) { sum((1 - Y)*(R - x)) }
  pi_hat <- uniroot(f0, c(0, 1))$root

  # weights
  W <- 1 / pi_hat*(1 - Y)*R + Y

  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star
        
  # Estimate identifiable parameters risk_z, p(0, 0), p(1, 0), and p(1, 1)
  # risk_0 = P(Y(0) = 1 | Y^tau(0) = Y^tau(1) = 0)
  f1 <- function(x) { sum( (1 - Y_tau)*(1 - Z)*(Y - x)) }
  risk_0 <- uniroot(f1, c(0, 1))$root
  
  # risk_1 = P(Y(1) = 1 | Y^tau(0) = Y^tau(1) = 0)
  f2 <- function(x) { sum((1 - Y_tau)*Z*(Y - x)) }
  risk_1 <- uniroot(f2, c(0, 1))$root
  
  # p(0, 0) = P(S^tau(1) = S^tau(0) = 0 | Y^tau(0) = Y^tau(1) = 0)
  # Recall Case CB: P(S^tau(0) = 0 | Y^tau(0) = Y^tau(1) = 0) = 1
  f3 <- function(x) { sum((1 - Y_tau)*Z*(1 - S_star - x)*W) }
  p_00 <- uniroot(f3, c(0, 1))$root						
  p_10 <- 1 - p_00
  # p(1, 1) = 0 under NEE scenario
  
  # Estimate identifiable parameters risk_1(0, 0) and risk_1(1, 0)
  # risk_z(s1, s0) = 
  # P(Y(z) = 1 | S^tau(1) = s1, S^tau(0) = s0, Y^tau(1) = Y^tau(0) = 0)
  f5 <- function(x) { sum(Z*S_0*(Y - x)*W) }
  risk_1_00 <- uniroot(f5, c(0, 1))$root	
  f6 <- function(x) { sum(Z*S_1*(Y - x)*W) }
  risk_1_10 <- uniroot(f6, c(0, 1))$root

  # Estimate risk_0(0, 0) using SACE method
  # This is the only partially identifiable term under NEE
  f7 <- function(x) {
  	risk_0_10 <- 1 / ( 1 + exp(min(brange))*(1 - x) / x )
  	risk_0 - x*p_00 - risk_0_10*p_10
	}
  risk_0_00_first <- uniroot(f7, c(0, 1))$root
  risk_0_10_first <- (risk_0 - risk_0_00_first*p_00) / p_10

  f8 <- function(x) {
  	risk_0_10 <- 1 / ( 1 + exp(max(brange))*(1 - x) / x )
  	risk_0 - x*p_00 - risk_0_10*p_10
	}
  risk_0_00_second <- uniroot(f8, c(0, 1))$root
  risk_0_10_second <- (risk_0 - risk_0_00_second*p_00) / p_10
  
  # Calculate CEP(1,0) - CEP(0,0) at the 2 boundary points of beta_0
  # Take the Ignorance Interval to be the max and min of these estimates
  II1 <- h(risk_1_10, risk_0_10_first) - h(risk_1_00, risk_0_00_first)
  II2 <- h(risk_1_10, risk_0_10_second) - h(risk_1_00, risk_0_00_second)
  II_low <- min(II1, II2)
  II_up <- max(II1, II2)
  
  # Output estimates
  data.frame(risk_0, risk_1, p_00, risk_1_00, risk_0_00_first,
             risk_0_00_second, risk_0_10_first, risk_0_10_second, II_low,
             II_up, pi_hat, risk_1_10)

}