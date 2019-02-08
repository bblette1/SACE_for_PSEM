##########################################################################
# analyze_data()
# Description: Function to analyze data with full cohort sampling under
# No Harm Early Risk assumption set

analyze_data <- function(data, brange0, brange1, brange2, brange3) {
  ########################################################################
  # Arguments:
  #      data: data frame containing the following variables
  #               Z     : indicator of treatment
  #               Y     : indicator of outcome
  #               Y_tau : indicator of early outcome
  #               S_star: intermediate biomarker value
  #   brange*: (2 x 1) vector containing the specified lower and upper 
  #            bounds of the range for sensitivity parameter beta_*
  
  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Y_tau
  S_star <- data$S_star

  S_0 <- 1*(Y_tau==0 & S_star==0)
  S_1 <- 1*(Y_tau==0 & S_star==1)

  # Estimate identifiable parameters risk_0 and p(1, 1)
  # risk_0 = P(Y(0) = 1 | Y^tau(0) = Y^tau(1) = 0)
  f1 <- function(x) { sum((1 - Y_tau)*(1 - Z)*(Y - x)) }
  risk_0 <- uniroot(f1, c(0, 1))$root
  # p(1, 1) = 0 under NHER assumptions
  p_11 <- 0
  				

  # Estimate risk_0(0, 0), risk_1(0, 0), p(1, 0), and risk_1 using SACE
  # methods. These are the partially identifiable terms under NHER
  # Estimate p(1, 0)
  f11 <- function(x) {
    sum((1 - Y_tau)*Z*(1*(!is.na(S_star) & S_star == 1) - x)) /
      length(Y[Y_tau == 0 & Z == 1])
  }
  p_1 <- uniroot(f11, c(0, 1))$root
  
  fv <- function(x) { sum(x*(1 - Y_tau)*Z - (1 - Y_tau)*(1 - Z)) }
  vhat <- uniroot(fv, c(0, 1))$root

  p10func1 <- function(x) {
	  theta5 <- 1 / ( 1 + exp(min(brange3))*(1 - x) / x )
	  x*vhat + theta5*(1 - vhat) - p_1
  }
  p_10_first <- uniroot(p10func1, c(0, 1))$root
  theta5_first <- (p_1 - p_10_first*vhat) / (1 - vhat)
  p_00_first <- 1 - p_10_first

  p10func2 <- function(x) {
	  theta5 <- 1 / ( 1 + exp(max(brange3))*(1 - x) / x )
	  x*vhat + theta5*(1 - vhat) - p_1
  }
  p_10_second <- uniroot(p10func2, c(0, 1))$root
  theta5_second <- (p_1 - p_10_second*vhat) / (1 - vhat)
  p_00_second <- 1 - p_10_second

  # Estimate risk_0(0, 0)
  f7 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_first <- uniroot(f7, c(0, 1))$root
  risk_0_10_first <- (risk_0 - risk_0_00_first*p_00_first) / p_10_first

  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_second <- uniroot(f8, c(0, 1))$root
  risk_0_10_second <- (risk_0 - risk_0_00_second*p_00_second) / p_10_second

  f7 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_third <- uniroot(f7, c(0, 1))$root
  risk_0_10_third <- (risk_0 - risk_0_00_third*p_00_second) / p_10_second

  f8 <- function(x){
    risk_0_10 <- 1 / ( 1 + exp(max(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_fourth <- uniroot(f8, c(0,1))$root
  risk_0_10_fourth <- (risk_0 - risk_0_00_fourth*p_00_first) / p_10_first

  # Estimate risk_1(0, 0)
  row1 <- mean(Y[Z == 1 & Y_tau == 0 & S_star == 0])
  theta11 <- 1 - mean(Y_tau[Z == 0])
  theta12 <- length(Y_tau[Y_tau == 0 & S_star == 0 & Z == 1]) /
             length(Y_tau[Z == 1])
  theta13_first <- p_00_first*theta11 / theta12
  theta13_second <- p_00_second*theta11 / theta12

  f7 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    row1 - (1 - theta13_first)*risk_1_0star - theta13_first*x
  }
  risk_1_00_first <- uniroot(f7, c(0, 1))$root
  risk_1_0star_first <- (row1 - theta13_first*risk_1_00_first) /
                        (1 - theta13_first)

  f8 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    row1 - (1 - theta13_second)*risk_1_0star - theta13_second*x
  }
  risk_1_00_second <- uniroot(f8, c(0, 1))$root
  risk_1_0star_second <- (row1 - theta13_second*risk_1_00_second) /
                         (1 - theta13_second)

  f7 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    row1 - (1 - theta13_second)*risk_1_0star - theta13_second*x
  }
  risk_1_00_third <- uniroot(f7, c(0, 1))$root
  risk_1_0star_third <- (row1 - theta13_second*risk_1_00_third) /
                        (1 - theta13_second)

  f8 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    row1 - (1 - theta13_first)*risk_1_0star - theta13_first*x
  }
  risk_1_00_fourth <- uniroot(f8, c(0, 1))$root
  risk_1_0star_fourth <- (row1 - theta13_first*risk_1_00_fourth) /
                         (1 - theta13_first)
  
  # Estimate risk_1(1, 0)
  row2 <- mean(Y[Z == 1 & Y_tau == 0 & S_star == 1])
  theta17 <- length(Y_tau[Y_tau == 0 & S_star == 1 & Z == 1]) /
             length(Y_tau[Z == 1])
  theta18_first <- p_10_first*theta11 / theta17
  theta18_second <- p_10_second*theta11 / theta17

  f7 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(min(brange2))*(1 - x) / x )
    row2 - (1 - theta18_first)*risk_1_1star - theta18_first*x
  }
  risk_1_10_first <- uniroot(f7, c(0,1))$root
  risk_1_1star_first <- (row2 - theta18_first*risk_1_10_first) /
                        (1 - theta18_first)

  f8 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(max(brange2))*(1 - x) / x )
    row2 - (1 - theta18_second)*risk_1_1star - theta18_second*x
  }
  risk_1_10_second <- uniroot(f8, c(0, 1))$root
  risk_1_1star_second <- (row2 - theta18_second*risk_1_10_second) /
                         (1 - theta18_second)

  f7 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(min(brange2))*(1 - x) / x )
    row2 - (1 - theta18_second)*risk_1_1star - theta18_second*x
  }
  risk_1_10_third <- uniroot(f7, c(0, 1))$root
  risk_1_1star_third <- (row2 - theta18_second*risk_1_10_third) /
                        (1 - theta18_second)

  f8 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(max(brange2))*(1 - x) / x )
    row2 - (1 - theta18_first)*risk_1_1star - theta18_first*x
  }
  risk_1_10_fourth <- uniroot(f8, c(0, 1))$root
  risk_1_1star_fourth <- (row2 - theta18_first*risk_1_10_fourth) /
                         (1 - theta18_first)

  # Calculate CEP(1,0) - CEP(0,0) at the 16 boundary points of beta_0,
  # beta_1, beta_2, and beta_3
  # Take the Ignorance Interval to be the max and min of these estimates
  # II7 should be the upper bound, and II2 should be the lower bound
  II1  <- risk_1_10_first  - risk_0_10_first  - (risk_1_00_first  - risk_0_00_first)   
  II2  <- risk_1_10_first  - risk_0_10_first  - (risk_1_00_fourth - risk_0_00_first)   
  II3  <- risk_1_10_first  - risk_0_10_fourth - (risk_1_00_first  - risk_0_00_fourth)   
  II4  <- risk_1_10_first  - risk_0_10_fourth - (risk_1_00_fourth - risk_0_00_fourth)   
  II5  <- risk_1_10_fourth - risk_0_10_first  - (risk_1_00_first  - risk_0_00_first)   
  II6  <- risk_1_10_fourth - risk_0_10_first  - (risk_1_00_fourth - risk_0_00_first)   
  II7  <- risk_1_10_fourth - risk_0_10_fourth - (risk_1_00_first  - risk_0_00_fourth)   
  II8  <- risk_1_10_fourth - risk_0_10_fourth - (risk_1_00_fourth - risk_0_00_fourth)   
  II9  <- risk_1_10_second - risk_0_10_second - (risk_1_00_second - risk_0_00_second)   
  II10 <- risk_1_10_second - risk_0_10_second - (risk_1_00_third  - risk_0_00_second)   
  II11 <- risk_1_10_second - risk_0_10_third  - (risk_1_00_second - risk_0_00_third)   
  II12 <- risk_1_10_second - risk_0_10_third  - (risk_1_00_third  - risk_0_00_third)   
  II13 <- risk_1_10_third  - risk_0_10_second - (risk_1_00_second - risk_0_00_second)   
  II14 <- risk_1_10_third  - risk_0_10_second - (risk_1_00_third  - risk_0_00_second)   
  II15 <- risk_1_10_third  - risk_0_10_third  - (risk_1_00_second - risk_0_00_third)   
  II16 <- risk_1_10_third  - risk_0_10_third  - (risk_1_00_third  - risk_0_00_third)   

  II <- c(II1, II2, II3, II4, II5, II6, II7, II8, II9, II10, II11, II12,
          II13, II14, II15, II16)
  II_low <- min(II)
  II_up <- max(II)
  
  # Output estimates
  data.frame(risk_0, vhat, p_1, p_10_first, theta5_first, p_10_second,
             theta5_second, risk_0_00_first, risk_0_10_first,
             risk_0_00_fourth, risk_0_10_fourth, row1, theta11, theta12,
             theta13_first, risk_1_00_first, risk_1_0star_first,
             risk_1_00_fourth, risk_1_0star_fourth, row2, theta17,
             theta18_first, risk_1_10_first, risk_1_1star_first,
             risk_1_10_fourth, risk_1_1star_fourth, II_low, II_up)

}
  
##########################################################################
# analyze_data_cc()
# Description: Function to analyze data with case-cohort sampling under
# No Harm Early Risk assumption set

analyze_data_cc <- function(data,brange0,brange1,brange2,brange3) {
  ########################################################################
  # Arguments:
  #      data: data frame containing the following variables
  #               Z     : indicator of treatment
  #               Y     : indicator of outcome
  #               Y_tau : indicator of early outcome
  #               S_star: intermediate biomarker value
  #               R     : indicator of measurement of intermediate biomarker
  #   brange*: (2 x 1) vector containing the specified lower and upper 
  #            bounds of the range for sensitivity parameter beta_*

  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Y_tau
  S_star <- data$S_star
  R <- data$R

  # Estimate probability control has S_star observed: P(R = 1 | Y = 0)
  f0 <- function(x) { sum((1 - Y)*(R - x)) }
  pi_hat <- uniroot(f0, c(0, 1))$root

  # weights
  W <- 1 / pi_hat*(1 - Y)*R + Y

  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star
    
  # Estimate identifiable parameters risk_0, p(1, 1)
  # risk_0 = P(Y(0) = 1 | Y^tau(0) = Y^tau(1) = 0)
  f1 <- function(x) { sum((1 - Y_tau)*(1 - Z)*(Y - x)) }
  risk_0 <- uniroot(f1, c(0, 1))$root
  # p(1, 1) = 0 under NHER assumptions
  p_11 <- 0

  # Estimate risk_0(0, 0), risk_1(0, 0), p(1, 0), and risk_1 using SACE 
  # methods. These are the partially identifiable terms under NHER
  # Estimate p(1, 0)
  f11 <- function(x) {
    sum((1 - Y_tau)*Z*(1*(!is.na(S_star) & S_star == 1) - x)*W) /
      length(Y[Y_tau == 0 & Z == 1])
  }
  p_1 <- uniroot(f11, c(0, 1))$root
  
  fv <- function(x) { sum(x*(1 - Y_tau)*Z - (1 - Y_tau)*(1 - Z)) }
  vhat <- uniroot(fv, c(0, 1))$root

  p10func1 <- function(x) {
	  theta5 <- 1 / ( 1 + exp(min(brange3))*(1 - x) / x )
	  x*vhat + theta5*(1 - vhat) - p_1
  }
  p_10_first <- uniroot(p10func1, c(0, 1))$root
  theta5_first <- (p_1 - p_10_first*vhat) / (1 - vhat)
  p_00_first <- 1 - p_10_first

  p10func2 <- function(x) {
  	theta5 <- 1 / ( 1 + exp(max(brange3))*(1 - x) / x )
  	x*vhat + theta5*(1 - vhat) - p_1
  }
  p_10_second <- uniroot(p10func2, c(0, 1))$root
  theta5_second <- (p_1 - p_10_second*vhat) / (1 - vhat)
  p_00_second <- 1 - p_10_second

  # Estimate risk_0(0, 0)
  f7 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_first <- uniroot(f7, c(0, 1))$root
  risk_0_10_first <- (risk_0 - risk_0_00_first*p_00_first) / p_10_first

  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_second <- uniroot(f8, c(0, 1))$root
  risk_0_10_second <- (risk_0 - risk_0_00_second*p_00_second) / p_10_second

  f7 <- function(x){
    risk_0_10 <- 1 / ( 1 + exp(min(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_third <- uniroot(f7, c(0, 1))$root
  risk_0_10_third <- (risk_0 - risk_0_00_third*p_00_second) / p_10_second

  f8 <- function(x){
    risk_0_10 <- 1 / ( 1 + exp(max(brange0))*(1 - x) / x )
    risk_0 - x*p_00_second - risk_0_10*p_10_second
  }
  risk_0_00_fourth <- uniroot(f8, c(0, 1))$root
  risk_0_10_fourth <- (risk_0 - risk_0_00_fourth*p_00_first) / p_10_first

  # Estimate risk_1(0, 0)
  f5 <- function(x) { sum((1 - Y_tau)*Z*(1 - S_star)*W*(Y - x)) }
  row1 <- uniroot(f5, c(0, 1))$root
  theta11 <- 1-mean(Y_tau[Z == 0])
  f6 <- function(x) { sum(Z*((1 - S_star)*(1 - Y_tau)*W - x)) }
  theta12 <- uniroot(f6, c(0, 1))$root
  theta13_first <- p_00_first*theta11 / theta12
  theta13_second <- p_00_second*theta11 / theta12

  f7 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    row1 - (1 - theta13_first)*risk_1_0star - theta13_first*x
  }
  risk_1_00_first <- uniroot(f7, c(0, 1))$root
  risk_1_0star_first <- (row1 - theta13_first*risk_1_00_first) /
                        (1 - theta13_first)

  f8 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    row1 - (1 - theta13_second)*risk_1_0star - theta13_second*x
  }
  risk_1_00_second <- uniroot(f8, c(0, 1))$root
  risk_1_0star_second <- (row1 - theta13_second*risk_1_00_second) /
                         (1 - theta13_second)

  f7 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    row1 - (1 - theta13_second)*risk_1_0star - theta13_second*x
  }
  risk_1_00_third <- uniroot(f7, c(0, 1))$root
  risk_1_0star_third <- (row1 - theta13_second*risk_1_00_third) /
                        (1 - theta13_second)

  f8 <- function(x) {
    risk_1_0star <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    row1 - (1 - theta13_first)*risk_1_0star - theta13_first*x
  }
  risk_1_00_fourth <- uniroot(f8, c(0, 1))$root
  risk_1_0star_fourth <- (row1 - theta13_first*risk_1_00_fourth) /
                         (1 - theta13_first)

  # Estimate risk_1_10
  f5 <- function(x) { sum((1 - Y_tau)*Z*S_star*W*(Y - x)) }
  row2 <- uniroot(f5, c(0, 1))$root
  f6 <- function(x) { sum(Z*(S_star*(1 - Y_tau)*W - x)) }
  theta17 <- uniroot(f6, c(0, 1))$root
  theta18_first <- p_10_first*theta11 / theta17
  theta18_second <- p_10_second*theta11 / theta17

  f7 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(min(brange2))*(1 - x) / x )
    row2 - (1 - theta18_first)*risk_1_1star - theta18_first*x
  }
  risk_1_10_first <- uniroot(f7, c(0, 1))$root
  risk_1_1star_first <- (row2 - theta18_first*risk_1_10_first) /
                        (1 - theta18_first)

  f8 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(max(brange2))*(1 - x) / x )
    row2 - (1 - theta18_second)*risk_1_1star - theta18_second*x
  }
  risk_1_10_second <- uniroot(f8, c(0, 1))$root
  risk_1_1star_second <- (row2 - theta18_second*risk_1_10_second) /
                         (1 - theta18_second)

  f7 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(min(brange2))*(1 - x) / x )
    row2 - (1 - theta18_second)*risk_1_1star - theta18_second*x
  }
  risk_1_10_third <- uniroot(f7, c(0, 1))$root
  risk_1_1star_third <- (row2 - theta18_second*risk_1_10_third) /
                        (1 - theta18_second)

  f8 <- function(x) {
    risk_1_1star <- 1 / ( 1 + exp(max(brange2))*(1 - x) / x )
    row2 - (1 - theta18_first)*risk_1_1star - theta18_first*x
  }
  risk_1_10_fourth <- uniroot(f8, c(0, 1))$root
  risk_1_1star_fourth <- (row2 - theta18_first*risk_1_10_fourth) /
                         (1 - theta18_first)

  # Calculate CEP(1,0) - CEP(0,0) at the 16 boundary points of beta_0,
  # beta_1, beta_2, and beta_3
  # Take the Ignorance Interval to be the max and min of these estimates
  # II7 should be the upper bound, and II2 should be the lower bound
  II1  <- risk_1_10_first  - risk_0_10_first  - (risk_1_00_first  - risk_0_00_first)   
  II2  <- risk_1_10_first  - risk_0_10_first  - (risk_1_00_fourth - risk_0_00_first)   
  II3  <- risk_1_10_first  - risk_0_10_fourth - (risk_1_00_first  - risk_0_00_fourth)   
  II4  <- risk_1_10_first  - risk_0_10_fourth - (risk_1_00_fourth - risk_0_00_fourth)   
  II5  <- risk_1_10_fourth - risk_0_10_first  - (risk_1_00_first  - risk_0_00_first)   
  II6  <- risk_1_10_fourth - risk_0_10_first  - (risk_1_00_fourth - risk_0_00_first)   
  II7  <- risk_1_10_fourth - risk_0_10_fourth - (risk_1_00_first  - risk_0_00_fourth)   
  II8  <- risk_1_10_fourth - risk_0_10_fourth - (risk_1_00_fourth - risk_0_00_fourth)   
  II9  <- risk_1_10_second - risk_0_10_second - (risk_1_00_second - risk_0_00_second)   
  II10 <- risk_1_10_second - risk_0_10_second - (risk_1_00_third  - risk_0_00_second)   
  II11 <- risk_1_10_second - risk_0_10_third  - (risk_1_00_second - risk_0_00_third)   
  II12 <- risk_1_10_second - risk_0_10_third  - (risk_1_00_third  - risk_0_00_third)   
  II13 <- risk_1_10_third  - risk_0_10_second - (risk_1_00_second - risk_0_00_second)   
  II14 <- risk_1_10_third  - risk_0_10_second - (risk_1_00_third  - risk_0_00_second)   
  II15 <- risk_1_10_third  - risk_0_10_third  - (risk_1_00_second - risk_0_00_third)   
  II16 <- risk_1_10_third  - risk_0_10_third  - (risk_1_00_third  - risk_0_00_third)   

  II <- c(II1, II2, II3, II4, II5, II6, II7, II8, II9, II10, II11, II12,
          II13, II14, II15, II16)
  II_low <- min(II)
  II_up <- max(II)

  # Output estimates
  data.frame(risk_0, vhat, p_1, p_10_first, theta5_first, p_10_second,
             theta5_second, risk_0_00_first, risk_0_10_first,
             risk_0_00_fourth, risk_0_10_fourth, row1, theta11, theta12,
             theta13_first,	risk_1_00_first, risk_1_0star_first,
             risk_1_00_fourth, risk_1_0star_fourth,	row2, theta17,
             theta18_first,	risk_1_10_first, risk_1_1star_first,
             risk_1_10_fourth, risk_1_1star_fourth,	II_low, II_up, pi_hat)

}