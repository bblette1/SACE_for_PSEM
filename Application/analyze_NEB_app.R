##########################################################################
# analyze_data_NEB()
# Description: Function to analyze data with at least one biomarker of
# interest under No Early Benefit assumption set and three contrast options

analyze_data_NEB <- function(data, S_star_choice, brange,
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
  #                  bounds of the range for sensitivity parameters beta_0,1
  #        contrast: string of value "logRR", "Difference", or "VE"
  #                  specifying scale of estimand of interest
  
  # Specify estimand contrast
  h <- function(x, y) {
    ans <- log(x) - log(y)
    if (contrast == "Difference") { ans <- x - y }
    if (contrast == "VE") { ans <- 1 - (x / y) }
    return(ans)
  }
  
  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Y_tau
  S_star <- data[ , S_star_choice]
  R <- data$R
  S_star[R == 0] <- 0
  brange0 <- brange1 <- brange
  
  # Estimate probability control has S_star observed: P(R = 1 | Y = 0)
  f0 <- function(x) { sum((1 - Y)*(R - x)) }
  pi_hat <- uniroot(f0, c(0, 1))$root
  
  # weights
  W <- 1 / pi_hat*(1 - Y)*R + Y
  
  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star
  
  # Estimate p(1, 0)
  f3 <- function(x) {
    sum((1 - Y_tau)*Z*(S_star - x)*W)
  }
  p_10 <- uniroot(f3, c(0, 1))$root
  p_00 <- 1 - p_10
  
  # Estimate identifiable parameters risk_1(0, 0) and risk_1(1, 0)
  # risk_z(s1, s0) = 
  # P(Y(z) = 1 | S^tau(1) = s1, S^tau(0) = s0, Y^tau(1) = Y^tau(0) = 0)
  f4 <- function(x) {
    sum((1 - Y_tau)*Z*1*(!is.na(S_star) & S_star == 0)*(Y - x)*W)
  }
  risk_1_00 <- uniroot(f4, c(0, 1))$root
  f5 <- function(x) {
    sum((1 - Y_tau)*Z*1*(!is.na(S_star) & S_star == 1)*(Y - x)*W)
  }
  risk_1_10 <- uniroot(f5, c(0, 1))$root
  
  # Estimate risk_0 with SACE method (first partially indentifiable term)
  mix_1 <- 1 - mean(Y_tau[Z == 1])
  mix_2 <- 1 - mean(Y_tau[Z == 0])
  mix <- mix_1 / mix_2
  prob <- mean(Y[Y_tau == 0 & Z == 0])
  
  f6 <- function(x) {
    risk_new <- 1 / ( 1 + exp(min(brange0))*(1 - x) / x )
    prob - x*mix - risk_new*(1 - mix)
  }
  risk_0_first <- uniroot(f6, c(0, 1))$root
  risk_new_first <- (prob - risk_0_first*mix) / (1 - mix)
  
  f7 <- function(x) {
    risk_new <- 1 / ( 1 + exp(max(brange0))*(1 - x) / x )
    prob - x*mix - risk_new*(1 - mix)
  }
  risk_0_second <- uniroot(f7, c(0, 1))$root
  risk_new_second <- (prob - risk_0_second*mix) / (1 - mix)
  
  # Estimate risk_0(0, 0) and risk_0(1, 0) with SACE method
  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    risk_0_first - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_first <- uniroot(f8, c(0, 1))$root
  risk_0_10_first <- (risk_0_first - risk_0_00_first*p_00) / p_10
  
  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    risk_0_first - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_second <- uniroot(f8, c(0, 1))$root
  risk_0_10_second <- (risk_0_first - risk_0_00_second*p_00) / p_10
  
  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    risk_0_second - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_third <- uniroot(f8, c(0, 1))$root
  risk_0_10_third <- (risk_0_second - risk_0_00_third*p_00) / p_10
  
  f8 <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    risk_0_second - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_fourth <- uniroot(f8, c(0, 1))$root
  risk_0_10_fourth <- (risk_0_second - risk_0_00_fourth*p_00) / p_10
  
  # Calculate CEP(1,0) - CEP(0,0) at the 4 boundary points of beta_0, beta_1
  # Take the Ignorance Interval to be the max and min of these estimates
  II1 <- h(risk_1_10, risk_0_10_first) - h(risk_1_00, risk_0_00_first)
  II2 <- h(risk_1_10, risk_0_10_second) - h(risk_1_00, risk_0_00_second)
  II3 <- h(risk_1_10, risk_0_10_third) - h(risk_1_00, risk_0_00_third)
  II4 <- h(risk_1_10, risk_0_10_fourth) - h(risk_1_00, risk_0_00_fourth)
  
  II_low <- min(II1, II2, II3, II4)
  II_up <- max(II1, II2, II3, II4)
  
  # Output estimates
  data.frame(risk_1_00, risk_1_10, mix_1, mix_2, mix, prob, risk_0_first,
             risk_0_second, risk_new_first, risk_new_second, p_10,
             risk_0_00_first, risk_0_00_second, risk_0_00_third,
             risk_0_00_fourth, risk_0_10_first, risk_0_10_second,
             risk_0_10_third, risk_0_10_fourth, II_low, II_up, pi_hat)
  
}