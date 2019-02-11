analyze_data_C <- function(data,S_star_choice,beta,contrast="logRR") {
  
  h <- function(x,y) {
    ans <- log(x) - log(y)
    if (contrast=="Difference") { ans <- x - y }
    if (contrast=="VE") { ans <- 1-x/y }
    return(ans)
  }
  
  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Ytau
  S_star <- data[,S_star_choice]
  R <- data$R
  S_star[R==0] <- 0
  brange0 <- c(-beta,beta); brange1 <- c(-beta,beta); 
  
  # estimate probability control has S_star observed
  f0 <- function(x){sum( (1-Y)*(R-x)  )}
  pi_hat <- uniroot(f0,c(0,1))$root			# P[R=1|Y=0]
  
  W <- 1/pi_hat*(1-Y)*R + Y	# weights
  
  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star
  
  
  
  
  # Estimate p_10
  
  
  f3 <- function(x) {
    sum((1-Y_tau)*Z*(S_star-x)*W)
  }
  p_10 <- uniroot(f3,c(0,1))$root
  p_00 <- 1 - p_10
  
  
  
  # Estimate risk_1_00, risk_1_10
  
  f4 <- function(x) {
    sum((1-Y_tau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-x)*W)
  }
  risk_1_00 <- uniroot(f4,c(0,1))$root
  
  f5 <- function(x) {
    sum((1-Y_tau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-x)*W)
  }
  risk_1_10 <- uniroot(f5,c(0,1))$root
  
  
  
  # Estimate risk_0
  
  mix_1 <- 1-mean(Y_tau[Z == 1]); mix_2 <- 1-mean(Y_tau[Z == 0])
  mix <- mix_1/mix_2
  prob <- mean(Y[Y_tau == 0 & Z == 0])
  
  f6 <- function(x) {
    risk_new <- 1/( 1 + exp(min(brange0))*(1-x)/x )
    prob - x*mix - risk_new*(1-mix)
  }
  risk_0_first <- uniroot(f6,c(0,1))$root
  risk_new_first <- (prob - risk_0_first*mix)/(1-mix)
  
  f7 <- function(x) {
    risk_new <- 1/( 1 + exp(max(brange0))*(1-x)/x )
    prob - x*mix - risk_new*(1-mix)
  }
  risk_0_second <- uniroot(f7,c(0,1))$root
  risk_new_second <- (prob - risk_0_second*mix)/(1-mix)
  
  
  
  
  # Estimate risk_0_00 and risk_0_10
  
  f8 <- function(x){
    risk_0_10 <- 1/( 1 + exp(min(brange1))*(1-x)/x )
    risk_0_first-x*p_00 -risk_0_10*p_10
  }
  risk_0_00_first <- uniroot(f8,c(0,1))$root
  risk_0_10_first <- (risk_0_first - risk_0_00_first*p_00)/p_10
  
  f8 <- function(x){
    risk_0_10 <- 1/( 1 + exp(max(brange1))*(1-x)/x )
    risk_0_first-x*p_00 -risk_0_10*p_10
  }
  risk_0_00_second <- uniroot(f8,c(0,1))$root
  risk_0_10_second <- (risk_0_first - risk_0_00_second*p_00)/p_10
  
  f8 <- function(x){
    risk_0_10 <- 1/( 1 + exp(max(brange1))*(1-x)/x )
    risk_0_second-x*p_00 -risk_0_10*p_10
  }
  risk_0_00_third <- uniroot(f8,c(0,1))$root
  risk_0_10_third <- (risk_0_second - risk_0_00_third*p_00)/p_10
  
  f8 <- function(x){
    risk_0_10 <- 1/( 1 + exp(min(brange1))*(1-x)/x )
    risk_0_second-x*p_00 -risk_0_10*p_10
  }
  risk_0_00_fourth <- uniroot(f8,c(0,1))$root
  risk_0_10_fourth <- (risk_0_second - risk_0_00_fourth*p_00)/p_10
  
  
  
  
  # Calculate CEP(1,0) - CEP(0,0) at the 4 boundary points of beta0,beta1 for Scenario C
  # Take the Ignorance Intervals for C to be the max and min of these estimates
  
  #II1 <- (risk_1_00/risk_0_00_first) - (risk_1_10/risk_0_10_first)
  #II2 <- (risk_1_00/risk_0_00_second) - (risk_1_10/risk_0_10_second)
  #II3 <- (risk_1_00/risk_0_00_third) - (risk_1_10/risk_0_10_third)
  #II4 <- (risk_1_00/risk_0_00_fourth) - (risk_1_10/risk_0_10_fourth)
  II1 <- h(risk_1_10,risk_0_10_first) - h(risk_1_00,risk_0_00_first)
  II2 <- h(risk_1_10,risk_0_10_second) - h(risk_1_00,risk_0_00_second)
  II3 <- h(risk_1_10,risk_0_10_third) - h(risk_1_00,risk_0_00_third)
  II4 <- h(risk_1_10,risk_0_10_fourth) - h(risk_1_00,risk_0_00_fourth)
  
  II_low <- min(II1,II2,II3,II4)
  II_up <- max(II1,II2,II3,II4)
  
  whichmin <- which.min(c(II1,II2,II3,II4)) 
  whichmax <- which.max(c(II1,II2,II3,II4))
  
  data.frame(risk_1_00,risk_1_10,mix_1,mix_2,mix,prob,risk_0_first,risk_0_second,
             risk_new_first,risk_new_second,p_10,
             risk_0_00_first,risk_0_00_second,risk_0_00_third,risk_0_00_fourth,
             risk_0_10_first,risk_0_10_second,risk_0_10_third,risk_0_10_fourth,
             II_low,II_up,pi_hat,whichmin,whichmax)
  
}