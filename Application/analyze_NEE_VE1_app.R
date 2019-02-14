analyze_data_NEE_VE1 <- function(data, S_star_choice, beta,
                                 contrast="logRR") {

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
  S_star[is.na(S_star)] <- 0
  brange <- c(-beta, beta)
  
  # estimate probability control has S_star observed
  f0 <- function(x){ sum((1 - Y)*(R - x)) }
  pi_hat <- uniroot(f0, c(0, 1))$root			# P[R=1|Y=0]

  W <- 1/pi_hat*(1-Y)*R + Y	# weights

  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star
        
  # Estimate identifiable parameters risk_z and p_00, p_10, p_11
  
  f1 <- function(x) {sum((1-Y_tau)*(1-Z)*(Y-x))}
  risk_0 <- uniroot(f1,c(0,1))$root				# Pr[Y(0)=1|Y^tau(0)=Y^tau(1)=0]
  
  f2 <- function(x) {sum((1-Y_tau)*Z*(Y-x))}
  risk_1 <- uniroot(f2,c(0,1))$root				# Pr[Y(1)=1|Y^tau(0)=Y^tau(1)=0]
  
  f3 <- function(x) {sum((1-Y_tau)*Z*(1-S_star-x)*W)}	# Pr[S_star(1)=S_star(0)=0|Y^tau(0)=Y^tau(1)=0], recall Case CB: Pr[S_star(0)=0|Y^tau(0)=Y^tau(1)=0]=1
  p_00 <- uniroot(f3,c(0,1))$root						
  p_10 <- 1 - p_00
  
  # Estimate identifiable parameters risk_1(0,0) and risk_1(1,0)
  # where risk_z(s1,s0) = Pr[Y(z)=1 | S*(1)=s1, S*(0)=s0, Ytau(1)=Ytau(0)=0 ]

  f5 <- function(x) {sum(Z*S_0*(Y-x)*W) }
  risk_1_00 <- uniroot(f5,c(0,1))$root	
  f6 <- function(x) {sum(Z*S_1*(Y-x)*W) }
  risk_1_10 <- uniroot(f6,c(0,1))$root

  # Estimate risk_0_00 using SACE method, the only partially identifiable term for Scenario B
  
  f7 <- function(x){
	risk_0_10 <- 1/( 1 + exp(min(brange))*(1-x)/x )
	risk_0-x*p_00 -risk_0_10*p_10
	}
  risk_0_00_first <- uniroot(f7,c(0,1))$root
  risk_0_10_first <- (risk_0 - risk_0_00_first*p_00)/p_10

  f8 <- function(x){
	risk_0_10 <- 1/( 1 + exp(max(brange))*(1-x)/x )
	risk_0-x*p_00 -risk_0_10*p_10
	}
  risk_0_00_second <- uniroot(f8,c(0,1))$root
  risk_0_10_second <- (risk_0 - risk_0_00_second*p_00)/p_10
  
  #risk_1_10 <- (risk_1 - p_00*risk_1_00) / p_10
  
  # Calculate CEP(1,0) - CEP(0,0) at the 2 boundary points of beta0 for Scenario B
  # Take the Ignorance Intervals for B to be the max and min of these estimates
    
#  II1 <- (1- risk_1_10/risk_0_10_first)
#  II2 <- (1- risk_1_10/risk_0_10_second)
  II1 <- h(risk_1_10,risk_0_10_first)
  II2 <- h(risk_1_10,risk_0_10_second)
  II_low <- min(II1,II2)


  II_low <- min(II1,II2)
  II_up <- max(II1,II2)

  whichmin <- which.min(c(II1,II2)); print(whichmin)
  whichmax <- which.max(c(II1,II2)); print(whichmax)
  
  data.frame(risk_0,risk_1,  p_00, risk_1_00,  
             risk_0_00_first,risk_0_00_second, 
             risk_0_10_first,risk_0_10_second, 
             II_low,II_up,pi_hat,whichmin,whichmax,risk_1_10)

}

