# Make plots in NEH Simulations in manuscript
rm(list=ls())

# Set working directory to where results are saved
# Change names of files loaded to match number of sims run in main_NEH.R

# Load simulation results
load("results2000sims.0.1cc.Rdata"); results.0.1cc <- results; rm(results)
load("results2000sims.0.25cc.Rdata"); results.0.25cc <- results; rm(results)
load("results2000sims.1cc.Rdata")

# results array
# 1 \hat CEP_l, 2 \hat sigma^2_l, 3 \hat CEP_u, 4 \hat sigma^2_u
# 5 EUI lower limit, 6 EUI upper limit

# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 0, S^tau(0) = 0) = a
varyingprobs1 <- seq(0.7, 0.1, by = -0.05)
# P(Y = 1 | Y^tau(1) = Y^tau(0) = 0, S^tau(1) = 1, S^tau(0) = 0) = b
varyingprobs2 <- seq(0.1, 0.7, by = 0.05)

# Function for panel titles in figures
maketitle <- function(h, i) {
	
  if (h==1 & i==1) title(expression( paste("n=",2000,", [",l[0],",",u[0],"]=[0,0]")))
  if (h==1 & i==2) title(expression( paste("n=",4000,", [",l[0],",",u[0],"]=[0,0]")))
  if (h==1 & i==3) title(expression( paste("n=",8000,", [",l[0],",",u[0],"]=[0,0]")))

  if (h==2 & i==1) title(expression( paste("n=",2000,", [",l[0],",",u[0],"]=[-0.5,0.5]")))
  if (h==2 & i==2) title(expression( paste("n=",4000,", [",l[0],",",u[0],"]=[-0.5,0.5]")))
  if (h==2 & i==3) title(expression( paste("n=",8000,", [",l[0],",",u[0],"]=[-0.5,0.5]")))

  if (h==3 & i==1) title(expression( paste("n=",2000,", [",l[0],",",u[0],"]=[-1,1]")))
  if (h==3 & i==2) title(expression( paste("n=",4000,", [",l[0],",",u[0],"]=[-1,1]")))
  if (h==3 & i==3) title(expression( paste("n=",8000,", [",l[0],",",u[0],"]=[-1,1]")))

}


############################################################################################################3

# Power plot
pdf("Power_NEH.pdf")
par(mfrow=c(3, 3))

for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- c(-1, 1)
  
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Initialize vectors containing x and y coordinates for plot points
    yy <- yy.0.25cc <- yy.0.1cc <- xx <- matrix(NA, 1, length(varyingprobs1))
  
    for (j in 1:length(varyingprobs1)) {	
      
      # Power results (y-axis)
    	yy[j] <- sum(results[h, i, j, , 5] > 0 | results[h, i, j, , 6] < 0) /
    	         length(results[h, i, j, , 5])
    	yy.0.25cc[j] <- sum(results.0.25cc[h, i, j, , 5] > 0 | results.0.25cc[h, i, j, , 6] < 0) /
    	                length(results.0.25cc[h, i, j, , 5])
    	yy.0.1cc[j] <- sum(results.0.1cc[h, i, j, , 5] > 0 | results.0.1cc[h, i, j, , 6] < 0) /
    	               length(results.0.1cc[h, i, j, , 5])
  
    	# Varying terms that determine x-axis
    	Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      
      # x-axis
      xx[j] <- TrueCEP_10 - TrueCEP_00      

	  } # end of 'j' loop
    
  plot(xx, yy, ylim=c(0, 1), ylab="Power", xlab="CEP(1,0)-CEP(0,0)", type="l")
  abline(0.05, 0, col="gray")
  lines(xx, yy.0.1cc, lty=2)
  lines(xx, yy.0.25cc, lty=3)
  maketitle(h, i)

  } # end of 'i' loop
  
} # end of 'h' loop
   
dev.off()

###########################################################################################################

# EUI width plot
pdf("EUI_NEH.pdf")
par(mfrow=c(3, 3))

for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- c(-1, 1)
  
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Initialize vectors containing x and y coordinates for plot points
    yy <- yy.0.25cc <- yy.0.1cc <- xx <- matrix(NA, 1, length(varyingprobs1))
  
    for (j in 1:length(varyingprobs1)) {	
      
      # EUI widths (y-axis)
    	yy[j] <- mean(results[h, i, j, , 6] - results[h, i, j, , 5])
    	yy.0.25cc[j] <- mean(results.0.25cc[h, i, j, , 6] - results.0.25cc[h, i, j, , 5])
    	yy.0.1cc[j] <- mean(results.0.1cc[h,i,j,,6] - results.0.1cc[h,i,j,,5])
    
    	# Varying terms that determine x-axis
    	Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      
      # x-axis
      xx[j] <- TrueCEP_10 - TrueCEP_00     
      
    } # end of 'j' loop
    
  plot(xx, yy, ylim=c(0, 2.2), ylab="EUI Width", xlab="CEP(1,0)-CEP(0,0)",
       type="l")
  lines(xx, yy.0.1cc, lty=2)
  lines(xx, yy.0.25cc, lty=3)
  maketitle(h, i)

  } # end of 'i' loop
  
} # end of 'h' loop
   
dev.off()


###########################################################################################################

# Coverage plot
pdf("Coverage_NEH.pdf")
par(mfrow=c(3, 3))

for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- c(-1, 1)

  
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Initialize vectors containing x and y coordinates for plot points
    yy <- yy.0.25cc <- yy.0.1cc <- xx <- matrix(NA, 1, length(varyingprobs1))

    for (j in 1:length(varyingprobs1)) {		

      # Varying terms that determine x-axis
  	  Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      
      # x-axis
      xx[j] <- TrueDiff <- TrueCEP_10 - TrueCEP_00   	

      # EUI coverage (y-axis)
    	yy[j] <- sum(results[h, i, j, , 5] <= TrueDiff & TrueDiff <= results[h, i, j, , 6]) /
    	         length(results[h, i, j, , 5])
    	yy.0.25cc[j]  <- sum(results.0.25cc[h, i, j, , 5] <= TrueDiff & TrueDiff <= results.0.25cc[h, i, j, , 6]) /
    	                 length(results.0.25cc[h, i, j, ,5])
    	yy.0.1cc[j]  <- sum(results.0.1cc[h, i, j, , 5] <= TrueDiff & TrueDiff <= results.0.1cc[h, i, j, , 6]) /
    	                length(results.0.1cc[h, i, j, , 5])
    
    } # end of 'j' loop
    
    plot(xx, yy, ylim=c(0.8, 1), ylab="EUI Coverage",
         xlab="CEP(1,0)-CEP(0,0)", type="l", yaxt="n")
    axis(2, at=c(0.8, 0.95, 1))
    lines(xx, yy.0.1cc, lty=2)
    lines(xx, yy.0.25cc, lty=3)
    maketitle(h, i)

  } # end of 'i' loop

} # end of 'h' loop
   
dev.off()

###########################################################################################################

# Bias plot
pdf("Bias_NEH.pdf")
par(mfrow=c(3, 3))

for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- c(-1, 1)
  
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Initialize vectors containing x and y coordinates for plot points
    yy.l <- yy.u <- yy.0.25cc.l <- yy.0.25cc.u <- yy.0.1cc.l <- yy.0.1cc.u <- xx <- matrix(NA, 1, length(varyingprobs1))

    for (j in 1:length(varyingprobs1)) {		

      # Varying terms that determine x-axis
  	  Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      
      # x-axis
      xx[j] <- TrueDiff <- TrueCEP_10 - TrueCEP_00   	

      # Bias (y-axis)
    	yy.l[j] <- mean(results[h, i, j, , 1]- TrueDiff)
    	yy.u[j] <- mean(results[h, i, j, , 3]- TrueDiff)
    	yy.0.25cc.l[j]  <- mean(results.0.25cc[h, i, j, , 1] - TrueDiff)
    	yy.0.25cc.u[j]  <- mean(results.0.25cc[h, i, j, , 3] - TrueDiff)
    	yy.0.1cc.l[j]  <- mean(results.0.1cc[h, i, j, , 1] - TrueDiff)
    	yy.0.1cc.u[j]  <- mean(results.0.1cc[h, i, j, , 3] - TrueDiff)
     	
    } # end of 'j' loop
    
    plot(xx, yy.l, ylim=c(-1, 1), ylab="Bias", xlab="CEP(1,0)-CEP(0,0)",
         type="l")
    lines(xx, yy.u)
    lines(xx, yy.0.1cc.l, lty=2)
    lines(xx, yy.0.25cc.l, lty=3)
    lines(xx, yy.0.1cc.u, lty=2)
    lines(xx, yy.0.25cc.u, lty=3)
    maketitle(h, i)

  } # end of 'i' loop
  
} # end of 'h' loop
   
dev.off()

###########################################################################################################

# SE plot
pdf("SE_NEH.pdf")
par(mfrow=c(3, 3))

for (h in 1:3) {

  if (h == 1) beta0range <- c(0, 0)
  if (h == 2) beta0range <- c(-0.5, 0.5)
  if (h == 3) beta0range <- c(-1, 1)
  
  for (i in 1:3) {

    if (i == 1) n <- 2000
    if (i == 2) n <- 4000
    if (i == 3) n <- 8000
   
    # Initialize vectors containing x and y coordinates for plot points
    yy.l <- yy.u <- yy.0.25cc.l <- yy.0.25cc.u <- yy.0.1cc.l <- yy.0.1cc.u <- xx <- matrix(NA, 1, length(varyingprobs1))

    for (j in 1:length(varyingprobs1)) {		

      # Varying terms that determine x-axis
  	  Y_1_probvector <- c(varyingprobs1[j], 0, varyingprobs2[j], 0)
      Y_0_probvector <- c(0.5, 0, 0.5, 0)
      TrueCEP_00 <- Y_1_probvector[1] - Y_0_probvector[1]
      TrueCEP_10 <- Y_1_probvector[3] - Y_0_probvector[3]
      
      # x-axis
      xx[j] <- TrueDiff <- TrueCEP_10 - TrueCEP_00   	

      # ESE/ASE (y-axis)
    	yy.l[j] <-  sd(results[h, i, j, , 1]) / mean(sqrt(results[h, i, j, , 2]))
    	yy.u[j] <-  sd(results[h, i, j, , 3]) / mean(sqrt(results[h, i, j, , 4]))
    	yy.0.25cc.l[j]  <- sd(results.0.25cc[h, i, j, , 1]) / mean(sqrt(results.0.25cc[h, i, j, , 2]))
    	yy.0.25cc.u[j]  <- sd(results.0.25cc[h, i, j, , 3]) / mean(sqrt(results.0.25cc[h, i, j, , 4]))
    	yy.0.1cc.l[j]  <- sd(results.0.1cc[h, i, j, , 1]) / mean(sqrt(results.0.1cc[h, i, j, , 2]))
    	yy.0.1cc.u[j]  <- sd(results.0.1cc[h, i, j, , 3]) / mean(sqrt(results.0.1cc[h, i, j, , 4]))
    
    } # end of 'j' loop
    
    plot(xx, yy.l, ylim=c(0, 2), ylab="ESE/ASE", xlab="CEP(1,0)-CEP(0,0)",
         type="l")
    lines(xx, yy.u)
    lines(xx, yy.0.1cc.l, lty=2)
    lines(xx, yy.0.25cc.l, lty=3)
    lines(xx, yy.0.1cc.u, lty=2)
    lines(xx, yy.0.25cc.u, lty=3)
    maketitle(h, i)

  } # end of 'i' loop

} # end of 'h' loop
   
dev.off()