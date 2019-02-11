# Perform analysis from the application section of the manuscript
rm(list=ls())

library(rootSolve)
library(numDeriv)
library(plyr)

source("analyze_NEE_app.R")
source("low_eefun_NEE_app.R")
source("up_eefun_NEE_app.R")

source("analyze_NEB_app.R")
source("low_eefun_NEB_app.R")
source("up_eefun_NEB_app.R")

source("../M Estimation Helper Functions/compute.R")
source("../M Estimation Helper Functions/utilities.R")

# Define helper function for Imbens-Manski interval computation
f_cstar <- function(c,low,up,maxsig){
  pnorm(c + (up-low) / maxsig) - pnorm(-c) - 0.95}

# Save results array (analogous to main sims files)
results_B <- array(NA,dim=c(3,3,6))
results_C <- array(NA,dim=c(3,3,6))

risks_B <- array(NA,dim=c(3,3,6))
risks_C <- array(NA,dim=c(3,3,6))

# Make sure data is in working directory and load
data <- read.csv("HVTN505.csv")

data$R[is.na(data$R)] <- 0
data$R[is.na(data$R)] <- 0
#data$Z[!is.na(data$Z)] <- 1-data$Z[!is.na(data$Z)]
Sstar_num <- c(7, 6, 8)   # Column numbers of data for each S variable
Sstar_list <- c("SstarCD8", "SstarAb", "SstarAbCD8")
# For the plotting of results
markerslabel <- "S1 = CD8 Env Polyfunctionality Score, S2 = IgG Env, S3 = High on Either S1 or S2"

brange_list <- c(0,0.5,1)

# Use a CEP surface equal to VE (original code):
#contrast <- "VE"
# Try it with a log RR contrast:
contrast <- "logRR"

# This code analyzes the 3 markers using the scenario B and C methods in the manuscript
for (i in 1:3) {
  for (j in 1:3) {

    B_ests <- analyze_data_B(data,Sstar_num[i],brange_list[j],contrast)
    data$S_star <- data[,Sstar_num[i]]
    data$S_star[data$R==0] <- 0

    wmin <- B_ests$whichmin; wmax <- B_ests$whichmax
    thetahat_low_B <- as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmin],B_ests[6+wmin],B_ests[9],B_ests[11]))
    thetahat_up_B <-  as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmax],B_ests[6+wmax],B_ests[10:11]))
    if (wmin==1) {low_fun <- ScenarioB_low_eefun_cc_min}
    if (wmin==2) {low_fun <- ScenarioB_low_eefun_cc_max}
    if (wmax==1) {up_fun <- ScenarioB_up_eefun_cc_min}
    if (wmax==2) {up_fun <- ScenarioB_up_eefun_cc_max}
    risks_B[i,j,1] <- B_ests[[14]]; risks_B[i,j,2] <- B_ests[[4]]; risks_B[i,j,3] <- B_ests[[7]]
    risks_B[i,j,4] <- B_ests[[8]]; risks_B[i,j,5] <- B_ests[[5]]; risks_B[i,j,6] <- B_ests[[6]]
    
    data_wide <- count(data,c('Y','Z','Ytau','S_star','R'))
    data_wide$Group <- 1:nrow(data_wide)

    mats <- compute_matrices(list(eeFUN = low_fun,
                          splitdt = split(data_wide, f = data_wide$Group)),
                          numDeriv_options = list(method = 'simple'),
                          theta = thetahat_low_B, beta0range = c(-brange_list[j],brange_list[j]), contrast = "logRR")

    A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
    B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)

    Sigma <- solve(A) %*% B %*% t(solve(A))

    var_Ignorance_low_B_M <- Sigma[7,7]

    if (j==0){
	  	var_Ignorance_up_B_M <-  var_Ignorance_low_B_M}
	  else{
		mats <- compute_matrices(list(eeFUN = up_fun,
                          splitdt = split(data_wide, f = data_wide$Group)),
                          numDeriv_options = list(method = 'simple'),
                          theta = thetahat_up_B, beta0range = c(-brange_list[j],brange_list[j]), contrast = "logRR")
	  	A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
	  	B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
	  	Sigma <- solve(A) %*% B %*% t(solve(A))
  	  	var_Ignorance_up_B_M <- Sigma[7,7]
	  }
   
    results_B[i,j,1] <- as.numeric(B_ests$II_low)
    results_B[i,j,2] <- var_Ignorance_low_B_M
    results_B[i,j,3] <- as.numeric(B_ests$II_up)
    results_B[i,j,4] <- var_Ignorance_up_B_M

    # compute Imbens-Manski interval
    maxsig <- max(sqrt(results_B[i,j,2]),sqrt(results_B[i,j,4]))
    cB <- uniroot(f_cstar, c(1.64, 1.96),
    low=results_B[i,j,1],up=results_B[i,j,3],maxsig=maxsig)$root       
    # Wald symmetric confidence intervals on the scale of the contrast
    results_B[i,j,5] <- results_B[i,j,1] - cB*sqrt(results_B[i,j,2])
    results_B[i,j,6] <- results_B[i,j,3] + cB*sqrt(results_B[i,j,4]) 
    

    C_ests <- analyze_data_C(data,Sstar_num[i],brange_list[j],contrast)
    
    wmin <- C_ests$whichmin; wmax <- C_ests$whichmax
    thetahat_low_C <- as.numeric(c(C_ests[1:6],C_ests[6+round(wmin/1.99,0)],C_ests[8+round(wmin/1.99,0)],C_ests[11],C_ests[11+wmin],C_ests[15+wmin],C_ests[20],C_ests[22]))
    thetahat_up_C <- as.numeric(c(C_ests[1:6],C_ests[6+round(wmax/1.99,0)],C_ests[8+round(wmax/1.99,0)],C_ests[11],C_ests[11+wmax],C_ests[15+wmax],C_ests[21:22]))
    if (wmin==1) {low_fun <- ScenarioC_low_eefun_cc_minmin}
    if (wmin==2) {low_fun <- ScenarioC_low_eefun_cc_minmax}
    if (wmin==3) {low_fun <- ScenarioC_low_eefun_cc_maxmax}
    if (wmin==4) {low_fun <- ScenarioC_low_eefun_cc_maxmin}
    if (wmax==1) {up_fun <- ScenarioC_up_eefun_cc_minmin}
    if (wmax==2) {up_fun <- ScenarioC_up_eefun_cc_minmax}
    if (wmax==3) {up_fun <- ScenarioC_up_eefun_cc_maxmax}
    if (wmax==4) {up_fun <- ScenarioC_up_eefun_cc_maxmin}
    risks_C[i,j,] <- c(C_ests[[2]],C_ests[[1]],C_ests[[15+wmin]],C_ests[[15+wmax]],C_ests[[11+wmin]],C_ests[[11+wmax]])

    mats <- compute_matrices(list(eeFUN = low_fun,
                                  splitdt = split(data_wide, f = data_wide$Group)),
                             numDeriv_options = list(method = 'simple'),
                             theta = thetahat_low_C, beta0range = c(-brange_list[j],brange_list[j]), beta1range = c(-brange_list[j],brange_list[j]), contrast="logRR")
    
    A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
    B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
    
    Sigma <- solve(A) %*% B %*% t(solve(A))
    
    var_Ignorance_low_C_M <- Sigma[12,12]
    
    if (j==0){
      var_Ignorance_up_C_M <-  var_Ignorance_low_C_M}
    else{
      mats <- compute_matrices(list(eeFUN = up_fun,
                                    splitdt = split(data_wide, f = data_wide$Group)),
                               numDeriv_options = list(method = 'simple'),
                               theta = thetahat_up_C, beta0range = c(-brange_list[j],brange_list[j]), beta1range = c(-brange_list[j],brange_list[j]), contrast="logRR")
      A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
      B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
      Sigma <- solve(A) %*% B %*% t(solve(A))
      var_Ignorance_up_C_M <- Sigma[12,12]
    } # else
    
    
    results_C[i,j,1] <- as.numeric(C_ests$II_low)
    results_C[i,j,2] <- var_Ignorance_low_C_M
    results_C[i,j,3] <- as.numeric(C_ests$II_up)
    results_C[i,j,4] <- var_Ignorance_up_C_M
    
    # compute Imbens-Manski interval
    maxsig <- max(sqrt(results_C[i,j,2]),sqrt(results_C[i,j,4]))
    cC <- uniroot(f_cstar, c(1.64, 1.96),
                  low=results_C[i,j,1],up=results_C[i,j,3],maxsig=maxsig)$root     
    results_C[i,j,5] <- results_C[i,j,1] - cC*sqrt(results_C[i,j,2])
    results_C[i,j,6] <- results_C[i,j,3] + cC*sqrt(results_C[i,j,4]) 

  } 
}

source("analyze_NEE_VE1_app.R")
source("low_eefun_NEE_VE1_app.R")
source("up_eefun_NEE_VE1_app.R")

results_B_VE1 <- array(NA,dim=c(3,3,6))
results_B_VE0 <- array(NA,dim=c(3,3,6))

for (i in 1:3) {
  for (j in 1:3) {
    
    B_ests <- analyze_data_B(data,Sstar_num[i],brange_list[j],contrast)
    
    data$S_star <- data[,Sstar_num[i]]
    data$S_star[data$R==0] <- 0
    
    wmin <- B_ests$whichmin; wmax <- B_ests$whichmax
    thetahat_low_B <- as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmin],B_ests[6+wmin],B_ests[9],B_ests[11]))
    thetahat_up_B <-  as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmax],B_ests[6+wmax],B_ests[10:11]))
    if (wmin==1) {low_fun <- ScenarioB_low_eefun_cc_min}
    if (wmin==2) {low_fun <- ScenarioB_low_eefun_cc_max}
    if (wmax==1) {up_fun <- ScenarioB_up_eefun_cc_min}
    if (wmax==2) {up_fun <- ScenarioB_up_eefun_cc_max}
    risks_B[i,j,1] <- B_ests[[14]]; risks_B[i,j,2] <- B_ests[[4]]; risks_B[i,j,3] <- B_ests[[7]]
    risks_B[i,j,4] <- B_ests[[8]]; risks_B[i,j,5] <- B_ests[[5]]; risks_B[i,j,6] <- B_ests[[6]]
    
    data_wide <- count(data,c('Y','Z','Ytau','S_star','R'))
    data_wide$Group <- 1:nrow(data_wide)
    
    mats <- compute_matrices(list(eeFUN = low_fun,
                                  splitdt = split(data_wide, f = data_wide$Group)),
                             numDeriv_options = list(method = 'simple'),
                             theta = thetahat_low_B, beta0range = c(-brange_list[j],brange_list[j]), contrast="logRR")
    
    A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
    B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
    
    Sigma <- solve(A) %*% B %*% t(solve(A))
    
    var_Ignorance_low_B_M <- Sigma[7,7]
    
    if (j==0){
      var_Ignorance_up_B_M <-  var_Ignorance_low_B_M}
    else{
      mats <- compute_matrices(list(eeFUN = up_fun,
                                    splitdt = split(data_wide, f = data_wide$Group)),
                               numDeriv_options = list(method = 'simple'),
                               theta = thetahat_up_B, beta0range = c(-brange_list[j],brange_list[j]), contrast="logRR")
      A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
      B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
      Sigma <- solve(A) %*% B %*% t(solve(A))
      var_Ignorance_up_B_M <- Sigma[7,7]
    }
    
    results_B_VE1[i,j,1] <- as.numeric(B_ests$II_low) # affected by contrast
    results_B_VE1[i,j,2] <- var_Ignorance_low_B_M
    results_B_VE1[i,j,3] <- as.numeric(B_ests$II_up)  # affected by contrast
    results_B_VE1[i,j,4] <- var_Ignorance_up_B_M
    
    # compute Imbens-Manski interval
    maxsig <- max(sqrt(results_B_VE1[i,j,2]),sqrt(results_B_VE1[i,j,4]))
    cB <- uniroot(f_cstar, c(1.64, 1.96),
                  low=results_B_VE1[i,j,1],up=results_B_VE1[i,j,3],maxsig=maxsig)$root       
    results_B_VE1[i,j,5] <- results_B_VE1[i,j,1] - cB*sqrt(results_B_VE1[i,j,2])
    results_B_VE1[i,j,6] <- results_B_VE1[i,j,3] + cB*sqrt(results_B_VE1[i,j,4]) 
    
  } 
}

source("analyze_NEE_VE0_app.R")
source("low_eefun_NEE_VE0_app.R")
source("up_eefun_NEE_VE0_app.R")

for (i in 1:3) {
  for (j in 1:3) {
    
    B_ests <- analyze_data_B(data,Sstar_num[i],brange_list[j],contrast)
    data$S_star <- data[,Sstar_num[i]]
    data$S_star[data$R==0] <- 0
    
    wmin <- B_ests$whichmin; wmax <- B_ests$whichmax
    thetahat_low_B <- as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmin],B_ests[6+wmin],B_ests[9],B_ests[11]))
    thetahat_up_B <-  as.numeric(c(B_ests[1],B_ests[14],B_ests[3:4],B_ests[4+wmax],B_ests[6+wmax],B_ests[10:11]))
    if (wmin==1) {low_fun <- ScenarioB_low_eefun_cc_min}
    if (wmin==2) {low_fun <- ScenarioB_low_eefun_cc_max}
    if (wmax==1) {up_fun <- ScenarioB_up_eefun_cc_min}
    if (wmax==2) {up_fun <- ScenarioB_up_eefun_cc_max}
    risks_B[i,j,1] <- B_ests[[14]]; risks_B[i,j,2] <- B_ests[[4]]; risks_B[i,j,3] <- B_ests[[7]]
    risks_B[i,j,4] <- B_ests[[8]]; risks_B[i,j,5] <- B_ests[[5]]; risks_B[i,j,6] <- B_ests[[6]]
    
    data_wide <- count(data,c('Y','Z','Ytau','S_star','R'))
    data_wide$Group <- 1:nrow(data_wide)
    
    mats <- compute_matrices(list(eeFUN = low_fun,
                                  splitdt = split(data_wide, f = data_wide$Group)),
                             numDeriv_options = list(method = 'simple'),
                             theta = thetahat_low_B, beta0range = c(-brange_list[j],brange_list[j]), contrast="logRR")
    
    A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
    B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
    
    Sigma <- solve(A) %*% B %*% t(solve(A))
    
    var_Ignorance_low_B_M <- Sigma[7,7]
    
    if (j==0){
      var_Ignorance_up_B_M <-  var_Ignorance_low_B_M}
    else{
      mats <- compute_matrices(list(eeFUN = up_fun,
                                    splitdt = split(data_wide, f = data_wide$Group)),
                               numDeriv_options = list(method = 'simple'),
                               theta = thetahat_up_B, beta0range = c(-brange_list[j],brange_list[j]), contrast="logRR")
      A <- apply(simplify2array(Map(`*`, mats$A_i, data_wide$freq)), 1:2, sum)
      B <- apply(simplify2array(Map(`*`, mats$B_i, data_wide$freq)), 1:2, sum)
      Sigma <- solve(A) %*% B %*% t(solve(A))
      var_Ignorance_up_B_M <- Sigma[7,7]
    }
    
    results_B_VE0[i,j,1] <- as.numeric(B_ests$II_low)
    results_B_VE0[i,j,2] <- var_Ignorance_low_B_M
    results_B_VE0[i,j,3] <- as.numeric(B_ests$II_up)
    results_B_VE0[i,j,4] <- var_Ignorance_up_B_M
    
    # compute Imbens-Manski interval
    maxsig <- max(sqrt(results_B_VE0[i,j,2]),sqrt(results_B_VE0[i,j,4]))
    cB <- uniroot(f_cstar, c(1.64, 1.96),
                  low=results_B_VE0[i,j,1],up=results_B_VE0[i,j,3],maxsig=maxsig)$root       
    results_B_VE0[i,j,5] <- results_B_VE0[i,j,1] - cB*sqrt(results_B_VE0[i,j,2])
    results_B_VE0[i,j,6] <- results_B_VE0[i,j,3] + cB*sqrt(results_B_VE0[i,j,4]) 

  } 
}

# Produce figure
if (contrast !="logRR") {
pdf("plot505AbCD8VEscale.pdf") }
if (contrast=="logRR") {
pdf("plot505AbCD8logRRtransformedbacktoVEscale.pdf") }
par(mfrow=c(2,3),cex.axis=1.4,cex.lab=1.2,cex.main=1.4,cex.legend=1.4,oma=c(0,1,0,0),las=1)
labels <- c(expression("S"[1]),expression("S"[2]),expression("S"[3]))
ylimstop <- c(-1,4)
if (contrast=="logRR") {
ylimstop <- c(0,2.1) }
ylimsbottom <- c(-1,1)

# CEP(1) - CEP(0) values:
# Contrast for VE scale:
h <- function(x) { return(x) }
# If contrast=="logRR", transform back to the VE scale: VE(1)-VE(0)
# changes to (1-VE(1))/(1-VE(0)) = RR(1)/RR(0)
contrastlabel <- "VE(1) - VE(0)"
if (contrast=="logRR") {
h <- function(x) {return(exp(x))} 
contrastlabel <- "RR(1)/RR(0)"}
# Reports RR(1)/RR(0)

a11 <- h(results_B[1,1,1]); b11 <- h(results_B[1,1,3]); c11 <- h(results_B[1,1,5]); d11 <- h(results_B[1,1,6])
a12 <- h(results_B[1,2,1]); b12 <- h(results_B[1,2,3]); c12 <- h(results_B[1,2,5]); d12 <- h(results_B[1,2,6])
a13 <- h(results_B[1,3,1]); b13 <- h(results_B[1,3,3]); c13 <- h(results_B[1,3,5]); d13 <- h(results_B[1,3,6])
a21 <- h(results_B[2,1,1]); b21 <- h(results_B[2,1,3]); c21 <- h(results_B[2,1,5]); d21 <- h(results_B[2,1,6])
a22 <- h(results_B[2,2,1]); b22 <- h(results_B[2,2,3]); c22 <- h(results_B[2,2,5]); d22 <- h(results_B[2,2,6])
a23 <- h(results_B[2,3,1]); b23 <- h(results_B[2,3,3]); c23 <- h(results_B[2,3,5]); d23 <- h(results_B[2,3,6])
a31 <- h(results_B[3,1,1]); b31 <- h(results_B[3,1,3]); c31 <- h(results_B[3,1,5]); d31 <- h(results_B[3,1,6])
a32 <- h(results_B[3,2,1]); b32 <- h(results_B[3,2,3]); c32 <- h(results_B[3,2,5]); d32 <- h(results_B[3,2,6])
a33 <- h(results_B[3,3,1]); b33 <- h(results_B[3,3,3]); c33 <- h(results_B[3,3,5]); d33 <- h(results_B[3,3,6])

A11 <- h(results_C[1,1,1]); B11 <- h(results_C[1,1,3]); C11 <- h(results_C[1,1,5]); D11 <- h(results_C[1,1,6])
A12 <- h(results_C[1,2,1]); B12 <- h(results_C[1,2,3]); C12 <- h(results_C[1,2,5]); D12 <- h(results_C[1,2,6])
A13 <- h(results_C[1,3,1]); B13 <- h(results_C[1,3,3]); C13 <- h(results_C[1,3,5]); D13 <- h(results_C[1,3,6])
A21 <- h(results_C[2,1,1]); B21 <- h(results_C[2,1,3]); C21 <- h(results_C[2,1,5]); D21 <- h(results_C[2,1,6])
A22 <- h(results_C[2,2,1]); B22 <- h(results_C[2,2,3]); C22 <- h(results_C[2,2,5]); D22 <- h(results_C[2,2,6])
A23 <- h(results_C[2,3,1]); B23 <- h(results_C[2,3,3]); C23 <- h(results_C[2,3,5]); D23 <- h(results_C[2,3,6])
A31 <- h(results_C[3,1,1]); B31 <- h(results_C[3,1,3]); C31 <- h(results_C[3,1,5]); D31 <- h(results_C[3,1,6])
A32 <- h(results_C[3,2,1]); B32 <- h(results_C[3,2,3]); C32 <- h(results_C[3,2,5]); D32 <- h(results_C[3,2,6])
A33 <- h(results_C[3,3,1]); B33 <- h(results_C[3,3,3]); C33 <- h(results_C[3,3,5]); D33 <- h(results_C[3,3,6])


#pdf("plot1.pdf")
plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimstop,xlim=c(0,6),xaxt="n",xlab="",ylab=contrastlabel,main=expression(paste("[",l[0],",",u[0],"]=[0,0]")))
arrows(0.6,a11,0.6,b11,length=0.05,angle=90,code=3)
arrows(0.6,c11,0.6,d11,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a21,2.6,b21,length=0.05,angle=90,code=3)
arrows(2.6,c21,2.6,d21,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a31,4.6,b31,length=0.05,angle=90,code=3)
arrows(4.6,c31,4.6,d31,length=0.05,angle=90,code=3,lty=2)

points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A11,1.4,B11,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C11,1.4,D11,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A21,3.4,B21,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C21,3.4,D21,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A31,5.4,B31,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C31,5.4,D31,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

#legend("topleft",legend=c("Scen. B","Scen. C"),col=c("black","dimgrey"),lty=1,cex=1,bty='n',y.intersp = 0.35,inset=-0.03)
legend("topleft",legend=c("Scen. B","Scen. C"),col=c("black","dimgrey"),lty=1,cex=1,bty='n')
axis(1, c(1,3,5), labels)

plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimstop,xlim=c(0,6),xaxt="n",xlab="",ylab=contrastlabel,main=expression(paste("[",l[0],",",u[0],"]=[-0.5,0.5]")))
arrows(0.6,a12,0.6,b12,length=0.05,angle=90,code=3)
arrows(0.6,c12,0.6,d12,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a22,2.6,b22,length=0.05,angle=90,code=3)
arrows(2.6,c22,2.6,d22,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a32,4.6,b32,length=0.05,angle=90,code=3)
arrows(4.6,c32,4.6,d32,length=0.05,angle=90,code=3,lty=2)

points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A12,1.4,B12,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C12,1.4,D12,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A22,3.4,B22,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C22,3.4,D22,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A32,5.4,B32,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C32,5.4,D32,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

legend("topleft",legend=c("Scen. B","Scen. C"),col=c("black","dimgrey"),lty=1,cex=1,bty='n')
axis(1, c(1,3,5), labels)

plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimstop,xlim=c(0,6),xaxt="n",xlab="",ylab=contrastlabel,main=expression(paste("[",l[0],",",u[0],"]=[-1,1]")))
arrows(0.6,a13,0.6,b13,length=0.05,angle=90,code=3)
arrows(0.6,c13,0.6,d13,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a23,2.6,b23,length=0.05,angle=90,code=3)
arrows(2.6,c23,2.6,d23,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a33,4.6,b33,length=0.05,angle=90,code=3)
arrows(4.6,c33,4.6,d33,length=0.05,angle=90,code=3,lty=2)

points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A13,1.4,B13,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C13,1.4,D13,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A23,3.4,B23,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C23,3.4,D23,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A33,5.4,B33,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C33,5.4,D33,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

legend("topleft",legend=c("Scen. B","Scen. C"),col=c("black","dimgrey"),lty=1,cex=1,bty='n')
axis(1, c(1,3,5), labels)


# If contrast=="logRR", transform back to the VE scale:
h <- function(x) { return(x) }
if (contrast=="logRR") {
h <- function(x) {return(1-exp(x))} }

a11 <- h(results_B_VE1[1,1,1]); b11 <- h(results_B_VE1[1,1,3]); c11 <- h(results_B_VE1[1,1,5]); d11 <- h(results_B_VE1[1,1,6])
a12 <- h(results_B_VE1[1,2,1]); b12 <- h(results_B_VE1[1,2,3]); c12 <- h(results_B_VE1[1,2,5]); d12 <- h(results_B_VE1[1,2,6])
a13 <- h(results_B_VE1[1,3,1]); b13 <- h(results_B_VE1[1,3,3]); c13 <- h(results_B_VE1[1,3,5]); d13 <- h(results_B_VE1[1,3,6])
a21 <- h(results_B_VE1[2,1,1]); b21 <- h(results_B_VE1[2,1,3]); c21 <- h(results_B_VE1[2,1,5]); d21 <- h(results_B_VE1[2,1,6])
a22 <- h(results_B_VE1[2,2,1]); b22 <- h(results_B_VE1[2,2,3]); c22 <- h(results_B_VE1[2,2,5]); d22 <- h(results_B_VE1[2,2,6])
a23 <- h(results_B_VE1[2,3,1]); b23 <- h(results_B_VE1[2,3,3]); c23 <- h(results_B_VE1[2,3,5]); d23 <- h(results_B_VE1[2,3,6])
a31 <- h(results_B_VE1[3,1,1]); b31 <- h(results_B_VE1[3,1,3]); c31 <- h(results_B_VE1[3,1,5]); d31 <- h(results_B_VE1[3,1,6])
a32 <- h(results_B_VE1[3,2,1]); b32 <- h(results_B_VE1[3,2,3]); c32 <- h(results_B_VE1[3,2,5]); d32 <- h(results_B_VE1[3,2,6])
a33 <- h(results_B_VE1[3,3,1]); b33 <- h(results_B_VE1[3,3,3]); c33 <- h(results_B_VE1[3,3,5]); d33 <- h(results_B_VE1[3,3,6])

A11 <- h(results_B_VE0[1,1,1]); B11 <- h(results_B_VE0[1,1,3]); C11 <- h(results_B_VE0[1,1,5]); D11 <- h(results_B_VE0[1,1,6])
A12 <- h(results_B_VE0[1,2,1]); B12 <- h(results_B_VE0[1,2,3]); C12 <- h(results_B_VE0[1,2,5]); D12 <- h(results_B_VE0[1,2,6])
A13 <- h(results_B_VE0[1,3,1]); B13 <- h(results_B_VE0[1,3,3]); C13 <- h(results_B_VE0[1,3,5]); D13 <- h(results_B_VE0[1,3,6])
A21 <- h(results_B_VE0[2,1,1]); B21 <- h(results_B_VE0[2,1,3]); C21 <- h(results_B_VE0[2,1,5]); D21 <- h(results_B_VE0[2,1,6])
A22 <- h(results_B_VE0[2,2,1]); B22 <- h(results_B_VE0[2,2,3]); C22 <- h(results_B_VE0[2,2,5]); D22 <- h(results_B_VE0[2,2,6])
A23 <- h(results_B_VE0[2,3,1]); B23 <- h(results_B_VE0[2,3,3]); C23 <- h(results_B_VE0[2,3,5]); D23 <- h(results_B_VE0[2,3,6])
A31 <- h(results_B_VE0[3,1,1]); B31 <- h(results_B_VE0[3,1,3]); C31 <- h(results_B_VE0[3,1,5]); D31 <- h(results_B_VE0[3,1,6])
A32 <- h(results_B_VE0[3,2,1]); B32 <- h(results_B_VE0[3,2,3]); C32 <- h(results_B_VE0[3,2,5]); D32 <- h(results_B_VE0[3,2,6])
A33 <- h(results_B_VE0[3,3,1]); B33 <- h(results_B_VE0[3,3,3]); C33 <- h(results_B_VE0[3,3,5]); D33 <- h(results_B_VE0[3,3,6])

#pdf("plot4.pdf")
plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimsbottom,xlim=c(0,6),xaxt="n",xlab="",ylab="VE(s)",main=expression(paste("[",l[0],",",u[0],"]=[0,0]")))
arrows(0.6,a11,0.6,b11,length=0.05,angle=90,code=3)
arrows(0.6,c11,0.6,d11,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a21,2.6,b21,length=0.05,angle=90,code=3)
arrows(2.6,c21,2.6,d21,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a31,4.6,b31,length=0.05,angle=90,code=3)
arrows(4.6,c31,4.6,d31,length=0.05,angle=90,code=3,lty=2)


points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A11,1.4,B11,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C11,1.4,D11,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A21,3.4,B21,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C21,3.4,D21,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A31,5.4,B31,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C31,5.4,D31,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

legend("bottomleft",legend=c("VE(1)","VE(0)"),col=c("black","dimgrey"),lty=1,cex=1)
axis(1, c(1,3,5), labels)

plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimsbottom,xlim=c(0,6),xaxt="n",xlab="",ylab="VE(s)",main=expression(paste("[",l[0],",",u[0],"]=[-0.5,0.5]")))
arrows(0.6,a12,0.6,b12,length=0.05,angle=90,code=3)
arrows(0.6,c12,0.6,d12,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a22,2.6,b22,length=0.05,angle=90,code=3)
arrows(2.6,c22,2.6,d22,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a32,4.6,b32,length=0.05,angle=90,code=3)
arrows(4.6,c32,4.6,d32,length=0.05,angle=90,code=3,lty=2)

points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A12,1.4,B12,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C12,1.4,D12,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A22,3.4,B22,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C22,3.4,D22,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A32,5.4,B32,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C32,5.4,D32,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

legend("bottomleft",legend=c("VE(1)","VE(0)"),col=c("black","dimgrey"),lty=1,cex=1)
axis(1, c(1,3,5), labels)

plot(c(.6,2.6,4.6),c((a11+b11)/2,(a21+b21)/2,(a31+b31)/2),pch=19,ylim=ylimsbottom,xlim=c(0,6),xaxt="n",xlab="",ylab="VE(s)",main=expression(paste("[",l[0],",",u[0],"]=[-1,1]")))
arrows(0.6,a13,0.6,b13,length=0.05,angle=90,code=3)
arrows(0.6,c13,0.6,d13,length=0.05,angle=90,code=3,lty=2)
arrows(2.6,a23,2.6,b23,length=0.05,angle=90,code=3)
arrows(2.6,c23,2.6,d23,length=0.05,angle=90,code=3,lty=2)
arrows(4.6,a33,4.6,b33,length=0.05,angle=90,code=3)
arrows(4.6,c33,4.6,d33,length=0.05,angle=90,code=3,lty=2)

points(c(1.4,3.4,5.4),c((A11+B11)/2,(A21+B21)/2,(A31+B31)/2),pch=19,col="dimgrey")
arrows(1.4,A13,1.4,B13,length=0.05,angle=90,code=3,col="dimgrey")
arrows(1.4,C13,1.4,D13,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(3.4,A23,3.4,B23,length=0.05,angle=90,code=3,col="dimgrey")
arrows(3.4,C23,3.4,D23,length=0.05,angle=90,code=3,lty=2,col="dimgrey")
arrows(5.4,A33,5.4,B33,length=0.05,angle=90,code=3,col="dimgrey")
arrows(5.4,C33,5.4,D33,length=0.05,angle=90,code=3,lty=2,col="dimgrey")

legend("bottomleft",legend=c("VE(1)","VE(0)"),col=c("black","dimgrey"),lty=1,cex=1)
axis(1, c(1,3,5), labels)

# mtext(markerslabel,outer=T,side=1,line=0)
dev.off()

# Values of interest for supplemental simulation

1-mean(data$Ytau[data$Z == 1])+mean(data$Ytau[data$Z == 0])  # P(Ytau(1)=1,Ytau(0)=0)

