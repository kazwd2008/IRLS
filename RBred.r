##############################################################################################
#  RBred:  Robust estimator for the generalized ratio model   
#          by iteratively re-weighted least squares (IRLS) algorithm for M-estimation  
#   Weight function:  	Tukey's biweight function
#   Scale: 		Average absolute deviation (AAD)
#---------------------------------------------------------------------------------------------#
#       Implemented by K. Wada (NSTAC, Japan)	
#---------------------------------------------------------------------------------------------#
#       Ver. 0.0 [2016.07.22]  
#       Ver. 0.1 [2018.05.21]  A bug to replace the value of scale parameter was fixed. #A#
#       Ver. 0.2 [2018.05.21]  The iteration count improved. 
#       Ver. 0.3 [2019.07.08]  b1 => g1
#       Ver. 0.4 [2019.07.25]  Estimation of g1 with weights corrected, and
#                              loop count inproved.                               
#---------------------------------------------------------------------------------------------#
#  Parameters
#   x1       single explanatory variable 
#   y1       objective variable to be imputed  
#   g1       initial value of power (default: g1=0.5) 
#   c1       tuning constant for the biweight function  
#            the value should be 4 to 8 (default: c=8)ÅB 
#               [ c=4 :   very robust,      c=8 :    less robust]
#   dat      name of the dataframe containing x1 and y1 if any 
#   rp.max   maximum number of iteration (default: rp.max=100) 
#   cg1      convergence condition to stop iteration (default: cg1=0.001)
#---------------------------------------------------------------------------------------------#
#  Return values
#   rt       robust estimation of the rate between x1 and y1 
#   g1       estimated power
#   wt       robust weight 
#   rp       number of iteration for each step
#   fg       convergence flag (fg=0: converged; fg=1: not convergedÅj
#   rt.cg    histry of the rate 
#   g1.cg    history of the power
#   s1.cg    history of the scale parameter
################################################################################################

RBred <- function(x1, y1, g1=0.5, c1=8, dat="", rp.max=100,  cg1=0.001){

  # initial estimation
  if (dat!="") attach(get(dat))	                    
  s1.cg <- rt.cg <- g1.cg  <- matrix(NA, nr=rp.max, nc=2)  # preserve histories
  Trp <- rp1     <- 1  		# count of iteration(1)
  rp2  <- rp3 <- rp4   <- 0  		# count of iteration(2)
  s0      <- 0	    	# initial value of scale parameter

####---------------------------------------------------------------------------------------(0)
   (Trt1 <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1))))	# ratio estimation
   rs1.x <- y1 - Trt1 * x1			        # heteroscedastic residuals 
   rs1   <- y1 / x1^g1- Trt1 * x1^(1-g1)	        # homoscedastic quasi residuals 

   g1.cg[rp1,1]  <- g1
   rt.cg[rp1,1]  <- Trt1
   s1 <- s1.cg[rp1,1] <- mean(abs(rs1))	      # average absolute deviation (AAD) of quasi residuals 
####---------------------------------------------------------------------------------------(1)
#### Estimate rate (Trt1) and power (g1) without weight
  for (i in 1:10) {
       rp1 <- rp1 + 1; Trp <- Trp + 1
       ix1 <- which((x1 !=0) & (rs1 !=0))
       x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
       g1<- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1])))[2]      # estimate power (g1)

       (Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1))))
       rs1.x  <- y1 - Trt1 * x1		                      # heteroscedastic residuals 
       rs1    <- y1 / x1^g1- Trt1 * x1^(1-g1)	              # homoscedastic quasi residuals 

       s0 <- s1
       s1 <- s1.cg[rp1,1] <- mean(abs(rs1))	              # AAD of quasi residuals
       g1.cg[rp1,1]  <- g1
       rt.cg[rp1,1]  <- Trt1
       if (abs(1-s1/s0) < cg1) break 
  }

  ####---------------------------------------------------------------------------------------(2)
  ####  Robust estimation 
  ####---------------------------------------------------------------------------------------
  for (i in 1:10) {   # robust estimation of beta
      rp2 <- rp2 + 1; Trp <- Trp + 1

      # renew weights 
      u1 <- rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0

      ix1 <- which((w1*x1) !=0)     # remove data of w1*x1=0, since they make Trt1 and rs1 NaN 
      if (length(ix1) < 2) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))		
          ix1 <- which((w1*x1) !=0)
      }
      (Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1)))) 
      rs1.x  <- y1 - Trt1 * x1		            # heteroscedastic residuals
      rs1    <- y1 / x1^g1- Trt1 * x1^(1-g1)	    # homoscedastic quasi residuals

      # M-estimation of ratio with fixed power (g1)
      s0 <- s1
      s1 <- s1.cg[rp2,2] <- mean(abs(rs1))	    # preserve history of AAD scale 
      g1.cg[rp2,2]  <- g1
      rt.cg[rp2,2]  <- Trt1
      if (abs(1-s1/s0) < cg1) break 
  }		

  for (i in 1:10) {     # robust estimation of gamma   #--------------------------------------(3)
      rp3 <- rp3 + 1; Trp <- Trp + 1

      # renew weights 
      u1 <- rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0

      ix1 <- which((w1*x1) !=0)     # remove data of w1*x1=0, since they make Trt1 and rs1 NaN 
      if (length(ix1) < 2) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))		
          ix1 <- which((w1*x1) !=0)
      }
      x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
      g1 <- (solve(t(x2*w1[ix1]) %*% (x2*w1[ix1])) %*% t(x2*w1[ix1]) %*% log(abs(rs1.x[ix1]*w1[ix1])))[2]   

      rs1    <- y1 / x1^g1- Trt1 * x1^(1-g1)	    # homoscedastic quasi residuals

      # M-estimation of ratio with fixed power (g1)
      s0 <- s1
      s1 <- s1.cg[rp2,2] <- mean(abs(rs1))	    # preserve history of AAD scale 
      g1.cg[rp2,2]  <- g1
      rt.cg[rp2,2]  <- Trt1
      if (abs(1-s1/s0) < cg1) break 
  }		

  # estimate power (g1) together with ratio     #---------------------------------------------(4)
  lrp <- rp.max - Trp
  for (i in 1:lrp) {
      rp4 <- rp4 + 1;   Trp <- Trp + 1
 
       # renew weights 
      u1 <- rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0

   #  ix1 <- which((x1*w1) !=0))
      ix1 <- which((x1*w1) !=0)
      if (length(ix1) < 2) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))		
         ix1 <- which((w1*x1) !=0)
      }
      x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
    #  g1<- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1])))[2]      # estimate power (g1)
      g1 <- (solve(t(x2*w1[ix1]) %*% (x2*w1[ix1])) %*% t(x2*w1[ix1]) %*% log(abs(rs1.x[ix1]*w1[ix1])))[2]     # estimate power (g1)

      ix1 <- which((w1*x1) !=0)     # remove data of w1*x1=0, since they make Trt1 and rs1 NaN
      if (length(ix1) < 2) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))		
          ix1 <- which((w1*x1) !=0)
      }
      (Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1)))) 
      rs1.x  <- y1 - Trt1 * x1		        # heteroscedastic residuals
      rs1    <- y1 / x1^g1- Trt1 * x1^(1-g1)	# homoscedastic quasi residuals
      s0 <- s1
      s1 <- s1.cg[rp2,2] <- mean(abs(rs1))	# preserve history of AAD scale
      s1 <- mean(abs(rs1))    

      u1 <-rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0
      ix1 <- which((w1*x1) !=0)     # removing data with w1*x1=0
      if (length(ix1) < 2) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))		
          ix1 <- which((w1*x1) !=0)
      }

      g1.cg[rp2,2]  <- g1
      rt.cg[rp2,2]  <- Trt1
      if (abs(1-s1/s0) < cg1) break 
  }
  if (abs(1-s1/s0) < cg1) fg1 <- 0 else fg1 <- 1
  return(list(rt=Trt1, g1=g1, wt=w1, rp=c(rp1, rp2, rp3, rp4, Trp), fg=fg1, rt.cg=rt.cg, g1.cg=g1.cg, s1.cg=s1.cg))
}
