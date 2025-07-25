###############################################################################################
#  RBred:  Robust estimator for the generalized ratio model   
#          	by iteratively re-weighted least squares (IRLS) algorithm for M-estimation  
#          Weight function:  		Tukey's biweight function
#          Scale of quasi-residuals:	Average absolute deviation (AAD)
#---------------------------------------------------------------------------------------------#
#       Implemented by K. Wada (NSTAC, Japan)	
#---------------------------------------------------------------------------------------------#
#       Ver. 0.0 [2016.07.22]  
#       Ver. 0.1 [2018.05.21]  A bug to replace the value of scale parameter was fixed. #A#
#       Ver. 0.2 [2018.05.21]  The iteration count improved.
#       Ver. 0.3 [2019.07.08]  b1 => g1
#       Ver. 0.4 [2019.07.25]  correct the last estimation of g1 with weights
# 	Ver. 1.0 [2019.09.30]  publish for JSM2019 Proceedings after a minor bug correction
#                              and improved loop counters [See Wada, Takata, and Tsubaki (2019)]
#       Ver. 1.1 [2025.07.26]  An internal bug (which does not affect outcomes) corrected
#---------------------------------------------------------------------------------------------#
#  Parameters
#   x1       single explanatory variable
#   y1       objective variable to be imputed
#   g1       initial value of power (default: g1=0)
#   c1       tuning constant for the biweight function
#            the value should be 4 to 8 (default: c=8) based on Bienias et al. (1997)
#               [ c=4 :   very robust,      c=8 :    less robust]
#   dat      name of the dataframe containing x1 and y1 if any 
#   rp.max   maximum number of iteration (default: rp.max=100) 
#   cg1      convergence condition to stop iteration (default: cg1=0.001)
#---------------------------------------------------------------------------------------------#
#  Return values
#   rt       estimated robust ratio of y1 to x1 	(beta)
#   g1       estimated robust power 			(gamma)
#   wt       robust weight 
#   rp       number of iterations (four elements) 
#   fg       convergence flag (fg=0: converged; fg=1: not converged）
#   rt.cg    histry of the rate 
#   g1.cg    history of the power
#   s1.cg    history of the scale parameter
################################################################################################

RBred <- function(x1, y1, g1=0, c1=8, dat="", rp.max=100,  cg1=0.001){

 # initial estimation
  if (dat!="") attach(get(dat))
  s1.cg <- rt.cg <- g1.cg  <- rep(NA, (rp.max*4))	# for preserving histories
  rp1     <- c(1, 0, 0, 0)  				# count of iteration for each step
  s0      <- 0	    					# initial value of scale parameter

####---------------------------------------------------------------------------------------(1)
   Trt1 <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))	# ratio estimation
  # rs1.x <- y1 - Trt1 * x1				# heteroscedastic residuals 
   rs1.x <- y1 * x1^g1 - Trt1 * x1^(1+g1)		# heteroscedastic residuals (2025.07.26)
   rs1   <- y1 / x1^g1 - Trt1 * x1^(1-g1)		# homoscedastic quasi-residuals

   s1 <- s1.cg[rp1[1]] <- mean(abs(rs1))	      　# AAD scale of quasi-residuals
   g1.cg[rp1[1]]  <- g1
   rt.cg[rp1[1]]  <- Trt1
####---------------------------------------------------------------------------------------
#### Non-robust Estimation of rate (Trt1) and power (g1) (without weight)
####---------------------------------------------------------------------------------------
  for (i in 2:rp.max) {
       rp1[1] <- rp1[1] + 1
       ix1 <- which((x1 !=0) & (rs1 !=0))
       x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
       g1 <- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1])))[2]   # estimate power

       Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))		   # estimate ratio
       # rs1.x  <- y1 - Trt1 * x1		# heteroscedastic residuals (conventional)
       rs1.x  <- y1 * x1^g1 - Trt1 * x1^(1+g1)	# heteroscedastic residuals (2025.07.26)
       rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals 

       s0 <- s1					# preserve the previous value (AAD)
       s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))	# AAD of current quasi residuals
       g1.cg[sum(rp1)]  <- g1
       rt.cg[sum(rp1)]  <- Trt1
       if (abs(1-s1/s0) < cg1) break 		# convergence condition
  }

####---------------------------------------------------------------------------------------(2)
####  Iterative robust estimation of the ratio "Trt1" 
####                               with a fixed power "g1" derived in step(1)
####---------------------------------------------------------------------------------------
  for (i in 1:rp.max) {
      rp1[2] <- rp1[2] + 1

    # renew weights 
　    u1 <-rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0

    # prevent all zero weights
      ix1 <- which((w1*x1) !=0)  & (rs1.x !=0))    # remove observations that make Trt1 and rs1 NaN
      if (length(ix1)==0) {         # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))
          ix1 <- which((w1*x1) !=0)
      }

    # estimate Trt1 with fixed value of g1 obtained in the previous step
      Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
      # rs1.x  <- y1 - Trt1 * x1		      # heteroscedastic residuals
      rs1.x <- y1 * x1^g1 - Trt1 * x1^(1+g1)	    # heteroscedastic residuals (2025.07.26)
      rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	    # homoscedastic quasi-residuals

      s0 <- s1					    # preserve previous AAD value
      s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))	    # preserve history of AAD scale
      g1.cg[sum(rp1)]  <- g1
      rt.cg[sum(rp1)]  <- Trt1
      if (abs(1-s1/s0) < cg1) break 		　　# convergence condition
  }

####---------------------------------------------------------------------------------------(3)
####  Iterative robust and simultaneous estimation of the ratio "Trt1" and the power "g1"
####---------------------------------------------------------------------------------------

  for (i in 1:rp.max) {
      rp1[3] <- rp1[3] + 1

    # renew weights 
　    u1 <-rs1/(c1*s1)			# added 2019.09.27
      w1 <- (1-u1**2)**2		# added 2019.09.27
      w1[which(abs(u1)>=1)] <- 0	# added 2019.09.27

      ix1 <- which(((w1*x1)!=0) & (rs1.x !=0))  # remove observations that make Trt1 and rs1 NaN
      if (length(ix1)==0) {                     # reset w1 as all 1 when all the weights become zero 
          w1 <- rep(1, length(x1))
          ix1 <- which((w1*x1) !=0)
      }
      x2  <- cbind(rep(1, length(ix1)), log(x1[ix1]*w1[ix1]))
      (g1 <- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1]*w1[ix1])))[2])     # estimate power (g1)

      Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
      # rs1.x  <- y1 - Trt1 * x1			# heteroscedastic residuals
      rs1.x <- y1 * x1^g1 - Trt1 * x1^(1+g1)	# heteroscedastic residuals (2025.07.26)
      rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals

      s0 <- s1
      s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))
      g1.cg[sum(rp1)]  <- g1
      rt.cg[sum(rp1)]  <- Trt1

      if (abs(1-s1/s0) < cg1) break 		# convergence condition
    }

####---------------------------------------------------------------------------------------(4)
####  If the estimation is not converge, reduce tuning constant slightly
####---------------------------------------------------------------------------------------

  if (abs(1-s1/s0) >= cg1) {　　# when convergence condition is not met
     c1 <- c1 - 0.1		# reduce tuning constant slightly
     for (i in 1:rp.max) {
         rp1[4] <- rp1[4] + 1

       # renew weights 
　       u1 <-rs1/(c1*s1)
         w1 <- (1-u1**2)**2
         w1[which(abs(u1)>=1)] <- 0

         ix1 <- which((x1 !=0) & (rs1.x !=0))
         if (length(ix1)==0) {         # reset w1 as all 1 when all the weights become zero 
            w1 <- rep(1, length(x1))		
            ix1 <- which((w1*x1) !=0)
         }
         x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]*w1[ix1]))
         (g1 <- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1]*w1[ix1])))[2])

         Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
         # rs1.x  <- y1 - Trt1 * x1			# heteroscedastic residuals
         rs1.x  <- y1 * x1^g1 - Trt1 * x1^(1+g1)	# heteroscedastic residuals (2025.07.26)
         rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals

         s0 <- s1
         s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))
         g1.cg[sum(rp1)]  <- g1
         rt.cg[sum(rp1)]  <- Trt1

      if (abs(1-s1/s0) < cg1) break			# convergence condition
    }
  }
    if (abs(1-s1/s0) < cg1) fg1 <- 0 else fg1 <- 1
    return(list(rt=Trt1, g1=g1, wt=w1, rp=rp1, fg=fg1, rt.cg=rt.cg, g1.cg=g1.cg, s1.cg=s1.cg))
}

###############################################################################################
#  Bred:  Estimator for the generalized ratio model   (not robust)
#          	by iteratively re-weighted least squares (IRLS) algorithm for M-estimation  
#          Scale of quasi-residuals:	Average absolute deviation (AAD)
#---------------------------------------------------------------------------------------------#
#       Implemented by K. Wada (NSTAC, Japan)
#---------------------------------------------------------------------------------------------#
#       Ver. 0.0 [2019.09.30]  Prepared for JSM2019 Proceedings 
# 				[See Wada, Takata, and Tsubaki (2019)]
#---------------------------------------------------------------------------------------------#
#  Parameters
#   x1       single explanatory variable
#   y1       objective variable to be imputed
#   g1       initial value of power (default: g1=0)
#   dat      name of the dataframe containing x1 and y1 if any 
#   rp.max   maximum number of iteration (default: rp.max=100) 
#   cg1      convergence condition to stop iteration (default: cg1=0.001)
#---------------------------------------------------------------------------------------------#
#  Return values
#   rt       estimated robust ratio of y1 to x1 	(beta)
#   g1       estimated robust power 			(gamma)
#   rp       number of iterations (four elements) 
#   fg       convergence flag (fg=0: converged; fg=1: not converged）
#   rt.cg    histry of the rate 
#   g1.cg    history of the power
#   s1.cg    history of the scale parameter
################################################################################################

Bred <- function(x1, y1, g1=0, dat="", rp.max=100,  cg1=0.001){

 # initial estimation
  if (dat!="") attach(get(dat))
  s1.cg <- rt.cg <- g1.cg  <- rep(NA, rp.max)	# for preserving histories
  rp1     <- 1  				# Number of iteration
  s0      <- 0	    				# Initial value for AAD scale

####---------------------------------------------------------------------------------------(1)
   (Trt1 <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1))))	# estimate ratio
   # rs1.x <- y1 - Trt1 * x1
   rs1.x <- y1 * x1^g1 - Trt1 * x1^(1+g1)		# heteroscedastic residuals (2025.07.26)
   rs1   <- y1 / x1^g1 - Trt1 * x1^(1-g1)

   g1.cg[rp1]  <- g1
   rt.cg[rp1]  <- Trt1
   s1 <- s1.cg[rp1] <- mean(abs(rs1))
####---------------------------------------------------------------------------------------(2)
#### Iterative estimation of Trt1 and g1
####---------------------------------------------------------------------------------------
  for (i in 2:rp.max) {
       rp1 <- rp1 + 1
       ix1 <- which((x1 !=0) & (rs1 !=0))
       x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
       g1 <- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1])))[2]

       (Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1))))
       # rs1.x  <- y1 - Trt1 * x1
       rs1.x  <- y1 * x1^g1 - Trt1 * x1^(1+g1)		# heteroscedastic residuals (2025.07.26)
       rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)

       s0 <- s1
       s1 <- s1.cg[rp1] <- mean(abs(rs1))
       g1.cg[rp1]  <- g1
       rt.cg[rp1]  <- Trt1
       if (abs(1-s1/s0) < cg1) break
  }
  if (abs(1-s1/s0) < cg1) fg1 <- 0 else fg1 <- 1
  return(list(rt=Trt1, g1=g1, rp=rp1, fg=fg1, rt.cg=rt.cg, g1.cg=g1.cg, s1.cg=s1.cg))
}
###############################################################################################
# References
#---------------------------------------------------------------------------------------------#
# Bienias JL, Lassman DM, Scheleur SA, Hogan H (1997). “Improving Outlier Detection
#    in Two Establishment Surveys.” In Statistical Data Editing Volume No.2 Methods and
#    Techniques, pp. 76?83. UNSC/UNECE. 
#
# Wada K, Takata S, Tsubaki H, (2019) An Algorithm of Generalized Robust Ratio Model 
#    Estimation for Imputation, to be appeared in JSM2019 Proceedings.
###############################################################################################
