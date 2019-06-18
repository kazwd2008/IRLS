#################################################################################################
#################################################################################################
#  Robust estimators for the generalised ratio model　
#         by iteratively re-weighted least squares (IRLS) algorithm 
#    Weight function:  Tukey's biweight function
#    Scale: Average absolute deviation (AAD) or median absolute deviation (MAD)
#------------------------------------------------------------------------------------------------#
#       Created from Tirls.r  Ver. 4.0 [2014/09/05] for regression model                      
#       Ver.0     [2017/03/27]  Functions with AAD scale are implemented. 
#       Ver.1     [2017/04/10]  - Functions with MAD scale are added.
#               		- An error flag added, which become 1 and abort the function 
#	                          when all weights (w1) become zero.
#------------------------------------------------------------------------------------------------#
#  Functions  
# 	RrTa: 	  Gamma = 1,   AAD scale  　　　　　　                     
# 	RrTb:     Gamma = 1/2, AAD scale (conventional ratio model)
# 	RrTc:     Gamma = 0,   AAD scale (single regression without intercept)
# 	RrTa.mad: Gamma = 1,   MAD scale* 
# 	RrTb.mad: Gamma = 1/2, MAD scale* (conventional ratio model)
# 	RrTc.mad: Gamma = 0,   MAD scale* (single regression without intercept)
#           ※ Since the mad function in R returns values corresponding to the standard deviation
#              (SD), tuning constants for SD are used instead of those for MAD.    
#------------------------------------------------------------------------------------------------#
#  Parameters 
#   x1      single explanatory variable 
#   y1      objective variable                               
#   c1      tuning parameter for Tukey's biweight function
#           (default setting : c1=8 for AAD)。 
#   c=4     very robust
#   c=8     less robust
#   dat     name of dataframe in which x1 and y1 are included
#   rp.max  maximum number of iteration (default setting : 100)
#------------------------------------------------------------------------------------------------#
#  Returned values
#   par     robustly estimated rate y1/x1 
#   wt      robust weights
#   rp      total number of iteration
#   s1      changes of the scale (AAD or MAD) 
#   efg	    error flag, 1: acalculia (all weights become zero)  0: successful termination
#################################################################################################
#################################################################################################
# 	RrTa :   gamma = 1
#------------------------------------------------------------------------------------------------#
RrTa <- function(x1, y1, c1=8, dat="", rp.max=10, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow
  if (dat!="") attach(get(dat))	        
  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0
  par <- mean(y1 / x1)	                     	 # initial estimation
  res <- y1 / x1 - par			         # homoscedastic quasi-residuals  
  rp1 <- 1				         # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale 

  #### calculating weights
  	u1 <-res/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0
	if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {                            
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)     	            
      res <- y1 / x1 - par			# homoscedastic quasi-residuals
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <- mean(abs(res))	# AAD scale
      u1 <-res/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0
      if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1	
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
RrTa.mad <- function(x1, y1, c1=10.03, dat="", rp.max=50, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow
  if (dat!="") attach(get(dat))	        
  s1.cg <- rep(0, rp.max)               	 # preserve changes in s1 (scale)
  efg <- 0
  par <- mean(y1 / x1)	                	 # initial estimation
  res <- y1 / x1 - par				 # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)    	 # MAD scale

  #### calculating weights
  	u1 <-res/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0
	if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {                            
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)     	            
      res <- y1 / x1 - par			# homoscedastic quasi-residuals 
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <- mad(res)	        # MAD scale
      u1 <-res/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0
      if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1	
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RetTb :   gamma = 1/2 (conventional ratio estimator)
#------------------------------------------------------------------------------------------------#
RrTb <- function(x1, y1, c1=8, dat="", rp.max=10, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1) # prevent overflow
  if (dat!="") attach(get(dat))	      
  s1.cg <- rep(0, rp.max)                	# preserve changes in s1 (scale)
  efg <- 0
  par <- sum(y1) / sum(x1)     	      
  res <- y1/sqrt(x1) - par*sqrt(x1)	 	# homoscedastic quasi-residuals 
  rp1 <- 1					# number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))      # AAD scale

  #### calculating weights
 　u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0)  return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * w1) / sum(x1 * w1 )    
        res   <-  y1 / sqrt(x1) - par * sqrt(x1)  # homoscedastic quasi-residuals 
        rp1 <- rp1 + 1				  # number of iteration
        s1 <- s1.cg[rp1] <- mean(abs(res))        # AAD scale 
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break           # convergence condition
        s0 <- s1	
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
RrTb.mad <- function(x1, y1, c1=10.03, dat="", rp.max=50, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow
  if (dat!="") attach(get(dat))	      
  s1.cg <- rep(0, rp.max)                # preserve changes in s1 (scale)
  efg <- 0
  par <- sum(y1) / sum(x1)     	         # initial estimation
  res <- y1/sqrt(x1) - par*sqrt(x1)	 # homoscedastic quasi-residuals 
  rp1 <- 1				 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)     # MAD scale

  #### calculating weights
   u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * w1) / sum(x1 * w1 )     
        res   <-  y1 / sqrt(x1) - par * sqrt(x1)	# homoscedastic quasi-residuals 
        rp1 <- rp1 + 1					# number of iteration
        s1 <- s1.cg[rp1] <- mad(res)                    # MAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                 # convergence condition
        s0 <- s1	
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RrTc :   gamma = 1 (a single regression without intercept)
#------------------------------------------------------------------------------------------------#
RrTc <- function(x1, y1, c1=8, dat="", rp.max=10, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow
  if (dat!="") attach(get(dat))	      
  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0
  par <- sum(y1 * x1) / sum(x1**2)               # initial estimation
  res <- y1 - par*x1	                         # homoscedastic quasi-residuals 
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale

  #### calculating weights
   u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )    
        res   <-  y1  - par * x1		          # homoscedastic quasi-residuals 
        rp1 <- rp1 + 1					  # number of iteration
        s1 <- s1.cg[rp1] <- mean(abs(res))          	  # AAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                   # convergence condition
        s0 <- s1	
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
RrTc.mad <- function(x1, y1, c1=10.03, dat="", rp.max=50, cg.rt=0.01) {
  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow
  if (dat!="") attach(get(dat))	      
  s1.cg <- rep(0, rp.max)                   	 # preserve changes in s1 (scale)
  efg <- 0
  par <- sum(y1 * x1) / sum(x1**2)     	    	 # initial estimation
  res <- y1 - par*x1	                    	 # homoscedastic quasi-residuals 
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)        	 # MAD scale

  #### calculating weights
　u1 <-res/(c1*s1)
  w1 <- (1-u1**2)**2
  w1[which(abs(u1)>=1)] <- 0
  if (sum(w1)==0)  return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration 
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )     
        res   <-  y1  - par * x1		        # homoscedastic quasi-residuals 
        rp1 <- rp1 + 1					# number of iteration
        s1 <- s1.cg[rp1] <- mad(res)                    # MAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                 # convergence condition
        s0 <- s1	
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }
#################################################################################################
#################################################################################################
