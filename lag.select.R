lag.select <- function(x, kmax = floor(20*(length(x)/100)^(1/4)), type = "none", tpvl = 0.05){

  out <- vector("numeric",kmax+1)
  kMAIC <- kBIC <- kAIC <- kHQC <- kSBC <- kT <- out
  
  x <- as.vector(x, mode = "double")
  y <- diff(x)
  n <- length(y)

  kmaxCS     <- min(8, floor((length(x)-3)/2))
  csMod_Adf  <- kmaxCS
  csMod_KPSS <- min(trunc(10*sqrt(n)/14), kmax)
  
  if(type == "none")  f <- "yt ~ xt1 + -1"
  if(type == "drift") f <- "yt ~ xt1 + 1"
  if(type == "trend") f <- "yt ~ xt1 + 1 + tt"
  
  strtT <- kmax + 2
  endT <- length(x)
  
  k <- kmax + 1
  z <- embed(y, k)
  yt <- z[, 1]
  xt1 <- x[k:n]
  tt <- k:n
  
  for(k in 0:kmax){

    #without drift and trend
    if (k > 0) {
      yt1 <- z[, 2:(k+1)]
      form <- paste0(f, "+ yt1")
    }else{
      form <- f
    }
    
    res <- lm(as.formula(form))
    
    res.sum <- summary(res)
    p  <- res.sum$df[1]
    TT <- n - kmax #p + res.sum$df[2] #n - kmax#
    e  <- residuals(res)
    b0 <- res$coefficients["xt1"]
    #stt = kmax - k + 1
    sigma2_k = sum(e^2)/TT
    
    if(type=="drift") xf <- xt1 - mean(xt1)
    if(type=="trend") xf <- lsfit(1:TT, xt1)$residuals
    else xf <- xt1
    
    tau <- b0^2 * sum(xf^2)/sigma2_k
    kAIC[k+1]  = log(sigma2_k) + 2*p/TT
    kBIC[k+1]  = log(sigma2_k) + p*log(TT)/TT
    kSBC[k+1]  = log(sigma2_k) + log(TT)*2*p/TT
    kHQC[k+1]  = log(sigma2_k) + 2*p*log(log(TT))/TT
    kMAIC[k+1] = log(sigma2_k) + 2*(tau+p)/TT

    if(k==0) kT[k+1] = 0 else kT[k+1] = res.sum$coefficients[p,4]
  } 
  
  if(all(kT<=tpvl)==TRUE){
    STG = GTS = kmax
  }else if(any(kT<=tpvl)==TRUE){
    STG = min(which((kT<=tpvl)==!TRUE))-2
    GTS = max(which((kT<=tpvl)==TRUE))-1
  }else{
    STG = GTS = 0
  }
    
  return(c(MAIC = which.min(kMAIC)-1,
           BIC =  which.min(kBIC)-1,
           AIC =  which.min(kAIC)-1,
           HQC =  which.min(kHQC)-1,
           STG = STG,
           GTS =  GTS,
           csMod_Adf = csMod_Adf,
           csMod_KPSS = csMod_KPSS))
  return(out)  
}
