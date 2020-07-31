baseSim <- function(model, n, I = 100, size = 0.05, drift = 0, trend = NULL,
                    burn = 10, kmax = floor(12*(n/100)^(1/4))){
  
  sizekp = if(size==0.1) 1 else if(size==0.05) 2 else if (size==0.01) 4
  size = if(size==0.1) 3 else if(size==0.05) 2 else if (size==0.01) 1
  
  #Create df for storing the results
  criterias <- c("MAIC","BIC","AIC","HQC","STG","GTS", "csMod_Adf","csMod_KPSS")
  tests <- c("ADF", "PP", "DF-GLS", "KPSS")
  
  #out <- matrix(nrow = length(tests), ncol = 8, dimnames = list(tests, criterias))
  out_none <- array(NA, dim= c(1, 8, I), dimnames = list("ADF", criterias, 1:I))
  out_drift <- array(NA, dim= c(4, 8, I), dimnames = list(tests, criterias, 1:I))
  out_trend <- array(NA, dim= c(4, 8, I), dimnames = list(tests, criterias, 1:I))
  out_CADF <- array(NA, dim= c(3, 1, I), dimnames = list(c("none","drift","trend"), "CADF_MAIC", 1:I))
  
  for(i in seq(I)){
    
    #Generate the data
    simData <- arima.sim(model, n = n + burn)
    simData <- simData[(burn+1):(n+burn)]
    if(drift!=0) simData <- simData + drift
    if(!is.null(trend)){t <- 1:n
    simData <- simData + eval(parse(text = trend))}
    
    lags_none <-  lag.select(simData, type = "none", kmax = 10)
    lags_drift <- lag.select(simData, type = "drift", kmax = 10)
    lags_trend <- lag.select(simData, type = "trend", kmax = 10)
    
    
    #perform test with no trend/const
    
    out_none[,,i] <-  sapply(lags_none, function(k){test <- urca::ur.df(simData, type = "none", selectlags = "Fixed", lags = as.numeric(k))
    test@teststat < test@cval[size]})
    #out_CADF[1,1,i] <-  CADFTEST(model = simData, type = "none", criterion = "BIC", max.lag.y=10) 
    
    # estmodel(model = simData, type = "none", criterion = "BIC", max.lag.y=10, X=NULL)
    #perform test with drift
    out_drift[1,,i] <- sapply(lags_drift, function(k){test <- urca::ur.df(simData, type = "drift", selectlags = "Fixed", lags = as.numeric(k))
    test@teststat[,"tau2"] < test@cval["tau2",size]})
    out_drift[2,,i] <- sapply(lags_drift, function(k){test <- urca::ur.pp(simData, model = "constant", type = "Z-tau", use.lag = as.numeric(k))
    test@teststat < test@cval[size]})
    out_drift[3,,i] <- sapply(lags_drift, function(k){test <- urca::ur.ers(simData, model = "constant",type = "DF-GLS", lag.max = as.numeric(k))
    test@teststat < test@cval[size]})
    out_drift[4,,i] <- sapply(lags_drift, function(k){test <- urca::ur.kpss(simData, type = "mu", use.lag =  as.numeric(k))
    test@teststat > test@cval[size]})
    #out_CADF[2,1,i] <- CADFtest::CADFtest(model = simData, type = "drift", criterias = "MAIC")$p.value < 0.05
    
    #perform test with trend
    out_trend[1,,i] <- sapply(lags_trend, function(k){test <- urca::ur.df(simData, type = "trend", selectlags = "Fixed", lags = as.numeric(k))
    test@teststat[,"tau3"] < test@cval["tau3",size]})
    out_trend[2,,i] <- sapply(lags_trend, function(k){test <- urca::ur.pp(simData, model = "trend", type = "Z-tau", use.lag = as.numeric(k))
    test@teststat < test@cval[size]})
    out_trend[3,,i] <- sapply(lags_trend, function(k){test <- urca::ur.ers(simData, model = "trend",type = "DF-GLS", lag.max = as.numeric(k))
    test@teststat < test@cval[size]})
    out_trend[4,,i] <- sapply(lags_trend, function(k){test <- urca::ur.kpss(simData, type = "tau", use.lag = as.numeric(k))
    test@teststat > test@cval[size]})
    #out_CADF[3,1,i] <- CADFtest::CADFtest(model = simData, type = "trend", criterias = "MAIC")$p.value < 0.05
    
  }
  
  resultTab <- rbind(apply(out_none, c(1,2), sum),apply(out_drift, c(1,2), sum),apply(out_trend, c(1,2), sum))
  apply(out_CADF, c(1,2), sum)
  return(resultTab)
}