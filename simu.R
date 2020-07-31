library(ggplot2)
library(csmodelmon)

model = list(order = c(0,1,1), ma=-0.8)
n=500
burn = 50

tests = c("adf.test","pp.test.mod","kpss.test.mod")

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
    
    lags_none <-  lag.select(simData, type = "none", kmax = 0.3*n)
    lags_drift <- lag.select(simData, type = "drift", kmax = 0.3*n)
    lags_trend <- lag.select(simData, type = "trend", kmax = 0.3*n)
    
      
    #perform test with no trend/const
    
    out_none[,,i] <-  sapply(lags, function(k){test <- urca::ur.df(simData, type = "none", selectlags = "Fixed", lags = as.numeric(k))
                                              test@teststat < test@cval[size]})
    #out_CADF[1,1,i] <- CADFtest::CADFtest(model = simData, type = "none", criterias = "MAIC")$p.value < 0.05
    
    #perform test with drift
    out_drift[1,,i] <- sapply(lags, function(k){test <- urca::ur.df(simData, type = "drift", selectlags = "Fixed", lags = as.numeric(k))
                                                test@teststat[,"tau2"] < test@cval["tau2",size]})
    out_drift[2,,i] <- sapply(lags, function(k){test <- urca::ur.pp(simData, model = "constant", type = "Z-tau", use.lag = as.numeric(k))
                                                test@teststat < test@cval[size]})
    out_drift[3,,i] <- sapply(lags, function(k){test <- urca::ur.ers(simData, model = "constant",type = "DF-GLS", lag.max = as.numeric(k))
                                                test@teststat < test@cval[size]})
    out_drift[4,,i] <- sapply(lags, function(k){test <- urca::ur.kpss(simData, type = "mu", use.lag =  as.numeric(k))
                                                test@teststat > test@cval[size]})
    #out_CADF[2,1,i] <- CADFtest::CADFtest(model = simData, type = "drift", criterias = "MAIC")$p.value < 0.05
    
    #perform test with trend
    out_trend[1,,i] <- sapply(lags, function(k){test <- urca::ur.df(simData, type = "trend", selectlags = "Fixed", lags = as.numeric(k))
                                                test@teststat[,"tau3"] < test@cval["tau3",size]})
    out_trend[2,,i] <- sapply(lags, function(k){test <- urca::ur.pp(simData, model = "trend", type = "Z-tau", use.lag = as.numeric(k))
                                               test@teststat < test@cval[size]})
    out_trend[3,,i] <- sapply(lags, function(k){test <- urca::ur.ers(simData, model = "trend",type = "DF-GLS", lag.max = as.numeric(k))
                                               test@teststat < test@cval[size]})
    out_trend[4,,i] <- sapply(lags, function(k){test <- urca::ur.kpss(simData, type = "tau", use.lag = as.numeric(k))
                                               test@teststat > test@cval[size]})
    #out_CADF[3,1,i] <- CADFtest::CADFtest(model = simData, type = "trend", criterias = "MAIC")$p.value < 0.05
    
  }

  resultTab <- rbind(apply(out_none, c(1,2), sum),apply(out_drift, c(1,2), sum),apply(out_trend, c(1,2), sum))
  apply(out_CADF, c(1,2), sum)
  return(resultTab)
}

I = 100
burn = 100
#Plain model (no drift, no trend)

  #Plain white noise
  sim25_000 <- baseSim(model = list(order = c(0,0,0)), n = 25, burn = burn, I = I)
  sim50_000 <- baseSim(model = list(order = c(0,0,0)), n = 50, burn = burn, I = I)
  sim100_000 <- baseSim(model = list(order = c(0,0,0)), n = 100, burn = burn, I = I)
  
  #Integrated model
  sim25_010 <- baseSim(model = list(order = c(0,1,0)), n = 25, burn = burn, I = I)
  sim50_010 <- baseSim(model = list(order = c(0,1,0)), n = 50, burn = burn, I = I)
  sim100_010 <- baseSim(model = list(order = c(0,1,0)), n = 100, burn = burn, I = I)
  
  #Integrated with large MA - positive
  sim25_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 25, burn = burn, I = I)
  sim50_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 50, burn = burn, I = I)
  sim100_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 100, burn = burn, I = I)
  #sim1000_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 1000, burn = burn, I = I)
  
  #Integrated with large MA - negative
  sim25_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 25, burn = burn, I = I)
  sim50_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 50, burn = burn, I = I)
  sim100_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 100, burn = burn, I = I)
  sim1000_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 1000, burn = burn, I = I)
  
  #Integrated with large AR
  sim25_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 25, burn = burn, I = I)
  sim50_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 50, burn = burn, I = I)
  sim100_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = burn, I = I)
  sim1000_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = burn, I = I)

#Model with drift 
  
  #Plain white noise
  sim25_000 <- baseSim(model = list(order = c(0,0,0)), n = 25, burn = burn, drift = 100, I = I)
  sim50_000 <- baseSim(model = list(order = c(0,0,0)), n = 50, burn = burn, drift = 100, I = I)
  sim100_000 <- baseSim(model = list(order = c(0,0,0)), n = 100, burn = burn, drift = 100, I = I)
  
  #Integrated model
  sim25_010 <- baseSim(model = list(order = c(0,1,0)), n = 25, burn = burn, drift = 100, I = I)
  sim50_010 <- baseSim(model = list(order = c(0,1,0)), n = 50, burn = burn, drift = 100, I = I)
  sim100_010 <- baseSim(model = list(order = c(0,1,0)), n = 100, burn = burn,drift = 100, I = I)
  
  #Integrated with large MA - positive
  sim25_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 25, burn = burn, drift = 100, I = I)
  sim50_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 50, burn = burn, drift = 100, I = I)
  sim100_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 100, burn = burn, drift = 100, I = I)
  sim1000_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 1000, burn = burn, drift = 100, I = I)
  
  #Integrated with large MA - negative
  sim25_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 25, burn = burn, drift = 100, I = I)
  sim50_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 50, burn = burn, drift = 100, I = I)
  sim100_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 100, burn = burn, drift = 100, I = I)
  sim1000_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 1000, burn = burn, drift = 100, I = I)
  
  #Integrated with large AR
  sim25_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 25, burn = burn, drift = 100, I = I)
  sim50_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 50, burn = burn, drift = 100, I = I)
  sim100_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = burn, drift = 100,  I = I)
  sim1000_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = burn, drift = 100, I = I)
  
#Model with drift + trend
  
  #Plain white noise
  sim25_000 <- baseSim(model = list(order = c(0,0,0)), n = 25, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim50_000 <- baseSim(model = list(order = c(0,0,0)), n = 50, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim100_000 <- baseSim(model = list(order = c(0,0,0)), n = 100, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  
  #Integrated model
  sim25_010 <- baseSim(model = list(order = c(0,1,0)), n = 25, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim50_010 <- baseSim(model = list(order = c(0,1,0)), n = 50, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim100_010 <- baseSim(model = list(order = c(0,1,0)), n = 100, burn = 10,drift = 100, trend = 1.1 * n, I = 1000)
  
  #Integrated with large MA - positive
  sim25_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 25, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim50_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 50, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim100_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 100, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  sim1000_010_MAP <- baseSim(model = list(order = c(0,1,1), ma = 0.8), n = 1000, burn = 10, drift = 100, trend = 1.1 * n, I = 1000)
  
  #Integrated with large MA - negative
  sim25_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 25, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  sim50_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 50, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  sim100_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 100, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  sim1000_010_MAN <- baseSim(model = list(order = c(0,1,1), ma = -0.8), n = 1000, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  
  #Integrated with large AR
  sim25_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 25, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  sim50_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 50, burn = 10, drift = 100, trend = 1.1 * n, I = 100)
  sim100_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = 10, drift = 100, trend = 1.1 * n,  I = 100)
  sim1000_010_ARP <- baseSim(model = list(order = c(1,1,0), ar = 0.9), n = 100, burn = 10, drift = 100, trend = 1.1 * n, I = 100)  
  
#Comparison between arima  
  theta <- c(0, -0.2, -0.5, -0.8, -0.99)
  phi   <- c(0, -0.2, -0.5, -0.8, -0.99)
  
  #whiten <- arima.sim(model = list(order = c(0,0,0)), n = 100)
  #ARIMA_010 <- arima.sim(model = list(order = c(0,1,0)), n = 100)
  
  ARIMA_010_MA <-  lapply(theta, function(i) sapply(1:100, function(x) arima.sim(model = list(order = c(0,1,1), ma = i), n = 100)) ) 
  ARIMA_010_AR <-  lapply(theta, function(i) sapply(1:100, function(x) arima.sim(model = list(order = c(1,1,0), ar = i), n = 100)) ) 
  ARIMA_010_MA_plots <-  ARIMA_010_AR_plots <- vector("list", length = length(theta))
 
  for(i in seq_along(ARIMA_010_MA)){
    
   meltdf <- reshape2::melt(cbind(as.data.frame(ARIMA_010_MA[[i]]), index = 1:101), id = "index")
   p <- ggplot(meltdf,aes(x=index,y=value,colour=variable,group=variable)) + geom_line() + 
        theme(legend.position = "none", 
              axis.title.y = element_blank(),
              axis.title.x = element_blank()) + 
     ggtitle(paste0("ARIMA(0,1,1) process, theta = ", theta[i]))
   print(p)
  }
  
  for(i in seq_along(ARIMA_010_AR)){
    
    meltdf <- reshape2::melt(cbind(as.data.frame(ARIMA_010_AR[[i]]), index = 1:101), id = "index")
    p <- ggplot(meltdf,aes(x=index,y=value,colour=variable,group=variable)) + geom_line() + 
      theme(legend.position = "none", 
            axis.title.y = element_blank(),
            axis.title.x = element_blank()) + 
      ggtitle(paste0("ARIMA(0,1,1) process, phi = ", theta[i]))
    print(p)
  }
  
plot(ARIMA_010)
plot(whiten)
plot(ARIMA_010_MA[[5]])
#Kmax
x <- 20:1000
kmaxCS         <- floor((x-3)/2)
kmax_Schwert4  <- 4 *(x/100)^(1/4)
kmax_Schwert12 <- 12 *(x/100)^(1/4)

kmax_log <- log(x)
plot(kmax_Schwert12)
plot(kmaxCS)
plot(kmax_log)
