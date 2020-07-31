MZtest <- function(x, k = 0, p=0){
  
  k <- k+1
  y <- diff(x)
  n <- length(y)
  z <- embed(y, k)
  yt <- z[,1]
  p=1
  tt <- (k:n)^p
  y_l1 <- x[k:n]
  
  
  if (k > 1){
    yt_j <- z[,2:k]
    
    if(p==0){
       rx <- lm(x[k:n] ~ 1)$residuals
       res  <- lm(yt ~ 1 + y_l1 + yt_j)
       B1   <- res$coefficients[3:length(res$coefficients)]
    }else{
      rx <- lm(x ~ 1:n)$residuals
      res  <- lm(yt ~ 1 + tt + y_l1 + yt_j)
      B1   <- res$coefficients[4:length(res$coefficients)]
     }
    B1   <- sum(B1)
  
  }else{
    if(p==0) 
      rx <- lm(x[1:n] ~ 1)$residuals
    res <- lm(yt ~ 1 + tt + y_l1)
    else res <- lm(yt ~ 1 + tt + y_l1)
    B1  <- 0
  }

  res.sum <- summary(res)
  e <- res$residuals
  TT <- n # sum(res.sum$df[1:2])
  sigma2_k <- sum(e^2)/TT
  s2_AR = sigma2_k/((1-B1)^2)
  p1 <- (x[n]^2)/n - s2_AR
  p2 <- (2/n^2) * sum(x[1:n]^2)
  STAT <-  p1/p2

  #STAT <- res.sum$coefficients[3, 1]/res.sum$coefficients[3, 2]
  
  table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), 
                 c(3.95, 3.8, 3.73, 3.69, 3.68, 3.66), 
                 c(3.6, 3.5, 3.45, 3.43,  3.42, 3.41), 
                 c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), 
                 c(1.14,1.19, 1.22, 1.23, 1.24, 1.25), 
                 c(0.8, 0.87, 0.9, 0.92, 0.93, 0.94), 
                 c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), 
                 c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
  
  table <- -table
  tablen <- dim(table)[2]
  tableT <- c(25, 50, 100, 250, 500, 1e+05)
  tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
  tableipl <- numeric(tablen)
  for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[, 
                                                            i], n, rule = 2)$y
  interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
  if (!is.na(STAT) && is.na(approx(tableipl, tablep, STAT, 
                                   rule = 1)$y)) 
    if (interpol == min(tablep)){ 
      warning("p-value smaller than printed p-value")
  }else{ warning("p-value greater than printed p-value")}

    PVAL <- interpol
  
  return(PVAL)

}

MZtest <- function(x, k = 0, p=0){
  
  # if ((NCOL(x) > 1) || is.data.frame(x)) 
  #   stop("x is not a vector or univariate time series")
  # if (any(is.na(x))) 
  #   stop("NAs in x")
  # if (k < 0) 
  #   stop("k negative")
  
  #alternative <- match.arg(alternative)
  #DNAME <- deparse(substitute(x))
  
  k <- k + 1
  x <- as.vector(x, mode = "double")
  y <- diff(x)
  n <- length(y)
  z <- embed(y, k)
  yt <- z[, 1]
  xt1 <- x[k:n]
  
  tt <- k:n
  rt <- 1:(n+1)
  
  rx <- lm(x ~ rt)
  rx <- rx$residuals
  
  if (k > 1) {
    yt1 <- z[, 2:k]
    res <- lm(yt ~ xt1 + 1 + tt + yt1)
  }
  
  res.sum <- summary(res)
  B1 <- sum(res.sum$coefficients[4:(k+2),1])
  
  e <- res$residuals
  sigma2 <- sum(e^2)/(n-k)
  s2_AR <- sigma2/(1-B1)^2
  term1 <- (rx[n]^2)/n - s2_AR
  term2 <- 2/(n^2) * sum(rx[k:n]^2)
  STAT <- term1/term2
  
  table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), c(3.95, 
                                                          3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 3.45, 3.43, 
                                                                                          3.42, 3.41), c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), c(1.14, 
                                                                                                                                                1.19, 1.22, 1.23, 1.24, 1.25), c(0.8, 0.87, 0.9, 0.92, 
                                                                                                                                                                                 0.93, 0.94), c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15, 
                                                                                                                                                                                                                                      0.24, 0.28, 0.31, 0.32, 0.33))
  table <- -table
  tablen <- dim(table)[2]
  tableT <- c(25, 50, 100, 250, 500, 1e+05)
  tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
  tableipl <- numeric(tablen)
  for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[, 
                                                            i], n, rule = 2)$y
  interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
  interpol
  
}

MSBtest <- function(x, k=0, p = 0){
  
  k <- k+1
  y <- diff(x)
  n <- length(y)
  z <- embed(y, k)
  yt <- z[,1]
  tt <- (k:n)^p
  y_l1 <- x[k:n]
  
  if (k > 1) {
    yt_j <- z[,2:k]
    res  <- lm(yt ~  1 + tt + y_l1 + yt_j)
    B1   <- res$coefficients[4:length(res$coefficients)]
    B1   <- sum(B1)
  }else{
    res <- lm(yt ~ 1 + tt + y_l1)
    B1 <- 0 
  }
  
  res.sum <- summary(res)
  e <- residuals(res)
  sigma2_k <- sum(e^2)/n
  
  s2_AR = sigma2_k/((1-B1)^2)
  
  STAT <- sqrt(1/n^2 * sum(x[1:n]^2)/s2_AR)
  STAT
  
}


kpss.test.mod(x, null="Trend")

simu <- do.call("rbind",lapply(1:1000, function(x) lag.select(arima.sim(model = list(order=c(0,1,0)),n=100), tpvl = 0.1)))
colMeans(simu)

apply(simul, 2, as.numeric)
lag.select(x, tpvl = 0.05)

#lapply(MAIC(x, 10), which.min)
MZtest(x,k=0)
MZtest(GLSdetrend(x,"constant"),k=0)
MSBtest(x,k=0)

MZtest(x,k=0)

MZtest(GLSdetrend(x,"trend"),k=1)
MZtest(GLSdetrend(x,"constant"),k=1)

sapply(1:1000, function(x) MSBtest(arima.sim(model = list(order=c(0,0,0)),n=1000)[50:550], k=0)) %>% mean()

sapply(1:1000, function(x) MZtest(arima.sim(model = list(order=c(0,1,0)),n=1000)[50:75], k=0) <0.05) %>% sum()
sapply(1:1000, function(x) MZtest(arima.sim(model = list(order=c(0,1,0)),n=1000)[50:100], k=0) <0.05) %>% sum()
sapply(1:1000, function(x) MZtest(arima.sim(model = list(order=c(0,1,0)),n=1000)[50:550], k=0) <0.05) %>% sum()
sapply(1:1000, function(x) MZtest(arima.sim(model = list(order=c(0,1,0)),n=10000)[50:10000], k=0) <0.05) %>% sum()
cat("-----------------------------\n")

sapply(1:1000, function(x) adf.test(arima.sim(model = list(order=c(0,1,0)),n=1000)[50:550], k=0)$p.value <0.05) %>% sum()

sapply(1:1000, function(x) adf.test(cumsum(rnorm(100)))$p.value <0.05) %>% sum()


test_res <- MZtest(y_fala, k=1)
x<-y_fala
k=3

ur.ers(x)
