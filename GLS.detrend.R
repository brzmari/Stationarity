GLSdetrend <- function(x,model = c("trend", "constant")){
  
  n <- length(x)
  z <- embed(x,2)
  y <- z[,2]
  y_l1 <- z[,1]
  
  if(model=="trend"){
    
    tt <- 1:n
    cp = -13.5
    alpha <- 1 + (cp/n)
    y_alpha_c <- c(y[1], y - alpha * y_l1)
    z_alpha_c <- c(1, rep(1 - alpha, n - 1))
    ztrend <- embed(tt,2)
    y_alpha_t <- c(1, ztrend[,1] - alpha * ztrend[,2])
    z_alpha_t <- c(1, ztrend[,1] - alpha * ztrend[,2])
    df.reg <- lm(y_alpha ~ -1 + z_alpha_c + z_alpha_t)
    
    y_detrend  <- x - coef(df.reg)[1] - coef(df.reg)[2] * tt
    
  }else if(model=="constant"){
    
    alpha <- 1 + (cp/n)
    y_alpha <- c(y[1], y - alpha * y_l1)
    z_alpha <- c(1, rep(1 - alpha, n - 1))
    df.reg <- lm(y_alpha ~ -1 + z_alpha)
    
    y_detrend <- x - alpha * df.reg$coefficient
    
  }
  
  return(y_detrend)
}
