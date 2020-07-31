kpss.test.mod <- function (x, null = c("Level", "Trend"), k = ifelse(length(x) > 100, trunc(12 * (n/100)^0.25), trunc(4 * (n/100)^0.25))){
  
  if ((NCOL(x) > 1) || is.data.frame(x)) 
    stop("x is not a vector or univariate time series")
  DNAME <- deparse(substitute(x))
  null <- match.arg(null)
  x <- as.vector(x, mode = "double")
  n <- length(x)
  if (null == "Trend") {
    t <- 1:n
    e <- residuals(lm(x ~ t))
    table <- c(0.216, 0.176, 0.146, 0.119)
  }
  else if (null == "Level") {
    e <- residuals(lm(x ~ 1))
    table <- c(0.739, 0.574, 0.463, 0.347)
  }
  tablep <- c(0.01, 0.025, 0.05, 0.1)
  s <- cumsum(e)
  eta <- sum(s^2)/(n^2)
  s2 <- sum(e^2)/n
  s2 <- .C("tseries_pp_sum", as.vector(e, mode = "double"), as.integer(n), 
           as.integer(k), s2 = as.double(s2), PACKAGE = "tseries")$s2
  STAT <- eta/s2
  PVAL <- approx(table, tablep, STAT, rule = 2)$y
  if (!is.na(STAT) && is.na(approx(table, tablep, STAT, rule = 1)$y)) 
    if (PVAL == min(tablep)) 
      warning("p-value smaller than printed p-value")
  else warning("p-value greater than printed p-value")
  PARAMETER <- k
  METHOD <- paste("KPSS Test for", null, "Stationarity")
  names(STAT) <- paste("KPSS", null)
  names(PARAMETER) <- "Truncation lag parameter"
  structure(list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
                 method = METHOD, data.name = DNAME), class = "htest")
}