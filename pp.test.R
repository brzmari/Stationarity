pp.test.mod <-function (x, alternative = c("stationary", "explosive"), type = c("Z(alpha)", 
                                                                  "Z(t_alpha)"), k = ifelse(length(x) > 100, trunc(12 * (n/100)^0.25), trunc(4 * (n/100)^0.25))) 
{
  if ((NCOL(x) > 1) || is.data.frame(x)) 
    stop("x is not a vector or univariate time series")
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  x <- as.vector(x, mode = "double")
  z <- embed(x, 2)
  yt <- z[, 1]
  yt1 <- z[, 2]
  n <- length(yt)
  tt <- (1:n) - n/2
  res <- lm(yt ~ 1 + tt + yt1)
  if (res$rank < 3) 
    stop("Singularities in regression")
  res.sum <- summary(res)
  u <- residuals(res)
  ssqru <- sum(u^2)/n
  ssqrtl <- .C(tseries_pp_sum, as.vector(u, mode = "double"), 
               as.integer(n), as.integer(k), ssqrtl = as.double(ssqru))$ssqrtl
  n2 <- n^2
  trm1 <- n2 * (n2 - 1) * sum(yt1^2)/12
  trm2 <- n * sum(yt1 * (1:n))^2
  trm3 <- n * (n + 1) * sum(yt1 * (1:n)) * sum(yt1)
  trm4 <- (n * (n + 1) * (2 * n + 1) * sum(yt1)^2)/6
  Dx <- trm1 - trm2 + trm3 - trm4
  if (type == "Z(alpha)") {
    alpha <- res.sum$coefficients[3, 1]
    STAT <- n * (alpha - 1) - (n^6)/(24 * Dx) * (ssqrtl - 
                                                   ssqru)
    table <- cbind(c(22.5, 25.7, 27.4, 28.4, 28.9, 29.5), 
                   c(19.9, 22.4, 23.6, 24.4, 24.8, 25.1), c(17.9, 19.8, 
                                                            20.7, 21.3, 21.5, 21.8), c(15.6, 16.8, 17.5, 
                                                                                       18, 18.1, 18.3), c(3.66, 3.71, 3.74, 3.75, 3.76, 
                                                                                                          3.77), c(2.51, 2.6, 2.62, 2.64, 2.65, 2.66), 
                   c(1.53, 1.66, 1.73, 1.78, 1.78, 1.79), c(0.43, 0.65, 
                                                            0.75, 0.82, 0.84, 0.87))
  }
  else if (type == "Z(t_alpha)") {
    tstat <- (res.sum$coefficients[3, 1] - 1)/res.sum$coefficients[3, 
                                                                   2]
    STAT <- sqrt(ssqru)/sqrt(ssqrtl) * tstat - (n^3)/(4 * 
                                                        sqrt(3) * sqrt(Dx) * sqrt(ssqrtl)) * (ssqrtl - ssqru)
    table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), 
                   c(3.95, 3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 
                                                           3.45, 3.43, 3.42, 3.41), c(3.24, 3.18, 3.15, 
                                                                                      3.13, 3.13, 3.12), c(1.14, 1.19, 1.22, 1.23, 
                                                                                                           1.24, 1.25), c(0.8, 0.87, 0.9, 0.92, 0.93, 0.94), 
                   c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15, 0.24, 
                                                           0.28, 0.31, 0.32, 0.33))
  }
  else stop("irregular type")
  table <- -table
  tablen <- dim(table)[2]
  tableT <- c(25, 50, 100, 250, 500, 1e+05)
  tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
  tableipl <- numeric(tablen)
  for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[, 
                                                            i], n, rule = 2)$y
  interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
  if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y)) 
    if (interpol == min(tablep)) 
      warning("p-value smaller than printed p-value")
  else warning("p-value greater than printed p-value")
  if (alternative == "stationary") 
    PVAL <- interpol
  else if (alternative == "explosive") 
    PVAL <- 1 - interpol
  else stop("irregular alternative")
  PARAMETER <- k
  METHOD <- "Phillips-Perron Unit Root Test"
  names(STAT) <- paste("Dickey-Fuller", type)
  names(PARAMETER) <- "Truncation lag parameter"
  structure(list(statistic = STAT, parameter = PARAMETER, alternative = alternative, 
                 p.value = PVAL, method = METHOD, data.name = DNAME), 
            class = "htest")
}