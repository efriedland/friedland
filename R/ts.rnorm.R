function (n, mean = 0, sd = 1, initial, allownegatives = TRUE) {
  now <- as.POSIXlt(Sys.time())
  r <- initial + cumsum(rnorm(n, mean, sd))
  if(!allownegatives) { r[r < 0] <- 0 }
  rts <- ts(r, 
            end = c(now$year + 1900, # 1900 is 0
                    now$mon + 1), # January is 0
            frequency = 12)
  rts
}