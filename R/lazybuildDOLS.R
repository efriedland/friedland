lazybuildDOLS <- function(results.coint, dat){
  mod.LagsnLeads <- mod.significance <- mod.adj.R2 <- mod.bic <- mod.Start <- mod.End <- mod.Nobs <- list()
  for(i in 1:length(levels(results.coint[,"Formula"]))){
    x <- levels(results.coint[,"Formula"])[i]
    dols <- buildDOLS(as.formula(x), data = dat)
    mod.LagsnLeads[[i]] <- dols$k
    mod.significance[[i]] <- tryCatch({
      coeftest(dols$model, vcov = NeweyWest(dols$model, lag = dols$k))[1:(length(unlist(strsplit(x,"[+]")))+1),"Pr(>|t|)"]
    }, error=function(e) return(NA))
    mod.adj.R2[[i]] <- summary(dols$model)$adj.r.squared
    mod.Start[[i]] <- as.yearmon(start(dols$model))[1]
    mod.End[[i]] <- as.yearmon(end(dols$model))[1]
    mod.Nobs[[i]] <- nrow(dols$model$model)
    mod.bic[[i]] <- bic(dols$model)
  }
  Yvar <- strsplit(as.character(results.coint[1,1]), " ~ ")[[1]][1]
  nVar <- colnames(dat)[-which(colnames(dat) == Yvar)]
  Compare_Cols <- c("(Intercept)", nVar, "adj.R2","Lags/Leads", "Start", "End", "Nobs", "bic", "AvgSignif")
  Compare_DOLS <- matrix(NA, nrow = length(mod.significance), ncol = length(Compare_Cols), dimnames = list(names(mod.significance), Compare_Cols))
  InCol <- lapply(mod.significance, function(x) match(names(x), Compare_Cols)) # where to put the significance values
  for(i in 1:length(mod.significance)){ Compare_DOLS[i, InCol[[i]]] <- round(mod.significance[[i]],3) }
  Compare_DOLS[,"adj.R2"] <- unlist(mod.adj.R2)
  Compare_DOLS[,"Lags/Leads"] <- unlist(mod.LagsnLeads)
  Compare_DOLS[,"Start"] <- unlist(mod.Start)
  Compare_DOLS[,"End"] <- unlist(mod.End)
  Compare_DOLS[,"Nobs"] <- unlist(mod.Nobs)
  Compare_DOLS[,"bic"] <- unlist(mod.bic)
  Compare_DOLS[,"AvgSignif"] <- rowSums(Compare_DOLS[,3:(1+length(nVar))], na.rm = T) /unlist(lapply(mod.significance, function(x) length(x)-1))
  rownames(Compare_DOLS) <- levels(results.coint[,"Formula"])
  as.data.frame(Compare_DOLS[order(Compare_DOLS[,"AvgSignif"]),]) # Comparison across all cointegrating buildDOLS results
}
