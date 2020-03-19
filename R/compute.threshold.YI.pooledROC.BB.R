compute.threshold.YI.pooledROC.BB <-
function(object) {
  if(class(object)[2] != "pooledROC.BB") {
    stop(paste0("This function can not be used for this object class: ", class(object)[2]))
  }
  B <- ncol(object$weights$h)

  weights.h <- object$weights$h
  weights.d <- object$weights$d

  y0 <- object$marker$h[!object$missing.ind$h]
  y1 <- object$marker$d[!object$missing.ind$d]

  thresholds.s <- YI.s <- FPF.s <- TPF.s <- numeric(B)

  #grid <- seq(min(c(y0,y1))-1, max(c(y0,y1))+1, length = 500)
  #ngrid <- length(grid)

  grid <- sort(unique(c(y0, y1)))
  ngrid <- length(grid)

  for(l in 1:B) {
    #F0bb <- apply(outer(y0, grid, "<="), 2, weighted.mean, w = weights.h[,l])
    #F1bb <- apply(outer(y1, grid, "<="), 2, weighted.mean, w = weights.d[,l])

    #difbb <- F0bb - F1bb
    #thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])  
    #YI.s[l] <- max(difbb)
    #TPF.s[l] <- 1 - apply(outer(y1, thresholds.s[l], "<="), 2, weighted.mean, w = weights.d[,l])
    #FPF.s[l] <- 1 - apply(outer(y0, thresholds.s[l], "<="), 2, weighted.mean, w = weights.h[,l])

    F0bb <- ewcdf(y0, weights.h[,l])
    F1bb <- ewcdf(y1, weights.d[,l])

    difbb <- F0bb(grid) - F1bb(grid)

    thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])  
    YI.s[l] <- max(difbb)
    TPF.s[l] <- 1 - F1bb(thresholds.s[l])
    FPF.s[l] <- 1 - F0bb(thresholds.s[l])

  }

  thresholds <- c(mean(thresholds.s), quantile(thresholds.s, c(0.025,0.975)))
  YI <- c(mean(YI.s), quantile(YI.s, c(0.025,0.975)))
  FPF <- c(mean(FPF.s), quantile(FPF.s, c(0.025,0.975)))
  TPF <- c(mean(TPF.s), quantile(TPF.s, c(0.025,0.975)))

  names(thresholds) <- names(YI) <- names(FPF) <- names(TPF) <- c("est","ql", "qh")

  res <- list()
  res$call <- match.call()
  res$thresholds <- thresholds
  res$YI <- YI
  res$FPF <- FPF
  res$TPF <- TPF
  res
}
