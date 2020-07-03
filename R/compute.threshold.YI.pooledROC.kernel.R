compute.threshold.YI.pooledROC.kernel <-
function(object, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
  if(class(object)[1] != "pooledROC.kernel") {
    stop(paste0("This function can not be used for this object class: ", class(object)[1]))
  }

  y0 <- object$marker$h[!object$missing.ind$h]
  y1 <- object$marker$d[!object$missing.ind$d]

  grid <- seq(min(c(y0,y1))-1, max(c(y0,y1))+1, length = max(c(length(y0), length(y1)), 500))
  ngrid <- length(grid)

  #grid <- sort(unique(c(y0, y1)))
  #ngrid <- length(grid)

  F0bb <- F1bb <- numeric(ngrid)
  for(i in 1:ngrid) {
    F0bb[i] <- Gk(grid[i], y = y0, h = object$bws$h)
    F1bb[i] <- Gk(grid[i], y = y1, h = object$bws$d)
  }

  difbb <- F0bb - F1bb

  thresholds <- mean(grid[which(difbb == max(difbb))])  
  YI <- max(difbb)
  FPF <- 1 - Gk(thresholds, y = y0, h = object$bws$h)
  TPF <- 1 - Gk(thresholds, y = y1, h = object$bws$d)

  res <- list()
  res$call <- match.call()
  res$thresholds <- thresholds
  res$YI <- YI
  res$FPF <- FPF
  res$TPF <- TPF
  res
}
