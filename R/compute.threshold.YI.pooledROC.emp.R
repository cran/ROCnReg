compute.threshold.YI.pooledROC.emp <-
function(object, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
  if(class(object)[1] != "pooledROC.emp") {
    stop(paste0("This function can not be used for this object class: ", class(object)[1]))
  }

  y0 <- object$marker$h[!object$missing.ind$h]
  y1 <- object$marker$d[!object$missing.ind$d]

  #grid <- seq(min(c(y0,y1))-1, max(c(y0,y1))+1, length = 500)
  #ngrid <- length(grid)

  grid <- sort(unique(c(y0, y1)))
  ngrid <- length(grid)

  F0bb <- ecdf(y0)(grid)
  F1bb <- ecdf(y1)(grid)

  difbb <- F0bb - F1bb
  thresholds <- mean(grid[which(difbb == max(difbb))])  
  YI <- max(difbb)
  TPF <- 1 - ecdf(y1)(thresholds)
  FPF <- 1 - ecdf(y0)(thresholds)

  res <- list()
  res$call <- match.call()
  res$thresholds <- thresholds
  res$YI <- YI
  res$FPF <- FPF
  res$TPF <- TPF
  res
}
