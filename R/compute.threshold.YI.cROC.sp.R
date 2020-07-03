compute.threshold.YI.cROC.sp <-
function(object, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
  if(class(object)[1] != "cROC.sp") {
    stop(paste0("This function can not be used for this object class: ", class(object)[1]))
  }
  names.cov.h <- all.vars(object$formula$h)[-1]
  names.cov.d <- all.vars(object$formula$d)[-1]
  names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])
  
  if(!missing(newdata) && !inherits(newdata, "data.frame"))
    stop("Newdata must be a data frame")
  if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
    stop("Not all needed variables are supplied in newdata") 

  if(missing(newdata)) {
    newdata <- cROCData(object$data, names.cov, object$group)
  } else {
    newdata <- na.omit(newdata[,names.cov,drop = FALSE])
  }  

  # Compute F_D|X and F_{\bar{D}}|X
  pred0 <- predict(object$fit$h, newdata = newdata)
  pred1 <- predict(object$fit$d, newdata = newdata)
  sigma0 <- summary(object$fit$h)$sigma
  sigma1 <- summary(object$fit$d)$sigma
  
  npred <- nrow(newdata)
  
  y0 <- (object$data[object$data[,object$group] == object$tag.h,])[!object$missing.ind$h,object$marker]
  y1 <- (object$data[object$data[,object$group] != object$tag.h,])[!object$missing.ind$d,object$marker]

  n0 <- length(y0)
  n1 <- length(y1)

  grid  <- seq(min(c(y0, y1), na.rm = TRUE) - 1, max(c(y0, y1), na.rm = TRUE) + 1, length = max(500, c(n0,n1)))
  
  #grid  <- seq(min(object$data[, object$marker], na.rm = TRUE) - 1, max(object$data[, object$marker], na.rm = TRUE) + 1, length = 500)
  ngrid <- length(grid)

  F0 <- F1 <- matrix(0, nrow = ngrid, ncol = npred)
  thresholds.s <- YI.s <- TPF.s <- FPF.s <- vector(length = npred)

  for(l in 1:npred) {
    if(object$est.cdf == "normal") {
      F0[,l] <- pnorm(grid, mean = pred0[l], sd = sigma0)
      F1[,l] <- pnorm(grid, mean = pred1[l], sd = sigma1)

      difbb <- abs(F0[,l] -  F1[,l])
      thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])

      YI.s[l] <- max(difbb)
      TPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = pred1[l], sd = sigma1)
      FPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = pred0[l], sd = sigma0)

    } else {
      res0p <- object$fit$h$residuals/sigma0
      res1p <- object$fit$d$residuals/sigma1

      F0[,l] <- ecdf(res0p)((grid - pred0[l])/sigma0)
      F1[,l] <- ecdf(res1p)((grid - pred1[l])/sigma1)

      difbb <- abs(F0[,l] -  F1[,l])
      thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])

      YI.s[l] <- max(difbb)
      TPF.s[l] <- 1 - ecdf(res1p)((thresholds.s[l] - pred1[l])/sigma1)
      FPF.s[l] <- 1 - ecdf(res0p)((thresholds.s[l] - pred0[l])/sigma0)
    }
  }
  res <- list()
  res$call <- match.call()
  res$newdata <- newdata
  res$thresholds <- thresholds.s
  res$YI  <- YI.s
  res$FPF <- FPF.s
  res$TPF <- TPF.s
  res
  
}
