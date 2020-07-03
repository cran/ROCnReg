compute.threshold.YI.pooledROC.BB <-
function(object, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
  
  doMCMCTH <- function(k, y.h, weights.h, y.d, weights.d, grid) {
    F0bb <- ewcdf(y.h, weights.h[,k])
    F1bb <- ewcdf(y.d, weights.d[,k])

    difbb <- F0bb(grid) - F1bb(grid)

    thresholds.s <- mean(grid[which(difbb == max(difbb))])  
    YI.s <- max(difbb)
    TPF.s <- 1 - F1bb(thresholds.s)
    FPF.s <- 1 - F0bb(thresholds.s)

    res <- list()
    res$thresholds.s <- thresholds.s
    res$YI.s <- YI.s
    res$TPF.s <- TPF.s
    res$FPF.s <- FPF.s
    res
  }

  if(class(object)[1] != "pooledROC.BB") {
    stop(paste0("This function can not be used for this object class: ", class(object)[1]))
  }

  parallel <- match.arg(parallel)

  B <- ncol(object$weights$h)

  weights.h <- object$weights$h
  weights.d <- object$weights$d

  y.h <- object$marker$h[!object$missing.ind$h]
  y.d <- object$marker$d[!object$missing.ind$d]

  grid <- sort(unique(c(y.h, y.d)))

  if(B > 0) {        
        do_mc <- do_snow <- FALSE
        if (parallel != "no" && ncpus > 1L) {
            if (parallel == "multicore") {
                do_mc <- .Platform$OS.type != "windows"
            } else if (parallel == "snow") {
                do_snow <- TRUE
            }
            if (!do_mc && !do_snow) {
                ncpus <- 1L
            }       
            loadNamespace("parallel") # get this out of the way before recording seed
        }
        # Seed
        #if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
        #seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        # Apply function
        resBoot <- if (ncpus > 1L && (do_mc || do_snow)) {
                if (do_mc) {
                    parallel::mclapply(seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, grid = grid  , mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, grid = grid)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, grid = grid)
                        }                         
                    }
                }
            } else {
                lapply(seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, grid = grid)
            }

        resBoot <- simplify2array(resBoot)    
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        YI.s <- simplify2array(resBoot["YI.s",])
        TPF.s <- simplify2array(resBoot["TPF.s",])
        FPF.s <- simplify2array(resBoot["FPF.s",])
  } else {
    stop("B should be larger than zero.")
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
