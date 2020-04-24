AROC.kernel <-
function(marker, covariate, group, tag.h, bw = c("LS","AIC"), regtype = c("LC","LL"), pauc = pauccontrol(), data, p = seq(0,1,l = 101), B = 1000) {
    compute.ROC <- function(marker, covariate, group, tag.h, bw, regtype, pauc, data, p = seq(0,1,l = 101)) {
        data.h <- data[data[,group] == tag.h,]
        data.d <- data[data[,group] != tag.h,]
        
        n0 <- nrow(data.h)
        n1 <- nrow(data.d)
        
        np <- length(p)
        
        x0 <- data.h[,covariate]
        y0 <- data.h[,marker]
        x1 <- data.d[,covariate]
        y1 <- data.d[,marker]
        
        # Fit the location-scale model in the healthy population
        bw.mean.h.p <-  npregbw(ydat = y0, xdat = x0, regtype = regtype, bwmethod = bw)
        fit.mean.h.p <- npreg(bw.mean.h.p, exdat = x0,residuals = TRUE)
        bw.var.h.p <- npregbw(ydat = (fit.mean.h.p$resid^2), xdat = x0, regtype = "lc", bwmethod = bw)
        fit.var.h.p <- npreg(bw.var.h.p, exdat = x0, residuals = TRUE)
        
        res0p <- fit.mean.h.p$resid/sqrt(fit.var.h.p$mean)
        F0res <- ecdf(res0p)
        
        # Evaluate the model in the diseased population, and compute the AROC
        fit.mean.d.p <- npreg(bw.mean.h.p, exdat = x1,residuals = TRUE)
        fit.var.d.p <- npreg(bw.var.h.p, exdat = x1, residuals = TRUE)
        u1 <- 1 - F0res((y1-fit.mean.d.p$mean)/sqrt(fit.var.d.p$mean))
        #arocp <- numeric(np)
        #for(i in 1:np){
        #    arocp[i] <- sum(u1<=p[i])/n1
        #}
        arocp <- apply(outer(u1, p, "<="), 2, mean)
        #aarocp <- simpson(arocp, p)
        aarocp <- mean(outer(res0p, (y1-fit.mean.d.p$mean)/sqrt(fit.var.d.p$mean), "<="))

        if(pauc$compute){
            if(pauc$focus == "FPF"){
                pu <- seq(0, pauc$value, len = np)
                #arocp_pauc <- numeric(np)
                #for(i in 1:np){
                #    arocp_pauc[i] <- sum(u1<=pu[i])/n1
                #}
                arocp_pauc <- apply(outer(u1, pu, "<="), 2, mean)
                #paarocp <- simpson(arocp_pauc, pu)
                paarocp <- pauc$value - mean(pmin(pauc$value, apply(outer(res0p, (y1-fit.mean.d.p$mean)/sqrt(fit.var.d.p$mean), ">="), 2, mean)))
            } else{
                arocp[1] <- 0
                arocp[np] <- 1
                rocapp <- approxfun(p, arocp, method = "linear")
                p1 <- uniroot(function(x) {rocapp(x) - pauc$value}, interval = c(0, 1))$root
                paarocp <- integrate(rocapp, lower = p1, upper = 1,
                stop.on.error = FALSE)$value - (1 - p1)*pauc$value
            }
        }
        res <- list()
        res$p <- p
        res$ROC <- arocp
        res$AUC <- aarocp
        if(pauc$compute){
            res$pAUC <- paarocp
        }
        res$data.h <- data.h
        res$data.d <- data.d
        res$fit <- list(bw.mean = bw.mean.h.p, bw.var = bw.var.h.p, fit.mean = fit.mean.h.p, fit.var = fit.var.h.p)
        res
    }
    
    pauc <- do.call("pauccontrol", pauc)
    
    np <- length(p)
    
    bw <- match.arg(bw)
    regtype <- match.arg(regtype)
    
    bw.aux <- switch(bw, "LS" = "cv.ls", "AIC" = "cv.aic")
    regtype.aux <- switch(regtype, "LC" = "lc", "LL" = "ll")
    
    # New data, removing missing values
    data.new <- data[,c(marker,group,covariate)]
    omit.h <- apply(data.new[data.new[,group] == tag.h, c(marker, group, covariate)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.h, c(marker, group, covariate)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.h,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.h,,drop = FALSE][!omit.d,,drop = FALSE])
    
    croc <- compute.ROC(marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, data = data.new, p = p)
    arocp <- croc$ROC
    aarocp <- croc$AUC
    if(pauc$compute){
        paarocp <- croc$pAUC
    }
    if(B > 0) {
        arocpb <- matrix(0, nrow = np, ncol = B)
        aarocpb <- numeric(B)
        if(pauc$compute){
            paarocpb <- numeric(B)
        }
        for(l in 1:B) {
            # Another option: healthy (residuals) - diseased (original sample)
            data.boot.d <- croc$data.d[sample(nrow(croc$data.d), replace=TRUE),]
            data.boot.h <- croc$data.h
            res.h.b <- sample(croc$fit$fit.mean$resid/sqrt(croc$fit$fit.var$mean), replace = TRUE)
            data.boot.h[,marker] <-croc$fit$fit.mean$mean + sqrt(croc$fit$fit.var$mean)*res.h.b
            data.boot <- rbind(data.boot.d, data.boot.h)
            
            res.boot <- compute.ROC(marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, data = data.boot, p = p)
            arocpb[,l] <- res.boot$ROC
            aarocpb[l] <- res.boot$AUC
            if(pauc$compute){
                paarocpb[l] <- res.boot$pAUC
            }
        }
    }
    columns <-switch(as.character(B>0),"TRUE" = 1:3,"FALSE"=1)
    col.names <-c("est","ql", "qh")[columns]
    
    poolROC <- matrix(0, ncol = length(columns), nrow = np, dimnames = list(1:np, col.names))
    poolROC[,1] <- arocp
    AUC <- aarocp
    if(pauc$compute){
        pAUC <- paarocp
    }
    if(B > 0) {
        poolROC[,2] <- apply(arocpb, 1, quantile, prob = 0.025)
        poolROC[,3] <- apply(arocpb, 1, quantile, prob = 0.975)
        AUC <- c(AUC, quantile(aarocpb, c(0.025, 0.975)))
        if(pauc$compute){
            pAUC <- c(pAUC, quantile(paarocpb, c(0.025, 0.975)))
        }
    }
    names(AUC) <- col.names
    if(pauc$compute){
        names(pAUC) <- col.names
    }
    res <- list()
    res$call <- match.call()
    res$data <- data
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$marker <- marker
    res$covariate <- covariate
    res$group <- group
    res$tag.h <- tag.h
    res$p <- p
    res$ROC <- poolROC
    res$AUC <- AUC
    if(pauc$compute){
        if(pauc$focus == "FPF"){
            res$pAUC <- pAUC/pauc$value
        } else{
            res$pAUC <- pAUC/(1 - pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    res$fit <- croc$fit
    class(res) <- c("AROC.kernel", "AROC")
    res
}
