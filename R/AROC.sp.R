AROC.sp <-
function(formula.h, group, tag.h, data, est.cdf.h = c("normal", "empirical"), pauc = pauccontrol(), p = seq(0,1,l = 101), B = 1000) {
    compute.AROC <- function(formula.h, data.h, data.d, est.cdf.h, pauc, p = seq(0,1,l = 101)) {
        np <- length(p)
        
        marker <- all.vars(formula.h)[1]
        
        # Fit the model in the healthy population
        fit0p <- lm(formula = formula.h, data = data.h)
        sigma0p <- summary(fit0p)$sigma
        pre.placement.values <- (data.d[,marker] - predict(fit0p, newdata = data.d))/sigma0p
        
        # Evaluate the model in the diseased population
        if(est.cdf.h == "normal") {
            u1 <- 1-pnorm(pre.placement.values)
        } else {
            res0p <- fit0p$residuals/sigma0p
            F0res <- ecdf(res0p)
            u1 <- 1 - F0res(pre.placement.values)
        }
        # Compute the AROC
        arocp <- apply(outer(u1, p, "<="), 2, mean)
        aarocp <- 1 - mean(u1)
        if(pauc$compute){
            if(pauc$focus == "FPF"){
                pu <- seq(0, pauc$value, len = np)
                arocp_pauc <- apply(outer(u1, pu, "<="), 2, mean)
                aarocp_pauc <- pauc$value - mean(pmin(pauc$value, u1))
            } else{
                arocp[1] <- 0
                arocp[np] <- 1
                rocapp <- approxfun(p, arocp, method = "linear")
                p1 <- uniroot(function(x) {rocapp(x) - pauc$value}, interval = c(0, 1))$root
                aarocp_pauc <- integrate(rocapp, lower = p1, upper = 1,
                stop.on.error = FALSE)$value - (1 - p1)*pauc$value
            }
        }
        
        res <- list()
        res$p <- p
        res$ROC <- arocp
        res$AUC <- aarocp
        if(pauc$compute){
            res$pAUC <- aarocp_pauc
        }
        res$fit <- fit0p
        res
    }
    est.cdf.h <- match.arg(est.cdf.h)
    np <- length(p)
    if(inherits(formula.h, "character")) {
        formula.h <- as.formula(formula.h)
    }
    # Marker variable
    tf <- terms.formula(formula.h)
    if (attr(tf, "response") > 0) {
        marker <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The formula should include the response variable (left hand side)")
    }
    
    # Variables in the model
    names.cov <- all.vars(formula.h)[-1]
    
    if(sum(is.na(match(c(marker, names.cov, group), names(data)))))
    stop("Not all needed variables are supplied in data")
    if(length(unique(data[,group])) != 2)
    stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep = "")
    
    pauc <- do.call("pauccontrol", pauc)
    
    # New data, removing missing values
    data.new <- data[,c(marker,group,names.cov)]
    omit.h <- apply(data.new[data.new[,group] == tag.h, c(marker, group, names.cov)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.h, c(marker, group, names.cov)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.h,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.h,,drop = FALSE][!omit.d,,drop = FALSE])
    data.new.h <- data.new[data.new[,group] == tag.h,]
    data.new.d <- data.new[data.new[,group] != tag.h,]

    res.fit <- compute.AROC(formula.h = formula.h, data.h = data.new.h, data.d = data.new.d, est.cdf.h = est.cdf.h, pauc = pauc, p = p)
    arocp  <- res.fit$ROC
    aarocp <- res.fit$AUC
    coeffp <- coefficients(res.fit$fit)

    if(pauc$compute){
        paarocp <- res.fit$pAUC
    }
    if(B > 0) {
        # Confidence intervals
        arocpb <- matrix(0, nrow = np, ncol = B)
        aarocpb <- numeric(B)
        coeffpb <- matrix(0, nrow = length(coeffp), ncol = B)
        if(pauc$compute){
            paarocpb <- numeric(B)
        }
        
        for(l in 1:B) {
            # Another option: healthy (residuals) - diseased (original sample)
            data.boot.d <- data.new.d[sample(nrow(data.new.d), replace=TRUE),]
            data.boot.h <- data.new.h
            res.h.b <- sample(res.fit$fit$residuals, replace = TRUE)
            data.boot.h[,marker] <- res.fit$fit$fitted + res.h.b
            
            res.boot <- compute.AROC(formula.h = formula.h, data.h = data.boot.h, data.d = data.boot.d, est.cdf.h = est.cdf.h, pauc = pauc, p = p)
            arocpb[,l] <- res.boot$ROC
            aarocpb[l] <- res.boot$AUC
            coeffpb[,l] <-  coefficients(res.boot$fit)

            if(pauc$compute){
                paarocpb[l] <- res.boot$pAUC
            }
        }
    }
    columns <-switch(as.character(B>0),"TRUE" = 1:3,"FALSE"=1)
    col.names <-c("est","ql", "qh")[columns]
    
    AROC <- matrix(0, ncol = length(columns), nrow = np, dimnames = list(1:np, col.names))
    AROC[,1] <- arocp

    AUC <- aarocp

    coeff <- matrix(0, ncol = length(columns), nrow = length(coeffp), dimnames = list(names(coeffp), col.names))
    coeff[,1] <- coeffp
    if(pauc$compute){
        pAUC <- paarocp
    }
    if(B > 0) {
        AROC[,2] <- apply(arocpb, 1, quantile, prob = 0.025)
        AROC[,3] <- apply(arocpb, 1, quantile, prob = 0.975)
        AUC <- c(AUC,quantile(aarocpb,c(0.025,0.975)))
        if(pauc$compute){
            pAUC <- c(pAUC,quantile(paarocpb,c(0.025,0.975)))
        }
        coeff[,2] <- apply(coeffpb, 1, quantile, prob = 0.025)
        coeff[,3] <- apply(coeffpb, 1, quantile, prob = 0.975)
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
    res$group <- group
    res$tag.h <- tag.h
    res$formula <- formula.h
    res$est.cdf.h <- est.cdf.h
    res$p <- p
    res$ROC <- AROC
    res$AUC <- AUC
     if(pauc$compute){
        if(pauc$focus == "FPF"){
            res$pAUC <- pAUC/pauc$value
        } else{
            res$pAUC <- pAUC/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    res$fit <- res.fit$fit
    res$coeff <- coeff
    class(res) <- c("AROC.sp", "AROC")
    res
}
