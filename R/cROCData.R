cROCData <-
function(data, names.cov, group) {
	if(length(names.cov)==0) {	
		return(data.frame())
	}   
	card.cont <- 30		 
	ind.cat <- unlist(lapply(data[ , names.cov, drop = FALSE], is.factor))
	names(ind.cat) <- names.cov 
	names.cont <- names.cov[!ind.cat]
	n.cont <- length(names.cont)
	n.cat <- sum(ind.cat)
	names.cat <- if(n.cat > 0) names.cov[ind.cat]	
	if(length(card.cont) == 1 & n.cont > 1) card.cont <- rep(card.cont, n.cont) 
	levels.cat <- if(n.cat > 0) lapply(data[, names.cat, drop = FALSE], levels) 
	
	if (n.cat > 0) {
		exp.cat <- expand.grid(levels.cat)
		if(n.cont == 0) return(exp.cat)
		dim.exp.cat <- nrow(exp.cat)
		vector.indicators <- rep(TRUE, dim.exp.cat)

		calculate.ind.comb <- function(cat, comb) {
			ind <- t(apply(data[ , cat, drop = FALSE], 1, function(x) x == comb))
			if(dim(ind)[1] == 1) ind <- t(ind)
			apply(ind, 1, all)
		}

		ind.comb <- list()
		for(i in 1: dim.exp.cat)
			ind.comb[[i]] <- calculate.ind.comb(names.cat, exp.cat[i, ])		
		
		calculate.set.cont <- function(cont, ind, card) {
			indicator <- TRUE			   
			if (all(!ind)) {
				warning("Observations missing for at least one combination of levels of categorical covariates") 
				indicator <- FALSE 
			} else {
				range <- c(min(data[ind,cont], na.rm = TRUE), max(data[ind,cont], na.rm = TRUE))
				if (range[1] == range[2]) {
					warning(paste("There must be at least two different values of", cont, "at each combination of levels of categorical covariates"))
					indicator <- FALSE
				}
			}
			list(indicator = indicator, result = if (indicator) seq(range[1], range[2], by = (range[2] - range[1])/(card - 1)))
		}

		set.cont <- comb.cat.cont <- list.indicators <- vector("list", dim.exp.cat)
		for(i in 1:dim.exp.cat) {
			list.indicators[[i]] <- vector(length = n.cont)
			comb.cat.cont[[i]] <- vector("list", n.cont)
			names(comb.cat.cont[[i]]) <- names.cont
			for(j in 1:n.cont) {
				calc.set.cont <- calculate.set.cont(names.cont[j], ind.comb[[i]], card.cont[j])
				list.indicators[[i]][j] <- calc.set.cont$indicator
				comb.cat.cont[[i]][[j]] <- calc.set.cont$result
			}
			if(all(list.indicators[[i]]))  set.cont[[i]] <- expand.grid(comb.cat.cont[[i]])
		}
		vector.indicators <- unlist(lapply(list.indicators, all))
			
		data.ROC <- data.frame()
		for(i in 1:dim.exp.cat)
			data.ROC <- rbind(data.ROC, set.cont[[i]])
		for(i in 1:n.cat)
				data.ROC[,names.cat[i]] <- rep(exp.cat[vector.indicators, names.cat[i]], rep(cumprod(card.cont)[n.cont], sum(vector.indicators)))
   } else {		 
		cont <- vector("list", n.cont)
		names(cont) <- names.cont
		for(i in 1:n.cont) {				
			range <- c(max(tapply(data[,names(cont)[i]],data[,group], min, na.rm = TRUE), na.rm = TRUE), min(tapply(data[,names(cont)[i]],data[,group],max, na.rm = TRUE), na.rm = TRUE))			 
			cont[[i]] <- seq(range[1], range[2], by = (range[2] - range[1])/(card.cont[i] - 1))
		}
		data.ROC <- expand.grid(cont)   
	}   
	data.ROC
}
