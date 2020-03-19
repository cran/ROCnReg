plot.cROC <-
function(x, ask = TRUE, ...) {
	change.ROC.format <- function(p, ROC) {
		temp <- reshape(ROC, varying = paste("p", round(p, 3), sep = ""), sep = "",
		v.names = "ROC", timevar = "p", times = p, idvar = "comb", direction = "long")
		temp[order(temp$comb),]
	}
	plot.accuracy <- function(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, accuracy, accuracy.main, dots, ask, ci.fit) {
		if(ask)
			readline("Press return for next page....")
		if(ci.fit) {
			accuracy.ci <- paste(accuracy,c("ql","qh"),sep="")
		}  
		if (n.cat == 0) {
			plot(ROC[ , names.cont], ROC[ , accuracy], xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type = "l", main = accuracy.main)			
			if(ci.fit) {
				lines(ROC[ , names.cont], ROC[ , accuracy.ci[1]], lty=2)
				lines(ROC[ , names.cont], ROC[ , accuracy.ci[2]], lty=2)
			}
			if (accuracy == "AUC") abline(h = 0.5, col = "grey")
		} else {
			if (n.cat == 1) {
				for(i in 1:dim.exp.cat) {
					if(ci.fit) {
						plot(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy], xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type="l", main = paste0(accuracy.main, " \n ", paste(names.cat, "=", exp.cat[i,,drop = TRUE])))			  
						lines(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy.ci[1]], lty=2)
						lines(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy.ci[2]], lty=2)
						if (accuracy == "AUC") abline(h = 0.5, col = "grey")
						if(ask & i < dim.exp.cat)
							readline("Press return for next page....")
					} else {
						if(i==1)
							plot(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy], xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type="l", main = accuracy.main)
						else 
							lines(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy], lty=i)
						if (accuracy == "AUC") abline(h = 0.5, col = "grey")
					}
				}
				if(!ci.fit)
					legend(if(!is.null(dots$pos.legend)) dots$pos.legend else "bottomright", legend = paste(names.cat, "=", exp.cat[, , drop = TRUE]),
					  lty=1:dim.exp.cat, cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1, y.intersp = if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 1)								  
			} else {
				for (i in 1:dim.exp.cat) {
					ind <- apply(t(apply(ROC[, names.cat], 1, function(x) x == exp.cat[i, ])), 1, all)				  
					if(ci.fit) {
						print(paste0(accuracy.main, " \n ", paste(paste(names.cat,"=",as.matrix(exp.cat)[i,]), collapse = ", ")))
						plot(ROC[ind, names.cont], ROC[ind, accuracy], xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type="l", main = paste0(accuracy.main, " \n ", paste(paste(names.cat,"=",as.matrix(exp.cat)[i,]), collapse = ", ")))
						lines(ROC[ind, names.cont], ROC[ind, accuracy.ci[1]], lty = i)
						lines(ROC[ind, names.cont], ROC[ind, accuracy.ci[2]], lty = i)					
						if (accuracy == "AUC") abline(h = 0.5, col = "grey")			
						if(ask & i < dim.exp.cat)
							readline("Press return for next page....")
					} else {
						if(i==1)
							plot(ROC[ind, names.cont], ROC[ind, accuracy], xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type="l", main = accuracy.main)
						else
							lines(ROC[ind, names.cont], ROC[ind, accuracy], lty = i)
						if (accuracy == "AUC") abline(h = 0.5, col = "grey")
					}
				}
				if(!ci.fit)
					legend(if(!is.null(dots$pos.legend)) dots$pos.legend else "bottomright", legend = sapply(1:dim.exp.cat, function(s) paste(paste(names.cat,"=",as.matrix(exp.cat)[s,]), collapse = ", ")),
					lty = 1:dim.exp.cat, cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1, y.intersp = if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 1)						  
			}
		}		   
	}

	dots <- list(...)
	ci.fit <- ifelse(is.null(x$ci.fit),FALSE,x$ci.fit)
	p <- x$p
	n.p <- length(p)
	colnames(x$ROC$est) <- paste("p", round(x$p, 3), sep = "")
	ROC <- cbind(x$newdata, x$ROC$est)	   
	set.accuracy <- c("AUC", "pAUC")
	if(is.null(x$pAUC)) {
		pAUC.legend <- ""
		pAUC.main <- ""
	} else {
		pAUC.legend <- ifelse(attr(x$pAUC, "focus") == "FPF", paste0(" (FPF = ", attr(x$pAUC, "value"), ")"), paste0(" (Se = ", attr(x$pAUC, "value"), ")"))
		pAUC.main <- ifelse(attr(x$pAUC, "focus") == "FPF", paste0("Partial area under the conditional ROC curve", pAUC.legend), paste0("Partial area under the specificity conditional ROC curve", pAUC.legend))
	}
	set.accuracy.legend <- c("AUC", paste0("pAUC", pAUC.legend))
	set.accuracy.main <- c("Area under the conditional ROC curve", pAUC.main) 
	ind.accuracy <- is.element(set.accuracy, set.accuracy[is.element(set.accuracy, names(x))])   
	if(any(ind.accuracy)) { 
		accuracy <- set.accuracy[ind.accuracy]
		accuracy.main <- set.accuracy.main[ind.accuracy]
		accuracy.legend <- set.accuracy.legend[ind.accuracy] 		  
	} else {
		accuracy <- NULL
	}   
	if(!is.null(accuracy)) {
		for (i in 1:length(accuracy)){
			aux <- names(ROC)		   
			ROC <- cbind(ROC, x[[accuracy[i]]])
			names(ROC) <- c(aux, colnames(x[[accuracy[i]]])) 
		}
	}
			
	range.marker <- range(x$data[, x$marker], na.rm = TRUE)	   
	
	names.cov <- names(ROC[, 1:(ncol(ROC) - n.p - (2*ci.fit+1)*sum(ind.accuracy)), drop = FALSE])
	ind.cat <- unlist(lapply(ROC[ ,names.cov, drop = FALSE], is.factor))
	names(ind.cat) <- names.cov	 
	names.cont <- names.cov[!ind.cat]		
	n.cont <- length(names.cont)
	n.cat <- sum(ind.cat)
	names.cat <- if(n.cat > 0) names.cov[ind.cat]
	
	delete.obs <- duplicated(x$newdata)
	ROC <- ROC[!delete.obs,,drop = FALSE]		 
	if (n.cont > 1) {
		#ROC <- ROC[order(ROC[, names.cont]),,drop = FALSE]
		ROC[, names.cont] <- apply(round(ROC[ , names.cont, drop = FALSE], 3), 2, factor)
	}
	if (n.cat > 0) {
		exp.cat <- unique(ROC[, names.cat, drop = FALSE])
		exp.cat.matrix <- as.matrix(exp.cat)
		dim.exp.cat <- nrow(exp.cat)
		levels.cat <- if(n.cat > 0) lapply(ROC[, names.cat, drop = FALSE], levels)			  
		n.levels <- as.numeric(unlist(lapply(levels.cat, length)))		  
		if(n.cont == 0) {
			ROC.long <- change.ROC.format(p, ROC)
			print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cat, collapse = "+"))),
			data = ROC.long,
			ylim = c(-0.1,1.05),
			xlab = "FPF",
			ylab = "TPF",
			strip = strip.custom(strip.names = TRUE, strip.levels = TRUE, sep = " = ",
			par.strip.text = list(cex = if(!is.null(dots$cex.par.strip.text)) dots$cex.par.strip.text else 0.75)),
			panel = function(x, y, subscripts) {
				panel.xyplot(x, y, type = "l")
				for (i in 1:length(set.accuracy))
					if(ind.accuracy[i]) { 
						acc.val <- round(unique(ROC.long[subscripts, set.accuracy[i]]),2)
						if(ci.fit) {
							acc.val <- paste(acc.val, "(", round(unique(ROC.long[subscripts, paste(set.accuracy[i],"ll",sep="")]),2),", ",round(unique(ROC.long[subscripts, paste(set.accuracy[i],"ul",sep="")]),2),")", sep="")
						}
						ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
						labels = paste(accuracy.legend[i],"=",acc.val), adj = c(1,0.5),
						cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
					}
			}))
		} else {
			if(n.cont == 1) {
				if (length(names.cat) == 1) {
					for(i in 1:dim.exp.cat) {
						if(i > 1) {
							if(ask)
								readline("Press return for next page....")
						}
						persp(p, ROC[ROC[ , names.cat] == exp.cat[i, ], names.cont], 
						t(as.matrix(ROC[ROC[ , names.cat] == exp.cat[i, ], -(c(1:(1 + n.cat), if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - (2*ci.fit+1)*length(accuracy))))])),
						xlab = "FPF", ylab = names.cont, zlab = "TPF",		  
						main = paste0("Conditional ROC surface \n ", exp.cat.matrix[i, ]),		  
						theta = if (!is.null(dots$theta))dots$theta else 20,
						phi = if (!is.null(dots$phi))dots$phi else 30,
						col = if(!is.null(dots$col))dots$col else "white",
						shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
						cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main, cex = dots$cex)							   
					}
					if(any(ind.accuracy))
						for(i in (1:length(set.accuracy))[ind.accuracy])
							plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], set.accuracy.main[i], dots, ask, ci.fit)					  
				} else {
					oldpar <- par(no.readonly = TRUE)
					on.exit(par(oldpar))
					par(mfrow = n.levels[1:2])   
					for (i in 1:(dim.exp.cat/prod(n.levels[1:2]))) {
						if(i > 1) {
							if(ask)
								readline("Press return for next page....")
						}
						k <- 0
						for (j in 1:(n.levels[1]*n.levels[2])) { 
							ind <- apply(t(apply(ROC[,names.cat], 1, function(x) x == exp.cat[(i-1)*prod(n.levels[1:2]) + j,])), 1, all)			
							persp(p, ROC[ind, names.cont],
							t(as.matrix(ROC[ind, -(c(1:(n.cont + n.cat),if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - (2*ci.fit+1)*length(accuracy))))])),
							xlab = "FPF", ylab = names.cont, zlab="TPF",
							main = paste(paste(names.cat, "=", c(exp.cat.matrix[j,1:2], exp.cat.matrix[1+(i-1)*prod(n.levels[1:2]),-(1:2)])), collapse = ", "),
							theta = if (!is.null(dots$theta))dots$theta else 20,
							phi = if (!is.null(dots$phi))dots$phi else 30,
							col = if(!is.null(dots$col))dots$col else "white",
							shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
							cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main,cex = dots$cex)						
						}
					}
					if(any(ind.accuracy))
						for(i in (1:length(set.accuracy))[ind.accuracy])
							plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], set.accuracy[i], dots, ask, ci.fit)
				}
			} else {				
				cat.cont <- vector("list", dim.exp.cat)					 
				for(i in 1:dim.exp.cat) {
					cat.cont[[i]] <- vector("list", n.cont)
					for(j in 1:n.cont) {
						ind <- t(apply(ROC[ , names.cat, drop = F], 1, function(x) x == exp.cat[i, ]))				  
						if(dim(ind)[1] == 1) ind <- t(ind)			  
						cat.cont[[i]][[j]] <- unique(ROC[apply(ind, 1, all), names.cont[j]])
					}
				}		 
				n.comb <- c(0, as.numeric(cumsum(unlist(lapply(cat.cont, function(x)cumprod(lapply(x, length))[n.cont]))))*n.p)			 
				ROC.long <- change.ROC.format(p, ROC)
				for (i in 1:dim.exp.cat) {
					if(i > 1 && ask) {
						readline("Press return for next page....")				  
					}
					print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
					data = ROC.long,
					ylim = c(-0.1,1.05),
					subset = (1 + n.comb[i]):n.comb[i + 1],
					strip = strip.custom(style = 3, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
					par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
					panel = function(x, y, subscripts) {
					panel.xyplot(x, y, type = "l")					  
						for (j in 1:length(set.accuracy))
							if(ind.accuracy[j]) { 
								acc.val <- round(unique(ROC.long[subscripts, set.accuracy[j]]),2)
								if(ci.fit) {
									acc.val <- paste(acc.val, "(", round(unique(ROC.long[subscripts, paste(set.accuracy[j],"ql",sep="")]),2),", ",round(unique(ROC.long[subscripts, paste(set.accuracy[j],"qh",sep="")]),2),")", sep="")
								}
								ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (j-1),
								labels = paste(accuracy.legend[j],"=",acc.val), adj = c(1,0.5),
								cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1)
							}
					},
					main = paste(names.cat, "=", exp.cat.matrix[i, ])))
				}
			}
		}
	} else {
		if(n.cont == 1) {   
			persp(p, ROC[ , names.cont], t(as.matrix(ROC[ ,-(c(1:(1 + n.cat), if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - (2*ci.fit+1)*length(accuracy))))])),
			xlab = "FPF", ylab = names.cont, zlab = "TPF",
			main = "Conditional ROC surface",	 
			theta = if (!is.null(dots$theta))dots$theta else 20,
			phi = if (!is.null(dots$phi))dots$phi else 30,
			col = if(!is.null(dots$col))dots$col else "white",
			shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
			cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main,cex = dots$cex)
			if(any(ind.accuracy))
				for(i in (1:length(set.accuracy))[ind.accuracy])
					plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], set.accuracy.main[i], dots, ask, ci.fit)
	   } else {
			ROC.long <- change.ROC.format(p, ROC)
			print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
			data = ROC.long,
			ylim = c(-0.1,1.05),
			strip = strip.custom(style = 3, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
			par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
			panel = function(x, y, subscripts) {
				panel.xyplot(x, y, type = "l")
				for (i in 1:length(set.accuracy))
					if(ind.accuracy[i]) { 
						acc.val <- round(unique(ROC.long[subscripts, set.accuracy[i]]),2)
						if(ci.fit) {
							acc.val <- paste(acc.val, "(", round(unique(ROC.long[subscripts, paste(set.accuracy[i],"ll",sep="")]),2),", ",round(unique(ROC.long[subscripts, paste(set.accuracy[i],"ul",sep="")]),2),")", sep="")
						}
						ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
						labels = paste(accuracy.legend[i],"=",acc.val), adj = c(1,0.5),
						cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
					}
			 }
			))
		}
	}
}
