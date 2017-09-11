
AICc<-function(n,k,LogLik){
  return(2*n*k/(n-k-1)+2*LogLik)
}

AkaikeWeight<-function(Delta.AICc.Array){
  return(exp(-Delta.AICc.Array/2) /sum(exp(-Delta.AICc.Array/2) ))
}

# se.function<-function(cov.matrix,var.name){
#   name.Index<-which(rownames(cov.matrix)==var.name)
#   if( length(name.Index)==1){
#     return( cov.matrix[name.Index,name.Index])
#   }else{return(0)}
# }
#
# var.model.Index.function<-
#   function(cov.matrix,var.name){
#     name.Index<-which(rownames(cov.matrix)==var.name)
#     if( length(name.Index)==1){
#       return( cov.matrix[name.Index,name.Index])
#     }else{return(0)}
#   }
#
# weight.para.value<-
#   function(para.vect,weight){
#     return(sum(para.vect*weight))
#   }
#

####################################################################
###################### MAIN PROGRAM ################################
####################################################################

#data is vector with names() matching taxon names
#phy is an ape phylo object
#flow is a data.frame with four columns
#donor = the taxon that is the gene flow source
#recipient = the taxon that is the gene flow recipient
#gamma = the fraction of the recipient trait that comes from the source. In the case of an equal hybridization between the recipient's sister on the tree and the donor, this is 0.5. In other cases where only, say, 10% of the recipient's quantitative trait
#	loci come from the donor, it would be 0.1
#time.from.root.donor = the time, counting forward FROM THE ROOT, when the gene flow happened from the donor. It may not be the same as time.from.root.recipient, as it may have spent time in a now extinct ghost lineage first (though time.from.root.donor <= time.from.root.recipient). It's treated as a one time event, which makes sense in the case of a single allopolyploid speciation event, probably less so in the case
#	of ongoing gene flow. Too bad.
#time.from.root.recipient = the time, counting forward FROM THE ROOT, when the gene flow happened from the donor
#If the gene flow happened to or from a lineage with multiple descendant species, use one row for each pair. For example, if lineage (A,B) had 20% of their genes coming in from lineage (C,D,E) at 14.5 MY since the root (not back in time), you would have
#	a flow data.frame of
#donor	recipient	gamma	time.from.root.donor	time.from.root.recipient
#C		A			0.2	14.5					14.5
#D		A			0.2	14.5					14.5
#E		A			0.2	14.5					14.5
#C		B			0.2	14.5					14.5
#D		B			0.2	14.5					14.5
#E		B			0.2	14.5					14.5
#We may write a utility function for dealing with this case in the future.
#Note the use of all updates of V.modified based on V.original; we don't want to add v_h to A three different times, for example, for one migration event (so we replace the variance three times based on transformations of the original variance)
#Note that we do not assume an ultrametric tree
BMhyb <- function(data, phy, flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=TRUE, get.se=TRUE, plot.se=TRUE, store.sims=FALSE, precision=2, auto.adjust=FALSE, likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000, measurement.error=0, do.kappa.check=FALSE, number.of.proportions=101, number.of.proportions.adaptive=101, allow.restart=TRUE, lower.bounds = c(0, -Inf, 0.000001, 0, 0), upper.bounds=c(10,Inf,100,100,100), badval.if.not.positive.definite=TRUE, attempt.deletion.fix=FALSE, starting.values=NULL, n.random.start.points=5000, do.Brissette.correction=FALSE, do.Higham.correction=TRUE, do.DE.correction=FALSE) {
  if(n.random.start.points>0 & is.null(starting.values)) {
    grid.results <- BMhybGrid(data=data, phy=phy, flow=flow, verbose=FALSE, precision=precision, n.points=n.random.start.points, attempt.deletion.fix=FALSE, measurement.error=measurement.error, get.se=FALSE, plot.se=FALSE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
    starting.values=grid.results$sims[which.min(grid.results$sims$AICc)[1],1:5]
    starting.names <- colnames(starting.values)
    starting.values <- as.numeric(starting.values)
    names(starting.values) <- starting.names
  }
	if(min(flow$gamma)<0) {
		stop("Min value of flow is too low; should be between zero and one")
	}
	if(max(flow$gamma)>1) {
		stop("Max value of flow is too high; should be between zero and one")
	}
	results<-list()
	#hessians <- list()
	results.summary <-data.frame()
	if(auto.adjust) {
		phy <- AdjustForDet(phy)
	}
	all.sims<-list()

  starting.from.geiger<-NA
  #starting.values <- NA
  if(is.null(starting.values)) {
    starting.values <- NA
    if(verbose) {
      print("Getting starting values from Geiger")
    }
    phy.geiger.friendly <- phy #geiger can't handle branch lengths near zero. Let's lengthen them if needed
    if(min(phy.geiger.friendly$edge.length)<0.00001) {
      phy.geiger.friendly$edge.length <- phy.geiger.friendly$edge.length + 0.00001
    }
    geiger.SE <- data*NA
    if(!is.null(measurement.error)) {
      if(length(measurement.error)==1) {
        geiger.SE <- rep(measurement.error, length(geiger.SE))
        names(geiger.SE) <- phy$tip.label
      } else {
        geiger.SE <- measurement.error
        names(geiger.SE) <- phy$tip.label

      }
      starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
      starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  0.01*starting.from.geiger$sigsq*max(branching.times(phy)), mean(measurement.error)) #sigma.sq, mu, beta, vh, SE

    } else {
  	    starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
  	    starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  0.01*starting.from.geiger$sigsq*max(branching.times(phy)), starting.from.geiger$SE) #sigma.sq, mu, beta, vh, SE
    }
  	if(verbose) {
  		print("Done getting starting values")
  	}
  }
  if(badval.if.not.positive.definite) {
    if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, actual.params= rep(TRUE,5)))) {
      if(attempt.deletion.fix) {
        phy <- AttemptDeletionFix(phy, flow, starting.values)
        tips <- tips[names(tips) %in% phy$tip.label]
      }
    #  stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
    }
  }
  for (model.index in sequence(length(models))) {
    do.run = TRUE
    preset.starting.parameters = NULL
    while(do.run) {
      do.run = FALSE
  		step.count <- 0
  		if(verbose) {
  			print(paste("Starting model", models[model.index], "of", length(models), "models"))
  		}
  		free.parameters<-rep(TRUE, 5)
  		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
  		model <- models[model.index]
  		if(model==1) {
  			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  		}
  		if(model==2) {
  			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  		}
  		if(model==3) {
  			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  		}
      if(!is.null(measurement.error)) {
        free.parameters[which(names(free.parameters)=="SE")]<-FALSE
      }

      if(is.null(preset.starting.parameters)) {
        preset.starting.parameters <- starting.values[free.parameters]
      }
      # if(badval.if.not.positive.definite) {
      #   if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, actual.params= free.parameters))) {
      #     stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
      #   }
      # }
  		best.run <- optim(par=preset.starting.parameters, fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
      best.run$par <- best.run$par
      attempts <- 1
  		while(best.run$convergence!=0){#want to get a convergence code 0
        if(verbose) {
          print(paste0("Initial search had a convergence code of ", best.run$convergence, ", indicating it did not converge. See ?optim for what the code may mean. Starting again, likely near that point. Negative log likelihood was ", best.run$value))
          print("Parameter estimates were")
          current.params <- best.run$par
          names(current.params) <- names(free.parameters)[which(free.parameters)]
          print(current.params)
        }
        attempts <- attempts+1
        if(attempts%%3!=1) {
    		    best.run<-optim(par=GenerateValues(best.run$par, lower=lower.bounds[which(free.parameters)], upper=upper.bounds[which(free.parameters)], examined.max=10*best.run$par, examined.min=0.1*best.run$par), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
        } else {
          best.run<-optim(par=GenerateRandomValues(data, free.parameters, lower.bounds[which(free.parameters)], upper.bounds[which(free.parameters)]), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
        }
          #best.run$par <- ConvertExpm1(best.run$par)

  		}


      # opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
      # best.run <- nloptr(x0=starting.values[free.parameters], eval_f=CalculateLikelihood, opts=opts, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], proportion.mix.with.diag=0, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions)
      # best.run$value = best.run$objective
      # best.run$par = best.run$solution
  		if(verbose) {
  			results.vector<-c(step.count, best.run$value, best.run$par)
  			names(results.vector) <- c("step","negloglik", names(free.parameters[which(free.parameters)]))
  			print(results.vector)
  		}
  		#this is to continue optimizing; we find that optim is too lenient about when it accepts convergence

      times.without.improvement <- 0
  		while(times.without.improvement<10) {
  			times.without.improvement <- times.without.improvement+1
        new.run <- optim(par=best.run$par, fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
        #new.run$par <- ConvertExpm1(new.run$par)
        attempts <- 1
        while(new.run$convergence!=0 & attempts < 20){#want to get a convergence code 0
          if(verbose) {
            print(paste0("This search had a convergence code of ", new.run$convergence, ", indicating it did not converge. See ?optim for what the code may mean. Starting again, likely near that point. Negative log likelihood was ", new.run$value))
            print("Parameter estimates were")
            current.params <- new.run$par
            names(current.params) <- names(free.parameters)[which(free.parameters)]
            print(current.params)
          }
          attempts <- attempts+1
          if(attempts%%3!=1) {
              new.run<-optim(par=GenerateValues(best.run$par, lower=lower.bounds[which(free.parameters)], upper=upper.bounds[which(free.parameters)], examined.max=10*best.run$par, examined.min=0.1*best.run$par), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
          } else {
            new.run<-optim(par=GenerateRandomValues(data, free.parameters, lower.bounds[which(free.parameters)], upper.bounds[which(free.parameters)]),  fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
          }
          #  new.run$par <- ConvertExpm1(new.run$par)

        }

  			#print("new.run best.run")
  			#print(c(new.run$value, best.run$value))

        if(new.run$value<best.run$value) {
  				if(best.run$value - new.run$value > likelihood.precision) {
  					times.without.improvement <- 0
  					if(verbose) {
  						print("New improvement found, resetting step counter")
  					}
  				} else {
  					if(verbose) {
  						print("New improvement found, but slight; taking the best value, but not resetting the step counter")
  					}
  				}
  				best.run <- new.run
  			}
  			if(verbose) {
  				step.count <- step.count+1
  				results.vector<-c(step.count, times.without.improvement, best.run$value, best.run$par)
  				names(results.vector) <- c("step", "steps.without.improvement","negloglik", names(free.parameters[which(free.parameters)]))
  				print(results.vector)
  			}
  		}#end of times.without.improvement<10
  		results[[model.index]] <- best.run
  		#try(hessians[[model.index]] <- hessian(func=CalculateLikelihood, x=new.run$par, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)]))
  		results.vector.full <- c(NA, NA, 1, 0, 0)
  		names(results.vector.full) <- names(free.parameters)
  		names(best.run$par) <- names(free.parameters[which(free.parameters)])
  		for (i in sequence(length(best.run$par))) {
  			results.vector.full[which(names(results.vector.full)==names(best.run$par)[i])] <- best.run$par[i]
  		}
  		#print(hessians[[model.index]])
  		#try(print(solve(hessians[[model.index]])))
  		ci.vector<-rep(NA,10)
  		for(parameter in sequence(5)) {
  			names(ci.vector)[1+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"lower", sep=".")
  			names(ci.vector)[2+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"upper", sep=".")
  		}
      weird.result <- FALSE
      if(best.run$value>1e100) {
        weird.result = TRUE
        warning(paste("It seems your likelihood for this run", best.run$value, "is invalid. This probably reflects a problem with numerical optimization for your tree. Your parameter estimates and AIC weights are likely meaningless. We will not calculate confidence: the confidence intervals for all parameters should be taken as c(-Inf, Inf) [unless the parameter is bounded by zero, in which case it is c(0, Inf)]. You could try a transformation of your branch lengths (and make sure to change the time of the flow matrix, too) -- note that parameter values (rate of evolution) will be based on this new timescale. Alternatively, you might allow SE for tip values to be estimated by setting measurement.error=NULL."))
        print(paste("It seems your likelihood for this run", best.run$value, "is invalid. This probably reflects a problem with numerical optimization for your tree. Your parameter estimates and AIC weights are likely meaningless. We will not calculate confidence: the confidence intervals for all parameters should be taken as c(-Inf, Inf) [unless the parameter is bounded by zero, in which case it is c(0, Inf)]. You could try a transformation of your branch lengths (and make sure to change the time of the flow matrix, too) -- note that parameter values (rate of evolution) will be based on this new timescale. Alternatively, you might allow SE for tip values to be estimated by setting measurement.error=NULL."))
      }
  		if(get.se & !weird.result) {
  			if(verbose) {
  				print("Now doing simulation to estimate parameter uncertainty")
  			}
  			interval.results <- AdaptiveConfidenceIntervalSampling(best.run$par, fn=CalculateLikelihood, lower=lower.bounds[which(free.parameters)], upper=upper.bounds[which(free.parameters)], data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], allow.extrapolation=allow.extrapolation, n.points=n.points,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions.adaptive, allow.restart=allow.restart, best.lnl = best.run$value, likelihood.precision=likelihood.precision, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
  			interval.results.in <- interval.results[which(interval.results[,1]-min(interval.results[,1])<=2),]
  			interval.results.out <- interval.results[which(interval.results[,1]-min(interval.results[,1])>2),]
        if(best.run$value - min(interval.results[,1]) > likelihood.precision) {
          print("The sampling to find confidence in parameters actually found a better part of the likelihood surface. Restarting the run for this model at that point")
          best.point <- interval.results[which.min(interval.results[,1]),]
          names(best.point) <- c("neglnL", names(free.parameters)[which(free.parameters)])
          preset.starting.parameters <- best.point[-1]
          do.run = TRUE
        }
  			if(plot.se) {
  				pdf(file=paste("Model",models[model.index], "_uncertainty_plot.pdf", sep=""), height=5, width=5*sum(free.parameters))
  				par(mfcol=c(1, sum(free.parameters)))
  				for(parameter in sequence(sum(free.parameters))) {
  					plot(x=interval.results[,parameter+1], y=interval.results[,1], type="n", xlab=names(free.parameters[which(free.parameters)])[parameter], ylab="NegLnL", bty="n", ylim=c(min(interval.results[,1]), min(interval.results[,1])+10))
  					points(x=interval.results.in[,parameter+1], y=interval.results.in[,1], pch=16, col="black")
  					points(x=interval.results.out[,parameter+1], y=interval.results.out[,1], pch=16, col="gray")
  					points(x= best.run$par[parameter], y= best.run$value, pch=1, col="red", cex=1.5)
  				}
  				dev.off()
  				if(verbose) {
  					print(paste("Uncertainty plot has been saved in Model",models[model.index], "_uncertainty_plot.pdf in ", getwd(), sep=""))
  				}
  			}
  			if(store.sims) {
  				colnames(interval.results) <- c("neglnL", names(free.parameters)[which(free.parameters)])
  				all.sims[[model.index]]<-interval.results
  			}
  			free.index=0
  			for(parameter in sequence(5)) {

  				if(free.parameters[parameter]) { #is estimated
  					free.index <- free.index + 1
  					ci.vector[1+2*(parameter-1)] <- min(interval.results.in[,free.index+1])
  					ci.vector[2+2*(parameter-1)] <- max(interval.results.in[,free.index+1])
  				} else {
  					ci.vector[1+2*(parameter-1)] <- results.vector.full[parameter]
  					ci.vector[2+2*(parameter-1)] <- results.vector.full[parameter]
  				}
  			}
  		}
  		local.df <- data.frame(matrix(c(models[model.index], results.vector.full, AICc(Ntip(phy),k=length(free.parameters[which(free.parameters)]), best.run$value), best.run$value, length(free.parameters[which(free.parameters)]), ci.vector), nrow=1), stringsAsFactors=FALSE)
  		colnames(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K", names(ci.vector))
      if(do.Higham.correction) {
        local.df$penalty=CalculateLikelihood(best.run$par,data=data, phy=phy, flow=flow,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, actual.params=free.parameters[which(free.parameters)], number.of.proportions=number.of.proportions.adaptive,  likelihood.precision=likelihood.precision, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction, return.penalty=TRUE)
      }
  		print(local.df)
  		results.summary <- rbind(results.summary, local.df)
    }
	}
	results.summary <- cbind(results.summary, deltaAICc=results.summary$AICc-min(results.summary$AICc))
	results.summary<-cbind(results.summary, AkaikeWeight = AkaikeWeight(results.summary$deltaAICc))
	if(store.sims) {
		return(list(results=results.summary, sims=all.sims))
	}
	return(results.summary)
}

BMhybGrid <- function(data, phy, flow, models=c(1,2,3,4), verbose=TRUE, get.se=TRUE, plot.se=TRUE, store.sims=TRUE, precision=2, auto.adjust=FALSE, likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000, measurement.error=0, do.kappa.check=FALSE, number.of.proportions=101, number.of.proportions.adaptive=101, allow.restart=TRUE, lower.bounds = c(0, -Inf, 0.000001, 0, 0), upper.bounds=c(10,Inf,100,100,100), badval.if.not.positive.definite=TRUE, attempt.deletion.fix=FALSE, starting.values=NULL, do.Brissette.correction=FALSE, do.Higham.correction=TRUE, do.DE.correction=FALSE) {
	if(min(flow$gamma)<0) {
		stop("Min value of flow is too low; should be between zero and one")
	}
	if(max(flow$gamma)>1) {
		stop("Max value of flow is too high; should be between zero and one")
	}
	results<-list()
	#hessians <- list()
	results.summary <-data.frame()
	phy.geiger.friendly <- phy #geiger can't handle branch lengths near zero. Let's lengthen them if needed
	if(min(phy.geiger.friendly$edge.length)<0.00001) {
		phy.geiger.friendly$edge.length <- phy.geiger.friendly$edge.length + 0.00001
	}
	if(auto.adjust) {
		phy <- AdjustForDet(phy)
	}
	all.sims<-list()
  if(is.null(starting.values)) {
    if(verbose) {
  		print("Getting starting values from Geiger")
  	}
    starting.from.geiger<-NA
    starting.values <- NA
    geiger.SE <- data*NA
    if(!is.null(measurement.error)) {
      if(length(measurement.error)==1) {
        geiger.SE <- rep(measurement.error, length(geiger.SE))
        names(geiger.SE) <- phy$tip.label
      } else {
        geiger.SE <- measurement.error
        names(geiger.SE) <- phy$tip.label

      }
      starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
      starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  starting.from.geiger$sigsq*max(branching.times(phy)), mean(measurement.error)) #sigma.sq, mu, beta, vh, SE

    } else {
  	    starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
  	    starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  starting.from.geiger$sigsq*max(branching.times(phy)), starting.from.geiger$SE) #sigma.sq, mu, beta, vh, SE
    }
  	if(verbose) {
  		print("Done getting starting values")
  	}
  }
#  if(badval.if.not.positive.definite) {
    if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, actual.params= rep(TRUE,5)))) {
      if(attempt.deletion.fix) {
        phy <- AttemptDeletionFix(phy, flow, starting.values)
        tips <- tips[names(tips) %in% phy$tip.label]
      }
  #    stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
    }
#  }
  for (model.index in sequence(length(models))) {
    do.run = TRUE
    preset.starting.parameters = NULL
    while(do.run) {
      do.run = FALSE


  		step.count <- 0
  		if(verbose) {
  			print(paste("Starting model", models[model.index], "of", length(models), "models"))
  		}
  		free.parameters<-rep(TRUE, 5)
  		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
  		model <- models[model.index]
  		if(model==1) {
  			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  		}
  		if(model==2) {
  			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  		}
  		if(model==3) {
  			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  		}
      if(!is.null(measurement.error)) {
        free.parameters[which(names(free.parameters)=="SE")]<-FALSE
      }

      starting.mins <- c(0, min(data)-(max(data)-min(data)), 0, 0, 0)
      for(i in sequence(length(starting.mins))) {
        starting.mins[i] <- max(starting.mins[i], lower.bounds[i])
      }
      starting.maxes <- c(10*starting.values[1], max(data)+(max(data)-min(data)), 10*starting.values[3], 10*starting.values[4], 10*starting.values[5])
      for (i in sequence(length(starting.maxes))) {
        starting.maxes[i] <- min(starting.maxes[i], upper.bounds[i])
      }
      names(starting.mins) <- names(free.parameters)
      names(starting.maxes) <- names(free.parameters)
      ci.vector<-rep(NA,10)
      for(parameter in sequence(5)) {
        names(ci.vector)[1+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"lower", sep=".")
        names(ci.vector)[2+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"upper", sep=".")
      }

      # if(badval.if.not.positive.definite) {
      #   if(!IsPositiveDefinite(GetVModified(preset.starting.parameters, phy, flow, actual.params= free.parameters))) {
      #     stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
      #   }
      # }

      if(model==1) {
        starting.mins["bt"] <- 1
        starting.maxes["bt"] <- 1
      }
      if(model==2) {
        starting.mins["vh"] <- 0
        starting.maxes["vh"] <- 0
      }
      if(model==3) {
        starting.mins["bt"] <- 1
        starting.maxes["bt"] <- 1
        starting.mins["vh"] <- 0
        starting.maxes["vh"] <- 0
  		}
      if(!is.null(measurement.error)) {
        starting.mins["SE"] <- 0
        starting.maxes["SE"] <- 0
      }

      #starting.mins <- starting.mins[free.parameters]
    #  starting.maxes <- starting.mins[free.parameters]
      grid.of.points <- lhs::randomLHS(n=n.points, k=length(starting.mins))
      for(parameter.index in sequence(ncol(grid.of.points))) {
        grid.of.points[,parameter.index] <- starting.mins[parameter.index] + grid.of.points[,parameter.index] * (starting.maxes[parameter.index] - starting.mins[parameter.index])
      }
      colnames(grid.of.points) <- names(free.parameters)
      likelihoods <- rep(NA, n.points)

      for (rep.index in sequence(n.points)) {
        local.likelihood <- try(CalculateLikelihood(as.numeric(grid.of.points[rep.index,]), data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], measurement.error=measurement.error, badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction))
        if(!is.numeric(local.likelihood)) {
          local.likelihood <- (0.5)*.Machine$double.xmax
        }
        likelihoods[rep.index] <- local.likelihood
        if(verbose & rep.index%%50==0) {
          print(paste0("Now done ", rep.index, " of ", n.points, " to analyze (", round(100*rep.index/n.points, 4), "% done)"))
        }
      }

      best.one <- which.min(likelihoods)[1]
      best.params <- grid.of.points[best.one,]
      best.likelihood <- min(likelihoods, na.rm=TRUE)


      results.vector.full <- c(NA, NA, 1, 0, 0)
      names(results.vector.full) <- names(free.parameters)
    #  names(best.run$par) <- names(free.parameters[which(free.parameters)])
      for (i in sequence(length(best.params))) {
        results.vector.full[which(names(results.vector.full)==names(best.params)[i])] <- best.params[i]
      }

    }

    if(get.se) {
      if(verbose) {
        print("Now doing simulation to estimate parameter uncertainty")
      }
      previous.results <- cbind(likelihoods, grid.of.points)
      colnames(previous.results)[1] <- "NegLogL"
      interval.results <- AdaptiveConfidenceIntervalSampling(best.params[free.parameters], fn=CalculateLikelihood, lower=lower.bounds[free.parameters], upper=upper.bounds[free.parameters], data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], allow.extrapolation=allow.extrapolation, n.points=n.points,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions.adaptive, allow.restart=allow.restart, best.lnl = min(likelihoods), likelihood.precision=likelihood.precision, lower.b=lower.bounds[free.parameters], upper.b=upper.bounds[free.parameters], do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
      colnames(interval.results) <- c("NegLogL", names(free.parameters[free.parameters]))
      if(!any(grepl("bt", names(interval.results)))) {
        interval.results$bt <- 1
      }
      if(!any(grepl("vh", names(interval.results)))) {
        interval.results$vh <- 0
      }
      if(!any(grepl("SE", names(interval.results)))) {
        interval.results$SE <- 0
      }
      interval.results <- interval.results[, colnames(previous.results)]
      interval.results.local <- interval.results
      interval.results <- rbind(previous.results, interval.results)
      interval.results <- interval.results[is.finite(interval.results[,1]),]
      interval.results.in <- interval.results[which(interval.results[,1]-min(as.numeric(interval.results[,1]))<=2),]
      interval.results.out <- interval.results[which(interval.results[,1]-min(as.numeric(interval.results[,1]))>2),]
      # if(best.run$value - min(interval.results[,1]) > likelihood.precision) {
      #   print("The sampling to find confidence in parameters actually found a better part of the likelihood surface. Restarting the run for this model at that point")
      #   best.point <- interval.results[which.min(interval.results[,1]),]
      #   names(best.point) <- c("neglnL", names(free.parameters)[which(free.parameters)])
      #   preset.starting.parameters <- best.point[-1]
      #   do.run = TRUE
      # }
      for(parameter in sequence(ncol(interval.results)-1)) {
        parameter.name <- names(interval.results)[parameter+1]
        ci.vector[paste0(parameter.name, ".upper")] <- max(interval.results.in[,parameter+1], na.rm=TRUE)
        ci.vector[paste0(parameter.name, ".lower")] <- min(interval.results.in[,parameter+1], na.rm=TRUE)
      }
      if(min(interval.results$NegLogL, na.rm=TRUE) <  min(likelihoods)) {
          best.params <- interval.results[which.min(interval.results$NegLogL), -1]
            results.vector.full <- c(NA, NA, 1, 0, 0)
            names(results.vector.full) <- names(free.parameters)
            for (i in sequence(length(best.params))) {
              results.vector.full[which(names(results.vector.full)==names(best.params)[i])] <- best.params[i]
            }
            best.likelihood <- min(interval.results$NegLogL, na.rm=TRUE)
      }
      if(plot.se) {
        pdf(file=paste("Model",models[model.index], "_uncertainty_plot.pdf", sep=""), height=5, width=5*sum(free.parameters))
        par(mfcol=c(1, sum(free.parameters)))
        for(parameter in sequence(sum(free.parameters))) {
          plot(x=interval.results[,parameter+1], y=interval.results[,1], type="n", xlab=names(free.parameters[which(free.parameters)])[parameter], ylab="NegLnL", bty="n", ylim=c(min(interval.results[,1]), min(interval.results[,1])+10))
          points(x=interval.results.in[,parameter+1], y=interval.results.in[,1], pch=16, col="black")
          points(x=interval.results.out[,parameter+1], y=interval.results.out[,1], pch=16, col="gray")
          points(x= best.params[parameter], y= best.likelihood, pch=1, col="red", cex=1.5)
        }
        dev.off()
        if(verbose) {
          print(paste("Uncertainty plot has been saved in Model",models[model.index], "_uncertainty_plot.pdf in ", getwd(), sep=""))
        }
      }
    }
    local.df <- data.frame(matrix(c(models[model.index], results.vector.full, AICc(Ntip(phy),k=length(free.parameters[which(free.parameters)]), best.likelihood), best.likelihood, length(free.parameters[which(free.parameters)]), ci.vector), nrow=1), stringsAsFactors=FALSE)
    local.df <- apply(local.df, 2, unlist)
    names(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K", names(ci.vector))

    # local.df <- data.frame(matrix(c(models[model.index], results.vector.full, AICc(Ntip(phy),k=length(free.parameters[which(free.parameters)]), likelihoods[best.one]), likelihoods[best.one], length(free.parameters[which(free.parameters)])), nrow=1))
    # colnames(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K")
    # print(local.df)
    if(do.Higham.correction) {
      param.estimates <- unlist(results.vector.full)
      names(param.estimates) <- names(free.parameters)
      local.df$penalty=CalculateLikelihood(param.estimates,data=data, phy=phy, flow=flow,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, actual.params=free.parameters[which(free.parameters)], number.of.proportions=number.of.proportions.adaptive,  likelihood.precision=likelihood.precision, lower.b=lower.bounds[which(free.parameters)], upper.b=upper.bounds[which(free.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction, return.penalty=TRUE)
    }
    all.points <- data.frame(grid.of.points, stringsAsFactors=FALSE)
    all.points$NegLogL <- likelihoods
    all.points$Model <- models[model.index]
    all.points$AICc <- AICc(Ntip(phy),k=length(free.parameters[which(free.parameters)]), all.points$NegLogL)
    all.points$K <- length(free.parameters[which(free.parameters)])
    results.summary <- rbind(results.summary, data.frame(t(local.df), stringsAsFactors=FALSE))
    all.sims <- rbind(all.sims, all.points)
    # if(plot.se) {
    #   pdf(file=paste("Model",models[model.index], "_uncertainty_plot.pdf", sep=""), height=5, width=5*sum(free.parameters))
    #   par(mfcol=c(1, sum(free.parameters)))
    #   for(parameter in sequence(length(free.parameters))) {
    #     if(free.parameters[parameter]) {
    #
    #       plot(x=all.points[,names(free.parameters)[parameter]], y=all.points[,"NegLogL"], type="n", xlab=names(free.parameters)[parameter], ylab="NegLnL", bty="n", ylim=c(min(all.points[,"NegLogL"]), min(all.points[,"NegLogL"])+10))
    #       points(x=all.points[,names(free.parameters)[parameter]], y=all.points[,"NegLogL"], pch=16, col=ifelse(all.points[,"NegLogL"] < (min(all.points[,"NegLogL"])+2), "black", "gray"))
    #       points(x= all.points[which.min(all.points[,"NegLogL"])[1],names(free.parameters)[parameter]], y= all.points$NegLogL[which.min(all.points[,"NegLogL"])[1]], pch=1, col="red", cex=1.5)
    #     }
    #   }
    #   dev.off()
    #   if(verbose) {
    #     print(paste("Uncertainty plot has been saved in Model",models[model.index], "_uncertainty_plot.pdf in ", getwd(), sep=""))
    #   }
    # }
	}
	results.summary <- cbind(results.summary, deltaAICc=as.numeric(results.summary$AICc)-min(as.numeric(results.summary$AICc)))
	results.summary<-cbind(results.summary, AkaikeWeight = AkaikeWeight(results.summary$deltaAICc))
  #save(list=ls(), file="~/Desktop/everything.rda")
	if(store.sims) {
    all.sims$deltaAICc <- all.sims$AICc - min(as.numeric(all.sims$AICc))
    all.sims$AkaikeWeight <- AkaikeWeight(all.sims$deltaAICc)
		return(list(results=results.summary, sims=all.sims))
	}
	return(results.summary)
}

PlotAICRegion <- function(sims, show.finite.only=TRUE, true.params=NULL, ...) {
  pairs.of.params <- utils::combn(colnames(sims)[1:5],2)
  sims.to.plot <- sims
  if(show.finite.only) {
    sims.to.plot <- sims.to.plot[which(sims.to.plot$NegLogL<1e300),]
  }
  best.one <- which.min(sims.to.plot$AICc)[1]
  par(mfcol=c(2, ceiling(ncol(pairs.of.params)/2)))
  for(pair.index in sequence(ncol(pairs.of.params))) {
    plot(sims.to.plot[,pairs.of.params[1,pair.index]], sims.to.plot[,pairs.of.params[2,pair.index]], pch=".", col=rgb(0,0,0,.5), xlim=range(sims[,pairs.of.params[1,pair.index]]), ylim=range(sims[,pairs.of.params[2,pair.index]]), xlab=paste0(pairs.of.params[1,pair.index],ifelse(0==max(sims[,pairs.of.params[1,pair.index]])-min(sims[,pairs.of.params[1,pair.index]]), " FIXED", "")), ylab=paste0(pairs.of.params[2,pair.index],ifelse(0==max(sims[,pairs.of.params[2,pair.index]])-min(sims[,pairs.of.params[2,pair.index]]), " FIXED", "")), ...)
    if(!is.null(true.params)) {
      points(x=true.params[pairs.of.params[1,pair.index]], y=true.params[pairs.of.params[2,pair.index]], pch=8, col="blue")
    }
    points(x=sims.to.plot[best.one,pairs.of.params[1,pair.index]], y=sims.to.plot[best.one,pairs.of.params[2,pair.index]], pch=5, col="red")
  }
}

# PlotUncertainty <- function(results, model.index, make.pdf=TRUE, region=2) {
#   model.sims <- results$sims[[model.index]]
#   best.likelihood <- min(model.sims$neglnL)
#   model.sims.in <- subset(model.sims, neglnL < best.likelihood+region)
#   model.sims.out <- subset(model.sims, neglnL >= best.likelihood+region)
#   best.run <- model.sims[which.min(model.sims$neglnL)[1],]
#   number.free.parameters <- dim(model.sims)[2] - 1
#   free.parameters <- colnames(model.sims)[-1]
#   if(make.pdf) {
#     pdf(file=paste("Model",model.index, "_uncertainty_plot.pdf", sep=""), height=5, width=5*number.free.parameters)
#   }
#   par(mfcol=c(1, number.free.parameters))
#   for(parameter in sequence(number.free.parameters)) {
#     plot(x=model.sims[,parameter+1], y=model.sims[,1], type="n", xlab=free.parameters[parameter], ylab="NegLnL", bty="n", ylim=c(min(model.sims[,1]), min(model.sims[,1])+10))
#     points(x=model.sims.in[,parameter+1], y=model.sims.in[,1], pch=16, col="black")
#     points(x=model.sims.out[,parameter+1], y=model.sims.out[,1], pch=16, col="gray")
#     points(x= best.run[parameter+1], y= best.run[1], pch=1, col="red", cex=1.5)
#   }
#   if(make.pdf) {
#     dev.off()
#   }
# }

DetPass <- function(phy) {
	det.pass <- TRUE
	vcv.result <- vcv.phylo(phy)
	det.tries <- c(det(vcv.result), det(1000*vcv.result), det(0.0001*vcv.result))
	if(min(det.tries)<0) {
		det.pass <- FALSE
	}
	if(sum(is.finite(det.tries))!=length(det.tries)) {
		det.pass <- FALSE
	}
	return(det.pass)
}

AdjustForDet <- function(phy, max.attempts=100) {
	attempts<-0
	if(!DetPass(phy)) {
		warning("Determininant of the phylogeny was difficulty to calculate, so the phylogeny needed to be adjusted. Your results may be approximate as a result")
		while(!DetPass(phy) && attempts<=max.attempts) {
			attempts <- attempts+1
			phy$edge.length <- phy$edge.length+0.00001*attempts
		}
	}
	if(attempts>max.attempts) {
		stop("Phylogeny could not be adjusted adequately")
	}
	return(phy)
}

IsPositiveDefinite <- function(V.modified) {
  eigenvalues <- eigen(V.modified)$values
  if(is.complex(eigenvalues)) {
    return(FALSE)
  } else {
    return(min(eigenvalues)>0)
  }
}

BrissetteEtAlCorrection <- function(V.modified, min.eigenvalue=1e-6, max.attempts=10) {
  V.corrected <- V.modified
#  if(!IsPositiveDefinite(V.modified)) {
#    correction.factor <- sqrt(diag(V.modified) %*% t(diag(V.modified)))
#    V.corrected <- V.modified / correction.factor
#  }
  attempt.count <- 0
  while(!IsPositiveDefinite(V.corrected) & attempt.count<max.attempts) {
    attempt.count <- attempt.count + 1
    V.eigen <- eigen(V.modified)
    V.eigen$values[which(V.eigen$values<=0)] <- min.eigenvalue
    V.corrected <- V.eigen$vectors %*% diag(V.eigen$values) %*% t(V.eigen$vectors)
    if(min(V.corrected)<0) {
      V.corrected <- V.modified # so that we start the loop again, but with a different min.eigenvalue
      min.eigenvalue <- min.eigenvalue * 1.5
    }
#    correction.factor <- sqrt(diag(V.corrected) %*% t(diag(V.corrected)))
#    V.corrected <- V.corrected / correction.factor
    warning(paste0("Needed to use Brissette et al. (2007) correction; maximum change was ", max(abs(V.corrected - V.modified)), " and maximum value in V.modified to start was ", max(V.modified)))
  }
  if(min(V.corrected) < 0) {
    V.corrected <- NULL
  }
  return(V.corrected)
}

GetVModified <- function(x, phy, flow, actual.params, measurement.error=NULL) {
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
  SE <- 0
  if(is.null(measurement.error)) {
	   SE <- x[length(x)]
}
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]
	}
	times.original <-vcv.phylo(phy, model="Brownian") #is the initial one based on the tree alone, so just time
	V.original <- sigma.sq * times.original
	V.modified <- V.original
	for (flow.index in sequence(dim(flow)[1])) {
		recipient.index <- which(rownames(V.modified)==flow$recipient[flow.index])
    gamma <- flow$gamma[flow.index]
		if(length(recipient.index)!=1) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but instead found ", paste(rownames(V.modified)[recipient.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		donor.index <- which(rownames(V.modified)==flow$donor[flow.index])
		if(length(donor.index)!=1) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but instead found ", paste(rownames(V.modified)[donor.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", sep=""))
		}
		V.modified[recipient.index, donor.index] <- (1-gamma) * V.original[recipient.index, donor.index] + (gamma) * (flow$time.from.root.recipient[flow.index]) * sigma.sq #covariance is the weighted sum of the covariance from evolution along the tree plus evolution along the migration path
		V.modified[donor.index, recipient.index] <- V.modified[recipient.index, donor.index]
		#covariance managed, now to manage the variance
		V.modified[recipient.index, recipient.index] <- (V.original[recipient.index, recipient.index] -  sigma.sq*flow$time.from.root.recipient[flow.index])+ (gamma^2 + (1- gamma)^2) * (flow$time.from.root.recipient[flow.index])*sigma.sq +2*gamma*(1-gamma)*V.original[recipient.index, donor.index]  + vh
        #this is variance for the hybrid. See math derivation at https://github.com/bomeara/bmhyb/issues/1
	}
  if(is.null(measurement.error)) {
	   diag(V.modified) <- diag(V.modified)+SE^2
  } else {
    diag(V.modified) <- diag(V.modified)+measurement.error^2
  }
	return(V.modified)
}

GetMeansModified <- function(x, phy, flow, actual.params) {
	badval<-(0.5)*.Machine$double.xmax
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
	#SE <- x[length(x)]
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]
	}
	times.original <-vcv.phylo(phy, model="Brownian") #is the initial one based on the tree alone, so just time
	V.original <- sigma.sq * times.original

	means.original <- rep(mu, Ntip(phy))
	names(means.original) <- rownames(V.original)
	means.modified <- means.original

	means.original <- rep(mu, Ntip(phy))
	names(means.original) <- rownames(V.original)
	means.modified <- means.original
	for (flow.index in sequence(dim(flow)[1])) {
		recipient.index <- which(rownames(V.original)==flow$recipient[flow.index])
		if(length(recipient.index)!=1) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but instead found ", paste(rownames(V.original)[recipient.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		donor.index <- which(rownames(V.original)==flow$donor[flow.index])
		if(length(donor.index)!=1) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but instead found ", paste(rownames(V.original)[donor.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", sep=""))
		}
		means.modified[recipient.index] <- means.original[recipient.index] + log(bt)
	}
	return(means.modified)
}

#precision is the cutoff at which we think the estimates become unreliable due to ill conditioned matrix
CalculateLikelihood <- function(x, data, phy, flow, actual.params, precision=2, proportion.mix.with.diag=0, allow.extrapolation=FALSE, measurement.error=NULL, do.kappa.check=FALSE, number.of.proportions=101, lower.b=c(0, -Inf, 0.000001, 0, 0), upper.b=c(10,Inf,100,100,100), badval.if.not.positive.definite=TRUE, do.Brissette.correction=FALSE, do.Higham.correction=TRUE, do.DE.correction=FALSE, return.penalty=FALSE, ...) {
	badval<-(0.5)*.Machine$double.xmax
#  x <- ConvertExpm1(x)
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
  SE <- 0
  likelihood.penalty <- 0
  if(is.null(measurement.error)) {
	   SE <- x[length(x)]
   }
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]
	}
  if(any(x<lower.b) | any(x>upper.b)) {
    return(badval)
  }
	V.modified <- GetVModified(x, phy, flow, actual.params, measurement.error=measurement.error)
  if(do.Brissette.correction) {
    V.modified <- BrissetteEtAlCorrection(V.modified)
    if(is.null(V.modified)) {
      return(badval)
    }
  }
  if(do.Higham.correction & !IsPositiveDefinite(V.modified)) {
    new.mat <- as.matrix(Matrix::nearPD(V.modified, corr=FALSE, posd.tol = 1e-16, eig.tol = 1e-16, conv.tol = 1e-16)$mat)
    if(any(new.mat!=V.modified)) {
      warning("Had to do Higham (2002) correction for not positive definite matrix")
      likelihood.penalty <- 10+dist(rbind(c(new.mat), c(V.modified)))
    }
    V.modified <- new.mat
  #  if(min(V.modified)<0) {
  #    warning("Had to remove negative values in VCV after Higham (2002) correction")
  #    V.modified[(V.modified<0)] <- 0
  #  }
  }
  if(return.penalty) {
    return(likelihood.penalty)
  }
  if(do.DE.correction & !IsPositiveDefinite(V.modified)) {
    warning("Have to modify variance covariance matrix to make it positive definite, so results are approximate and the analysis will be slow.")
    V.modified <-  AlterMatrixUsingDE(V.modified)
  }
  if(badval.if.not.positive.definite) {
    if(!IsPositiveDefinite(V.modified)) {
      return(badval)
    }
  }
	means.modified <- GetMeansModified(x, phy, flow, actual.params)
	if(sigma.sq <0 || vh<0 || bt <= 0.0000001 || SE < 0) {
    	return(badval)
	}
  #if(min(eigen(V.modified)$values) < 0.00000001) {
  #  return(0.5*badval) #b/c it might not be very bad, just hard to calculate
  #}
	data <- data[match(names(means.modified), names(data))]
	if(length(data)!=length(means.modified)) {
		stop("Mismatch between names of taxa in data vector and on phy")
	}
	#if(length(which(eigen(V.modified)$values<0))>0) {
	#	last.bad <- V.modified
	#	return(badval)
	#}
	#NegLogML <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified)))
  NegLogML <- NULL
  try(NegLogML <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified) + (1/2)*determinant(V.modified, logarithm=TRUE)$modulus + likelihood.penalty, silent=TRUE)
  if(is.null(NegLogML)) {
    NegLogML <- badval
  }
  #print(paste0("NegLogML = ", NegLogML-likelihood.penalty, ", penalty=", likelihood.penalty))
  # NegLogML.dmvnorm <- -dmvnorm(x=data, mean=means.modified, sigma=V.modified, log=TRUE)
  # if(!is.finite(NegLogML.dmvnorm)) {
  #
  #   print("sigma.sq, vh, bt, SE, mu")
  #   print(c(sigma.sq, vh, bt, SE, mu))
  #   print("us vs dmvnorm")
  #   print(c(NegLogML, NegLogML.dmvnorm))
  # }
	#if(min(V.modified)<0 || sigma.sq <0 || vh<0 || bt <= 0.0000001 || !is.finite(NegLogML) || SE<0) {
  if( sigma.sq <0 || vh<0 || bt <= 0.0000001 || !is.finite(NegLogML) || SE<0) {
  	NegLogML<-badval
      print("badval")
      print(paste0("min(V.modified) ",min(V.modified)))
      print(paste0("sigma.sq ", sigma.sq))
      print(paste0("vh ", vh))
      print(paste0("bt ", bt))
	}
	matrix.condition <- kappa(V.modified, exact=TRUE)
	#print("condition")
	#print(kappa(V.modified, exact=TRUE))
	#print("log(condition)")
	#print(log(kappa(V.modified, exact=TRUE)))

	#pretty<-c(NegLogML, log(matrix.condition))
	#names(pretty) <- c("NegLogL", "log(matrix.condition")
	#print(pretty)
	#The ratio  of the largest to smallest singular value in the singular value decomposition of a matrix. The base- logarithm of  is an estimate of how many base- digits are lost in solving a linear system with that matrix. In other words, it
	#estimates worst-case loss of precision. A system is said to be singular if the condition number is infinite, and ill-conditioned if it is too large, where "too large" means roughly  the precision of matrix entries.
	#if(rcond(V.modified) < .Machine$double.eps^2){
	if(log(matrix.condition) > precision & do.kappa.check) {
    proportions <- seq(from=1, to=0, length.out=number.of.proportions)
		lnl.vector <- rep(NA, length(proportions))
		max.diff <- 0
    kappa.vector <- rep(NA, length(proportions))
		for(i in sequence(length(proportions))) {
			V.modified.by.proportions<-(1-proportions[i]) * V.modified + proportions[i] * diag(dim(V.modified)[1]) * diag(V.modified)
      #local.lnl <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified.by.proportions)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified.by.proportions)))
      local.lnl <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified.by.proportions)%*%(data-means.modified) + (1/2)*determinant(V.modified.by.proportions, logarithm=TRUE)$modulus
      lnl.vector[i] <- local.lnl
      kappa.vector[i] <- kappa(V.modified.by.proportions, exact=TRUE)
			# if(i>6) {
			# 	very.local.lnl <- lnl.vector[(i-6):(i-1)]
			# 	max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)])) #looking locally for jumps in the likelihood
			# 	current.diff <- abs(local.lnl - lnl.vector[i-1])
			# 	if(current.diff > 2 * max.diff) {
			# 		#print(paste("breaking after ", i))
			# 		#break() #the modified matrix is still poorly conditioned, so stop here
      #     lnl.vector[i] <- NA
      #     break()
			# 	}
			# }
		}
    # small.kappas <- which(kappa.vector < 2)
    # if(length(small.kappas)<10) {
    #   small.kappas<-order(kappa.vector, decreasing=FALSE)[1:10]
    # }
    #smoothing.estimate <- smooth.spline(x=proportions[small.kappas], y=lnl.vector[small.kappas], w=1/kappa.vector[small.kappas], df=4)$spar
		proportions<-proportions[which(!is.na(lnl.vector))]
		lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
		NegLogML.4 <- predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=4), data.frame(proportions =0.000))$y #df of 4: don't want overfitting
    NegLogML.6 <- predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=6), data.frame(proportions =0.000))$y #df of 4: don't want overfitting
    if(abs(NegLogML.4 - NegLogML.6)>0.1) { #the result is very sensitive to the degrees of freedom.
      NegLogML <- badval
    } else {
      NegLogML <- NegLogML.4
    }

    #plot(proportions, lnl.vector)
    #lines(proportions, predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=4), data.frame(proportions))$y[,1])
		#plot(c(0, proportions), c(NegLogML, lnl.vector), type="n")
		#points(proportions, lnl.vector, pch=20)
		#points(0, NegLogML, col="red")
		if(abs(NegLogML - lnl.vector[length(lnl.vector)]) > 0.001) {	#means this point was extrapolated b/c the likelihood surface got strange
			if(allow.extrapolation) {
				warning("VCV matrix was ill-conditioned, so used splines to estimate its likelihood (allow.extrapolation=TRUE). This could lead to very bad estimates of the likelihood")
			} else {
				NegLogML <- badval
				warning("VCV matrix was ill-conditioned at some points in parameter space; treated the likelihood at these points as a bad value rather than estimating it. If you think the MLE is in this region, you could try allow.extrapolation=TRUE")
			}
		}

		#print(paste("Did interpolation, got ", NegLogML))
	}


	#print("datadiff")
	#print(quantile(data-means.modified))
	#print("middle")
	#print((1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified))
	#print("end")
	#print((1/2)*log(abs(det(V.modified))) )
	#print(x)
	#print(V.modified[1:10,1:10])
	#print(means.modified)
  # if(NegLogML < -1000 | NegLogML > 100000) {
  #   stop("bad")
  # } else {
  #   stop("good")
  # }
	#if(NegLogML< (0)) {
	#	NegLogML <- badval #since something seems off.
	#}

	return(NegLogML[1])

}

ConvertLog1P <- function(x) {
  x[which(names(x)=="mu")] <- log1p(x[which(names(x)=="mu")])
}

ConvertExpm1 <- function(x) {
  x[which(names(x)=="mu")] <- expm1(x[which(names(x)=="mu")])
}

AdaptiveConfidenceIntervalSampling <- function(par, fn, lower=-Inf, upper=Inf, desired.delta = 2, n.points=5000, verbose=TRUE, measurement.error=NULL, do.kappa.check=FALSE, allow.restart=TRUE,  best.lnl = -Inf, likelihood.precision=0.001, restart.mode=FALSE, ...) {
	starting<-fn(par, measurement.error=measurement.error,  ...)
	if(length(lower) < length(par)) {
		lower<-rep(lower, length(par))
	}
	if(length(upper) < length(par)) {
		upper<-rep(upper, length(par))
	}
	min.multipliers <- rep(1, length(par))
	max.multipliers <- rep(1, length(par))
	results<-data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
	results[1,]<-unname(c(starting, par))
	for (i in sequence(n.points)) {
		sim.points<-NA
		while(is.na(sim.points[1]) | !is.numeric(sim.points[1])) {
			sim.points<-GenerateValues(par, lower, upper, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
		}
		results[i+1,] <- c(fn(sim.points, measurement.error=measurement.error, do.kappa.check=do.kappa.check, ...), sim.points)
    if(i>5 & restart.mode) {
      if((best.lnl - min(results[,1], na.rm=TRUE) > likelihood.precision ) & allow.restart) {
        results <- results[sequence(i+1),] #stop here and restart
        return(results)
      }
    }
		if (i%%20==0) {
			for (j in sequence(length(par))) {
				returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
				total.range <- range(results[,j+1], na.rm=TRUE)
				width.ratio <- diff(returned.range)/diff(total.range)
				if(is.na(width.ratio)) {
					width.ratio=1
				}
				if(width.ratio > 0.5) { #we are not sampling widely enough
					min.multipliers[j] <- min.multipliers[j] * 0.9
					max.multipliers[j] <- max.multipliers[j] * 1.1 #expand the range
				} else {
					min.multipliers[j] <- 1
					max.multipliers[j] <- 1
				}
			}
		}
		if (verbose && i%%100==0) {
			print(paste(i, "of", n.points, "done"))
		}
	}
	return(results)
}

GenerateRandomValues <- function(data, free.parameters, lower, upper) {
  new.vals <- c(rexp(1,10), runif(1, min=min(data), max=max(data)), rexp(1,1), rexp(1,10), rexp(1,10))[which(free.parameters)]
  while(any(new.vals < lower) | any(new.vals>upper)) {
    new.vals <- c(rexp(1,10), runif(1, min=min(data), max=max(data)), rexp(1,1), rexp(1,10), rexp(1,10))[which(free.parameters)]
  }
  return(new.vals)
}

GenerateValues <- function(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) {
  if(is.null(lower)) {
    lower <- 0.1*par
  }
  if(is.null(upper)) {
    upper <- 10*par
  }
	pass=FALSE
	tries=0
	while(!pass && tries<=max.tries) {
		tries <- tries+1
		pass=TRUE
		new.vals <- rep(NA, length(par))
		for(i in sequence(length(par))) {
			examined.max[i]<-max(0.001, examined.max[i])
      new.vals.bounds <- sort(c(max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i])), decreasing=FALSE)
			new.vals[i]<-runif(1, min=ifelse(is.finite(new.vals.bounds[1]),new.vals.bounds[1], 0.000001) , max=ifelse(is.finite(new.vals.bounds[2]), new.vals.bounds[2], 10000))

			if(new.vals[i]<lower[i]) {
				pass=FALSE
			}
			if(new.vals[i]>upper[i]) {
				pass=FALSE
			}
		}
	}
	if(tries>max.tries) {
		return(NA)
	}
	return(new.vals)
}

GetClade <- function(phy, clade.size) {
	nodes <- phy$edge[,1]
	subtrees <- lapply(nodes, extract.clade, phy=phy)
	counts <- sapply(subtrees, Ntip)
	matches<-subtrees[which(counts==clade.size)]
	if(length(matches)==0) {
		return(NA)
	}
	lucky <- matches[sample.int(length(matches),1)][[1]]
	return(findMRCA(phy, tips=lucky$tip.label, type="node"))
}

GetAncestor <- function(phy, node) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}


AttemptDeletionFix <- function(phy, flow, params=c(1,0,0.1, 0, 0), m.vector = c(1,2)) {
  if(is.null(names(params))) {
    names(params) <- c("bt", "vh", "sigma.sq", "mu", "SE")
  }
  taxa.to.try.deleting <- phy$tip.label
  taxa.to.try.deleting <- taxa.to.try.deleting[!(taxa.to.try.deleting %in% flow$recipient)]
  taxa.to.try.deleting <- taxa.to.try.deleting[!(taxa.to.try.deleting %in% flow$donor)]
  taxa.to.try.deleting <- taxa.to.try.deleting[sample.int(length(taxa.to.try.deleting), length(taxa.to.try.deleting), replace=FALSE)]
  if(length(taxa.to.try.deleting)==0) {
    stop("There are no taxa to delete that aren't involved in hybridization.")
  }
  phy.pruned <- phy
  current.m.index <- 1
  current.index <- 1
  combos.to.delete <- utils::combn(taxa.to.try.deleting,m.vector[current.m.index])
  while(!IsPositiveDefinite(GetVModified(params, phy.pruned, flow, actual.params= rep(TRUE,length(params))))) {
    #print(current.index)
    #print(paste0("Dropping ", paste(combos.to.delete[,current.index], collapse=" ")))
    phy.pruned <- ape::drop.tip(phy, combos.to.delete[,current.index])
    current.index <- current.index + 1
    if(current.index > ncol(combos.to.delete)) {
      if(current.m.index < length(m.vector)) {
        current.m.index <- current.m.index + 1
        combos.to.delete <- utils::combn(taxa.to.try.deleting,m.vector[current.m.index])
        current.index <- 1
      } else {
        stop(paste0("Correction by removing random taxa up to ", max(m.vector), " taxa at a time failed."))
      }
    }
  }
  return(phy.pruned)
}

#allow.ghost allows ghost lineage: something that persists for awhile, hybridizes, goes extinct. Otherwise, hybridization events must between coeval edges with extant descendants
SimulateNetwork <- function(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 1, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE) {
	done = FALSE
	used.recipients <- c()
	available.recipient.ids <- sequence(ntax.nonhybrid + ntax.hybrid)
	flow <- data.frame()
	phy<-NA
	phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=death, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
	if(origin.type=="clade" && ntax.hybrid==1) {
		warning("For ntax.hybrid = 1 and clade sampling, this will do individual sampling instead (which is equivalent in this case)")
		origin.type<-"individual"
	}
	if(origin.type=="clade") {
		while(is.na(GetClade(phy, ntax.hybrid))) { #not all trees of a given size have a clade of a given size, so may need to resimulate it
			phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=0.5, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
		}
	}
	while(!done) {
		donors <- c()
		recipients <- c()
		recipient.ids <- c()
		focal.node <- c()
		if (origin.type=="clade") {
			focal.node <- GetClade(phy, ntax.hybrid)
			if(is.na(focal.node)) {
				done=FALSE
				break()
			}
			recipients <- phy$tip.label[getDescendants(phy, node=focal.node)]
			recipients <- recipients[!is.na(recipients)] #since we want just the tips
			recipient.ids <- which(phy$tip.label %in% recipients)
			used.recipients <- append(used.recipients, recipients)
		} else {
			focal.node<-sample(available.recipient.ids, 1, replace=FALSE)
			recipient.ids <- focal.node
			recipients <- phy$tip.label[focal.node]
			used.recipients <- append(used.recipients, recipients)
		}
		available.recipient.ids <- available.recipient.ids[!available.recipient.ids %in% recipient.ids]
		longest.from.root <- nodeheight(phy, node=focal.node)
		shortest.from.root <- nodeheight(phy, node=GetAncestor(phy, focal.node))
		all.heights <- nodeHeights(phy)
		#idea here: take a recipient clade. The gene flow must happen on its stem edge, which starts at shortest.from.root and goes up to longest.from.root. Gene flow can't go back in time
		qualifying.lower <- which(all.heights[,1]<longest.from.root) #if this is false, gene flow goes back in time
		qualifying.upper <- sequence(dim(all.heights)[1]) #in general, gene flow can go forward in time via ghost lineages
		if(!allow.ghost) {
			qualifying.upper <- which(all.heights[,2]>shortest.from.root) #if no ghost lineages, then there must be temporal overlap between the donor and recipient lineages. So the tipward end of the donor edge must be later than the rootward end of the recipient edge
		}
		qualifying.upper <- qualifying.upper[which(phy$edge[qualifying.upper,2]!=focal.node)] #let's not hybridize with ourselves
		qualifying.all <- qualifying.upper[qualifying.upper %in% qualifying.lower]
		if(length(qualifying.all)==0) {
			break()
		}
		donor.edge <- sample(qualifying.all, 1)
		donors <- phy$tip.label[getDescendants(phy, phy$edge[donor.edge,2])]
		donors <- donors[!is.na(donors)] #getDescendants includes all descendant nodes, including internal ones. We just want the terminal taxa
		time.in <- runif(1, min=max(all.heights[donor.edge,1],shortest.from.root), max=longest.from.root)
		time.out <- runif(1, min=all.heights[donor.edge,1], max=min(time.in, all.heights[donor.edge,2]))
		if (!allow.ghost) {
			time.in <- runif(1, min=max(shortest.from.root, all.heights[donor.edge,1]), max=min(longest.from.root, all.heights[donor.edge,2])) #if no ghost lineages, must move from the overlapping interval
			time.out <- time.in
		}
		pairs <- expand.grid(donors, recipients)
		for (pairs.index in sequence(dim(pairs)[1])) {
			flow <- rbind(flow, data.frame(donor=pairs[pairs.index,1], recipient=pairs[pairs.index,2], gamma=flow.proportion, time.from.root.donor=time.out, time.from.root.recipient=time.in, stringsAsFactors=FALSE))
		}
		if(length(used.recipients)==ntax.hybrid) {
			done=TRUE
		}
	}
	flow$donor <- as.character(flow$donor)
	flow$recipient <- as.character(flow$recipient)
	flow$gamma <- as.numeric(as.character(flow$gamma))
	flow$time.from.root.donor <-as.numeric(as.character(flow$time.from.root.donor))
	flow$time.from.root.recipient <-as.numeric(as.character(flow$time.from.root.recipient))
	return(list(phy=phy, flow=flow))
}

PlotNetwork <- function(phy, flow, col.non="black", col.hybrid="red", col.donor="blue", name.padding=1.5, cex=1, xlab="", bty="n", head.length=0.2, edge.width=2, col.tree="darkgray", col.arrow="red", arrow.width=1, try.rotations=FALSE, ...) {
	phy<-reorder(phy, "pruningwise")
  if(try.rotations) {
    phy <- phytools::rotateNodes(phy, nodes="all")
  }
	phy4 <- as(phy, "phylo4")
	xxyy <- phyloXXYY(phy4)
	#plot(phy4)
	plot(x=c(min(xxyy$xx), name.padding*max(xxyy$xx)), y=range(xxyy$yy), type="n", xaxt="n", xlab=xlab, yaxt="n", ylab="", bty=bty, ...)
	arrows(x0=xxyy$segs$v0x, x1=xxyy$segs$v1x, y0=xxyy$segs$v0y, y1=xxyy$segs$v1y, length=0, lwd=edge.width, col=col.tree)
	arrows(x0=xxyy$segs$h0x, x1=xxyy$segs$h1x, y0=xxyy$segs$h0y, y1=xxyy$segs$h1y, length=0, lwd=edge.width, col=col.tree)
	label.colors <- rep(col.non, Ntip(phy))
	for (i in sequence(Ntip(phy))) {
		if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$donor) {
			label.colors[i]<-col.donor
		}
		if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$recipient) {
			label.colors[i]<-col.hybrid
		}

	}
	text(x=rep(max(xxyy$xx), Ntip(phy)), y=xxyy$yy[which(edges(phy4)[xxyy$eorder,2] %in% sequence(Ntip(phy)))], names(getNode(phy4, xxyy$torder)), col=label.colors, pos=4, cex=cex)
	for (i in sequence(dim(flow)[1])) {
		recipient.node <- getNode(phy4, flow$recipient[i])
		recipient.path <- c(recipient.node, ancestors(phy4, recipient.node))
		recipient.path.heights <- nodeHeight(phy4, recipient.path, from="root")
		valid.recipients <- recipient.path[which(recipient.path.heights > flow$time.from.root.recipient[i])]
		recipient <- valid.recipients[length(valid.recipients)] #do it from the earliest qualifying tipward node
		if(length(valid.recipients)>1 && length(which(recipient.path.heights==flow$time.from.root.recipient[i]))>0) { #the latter condition means we're moving to an existing node
				recipient <- valid.recipients[length(valid.recipients)-1]
		}
		y1 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == recipient)]
		donor.node <- getNode(phy4, flow$donor[i])
		donor.path <- c(donor.node, ancestors(phy4, donor.node))
		donor.path.heights <- nodeHeight(phy4, donor.path, from="root")
		valid.donors <- donor.path[which(donor.path.heights > flow$time.from.root.donor[i])]
		donor <- valid.donors[length(valid.donors)] #do it from the earliest qualifying tipward node
		y0 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == donor)]
		arrows(x0=flow$time.from.root.donor[i]/max(branching.times(phy)), x1=flow$time.from.root.recipient[i]/max(branching.times(phy)), y1=y1, y0=y0, col=col.arrow, lwd=arrow.width, length=head.length) #rescale since it goes from zero to 1 in height
		#grid.arrows(x=c(flow$time.from.root[i],flow$time.from.root[i]), y=c(y0, y1))
	}
}

LumpIntoClades <- function(phy, flow) {
	flow.string <- paste(flow$time.from.root.donor, flow$time.from.root.recipient, flow$gamma)
	flow.together <- cbind(flow.string=flow.string, flow)
	unique.flow.strings <- unique(flow.string)
	donor.clades <- c()
	recipient.clades <- c()
	gamma.clades <-c()
	time.from.root.donor.clades <- c()
	time.from.root.recipient.clades <- c()
	for (i in sequence(length(unique.flow.strings))) {
		flow.local <- flow.together[which(flow.together$flow.string == unique.flow.strings[i]), ]
		donor.clades <- append(donor.clades, paste(unique(flow.local$donor), collapse=","))
		recipient.clades <- append(recipient.clades, paste(unique(flow.local$recipient), collapse=","))
		gamma.clades <-append(gamma.clades, flow.local$gamma[1])
		time.from.root.donor.clades <- append(time.from.root.donor.clades, flow.local$time.from.root.donor[1])
		time.from.root.recipient.clades <- append(time.from.root.recipient.clades, flow.local$time.from.root.recipient[1])

	}
	return(data.frame(donor.clades = donor.clades, recipient.clades=recipient.clades, gamma=gamma.clades, time.from.root.donor=time.from.root.donor.clades, time.from.root.recipient = time.from.root.recipient.clades, stringsAsFactors=FALSE))
}

AttachHybridsToDonor <- function(phy, flow, suffix="_DUPLICATE") {
	flow.clades <- LumpIntoClades(phy, flow)
	phy.merged <- phy
	for (i in sequence(dim(flow.clades)[1])) {
		taxa.to.retain <- strsplit(flow.clades$recipient.clades[i], ",")[[1]]
		donor.taxa <- strsplit(flow.clades$donor.clades[i], ",")[[1]]
		pulled.clade <- c()
		if(length(taxa.to.retain)>1) {
			pulled.clade <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% taxa.to.retain])
		} else {
			pulled.clade <- structure(list(edge = structure(c(3L, 3L, 1L, 2L), .Dim = c(2L,
2L)), edge.length = rep(flow.clades$time.from.root.recipient[i],2),
    tip.label = c(taxa.to.retain, "DUMMY"), Nnode = 1L), .Names = c("edge",
"edge.length", "tip.label", "Nnode"), class = "phylo")
		}
		pulled.clade$tip.label <- paste(pulled.clade$tip.label, suffix, sep="")
		attachment.crown.node <- which(phy.merged$tip.label==donor.taxa[1])
		if(length(donor.taxa)>1) {
			attachment.crown.node <- findMRCA(phy.merged, tips=donor.taxa, type=c("node"))
		}
		attachment.stem.node <- GetAncestor(phy.merged, attachment.crown.node)

		pulled.clade$root.edge<-max(branching.times(phy)) - max(branching.times(pulled.clade)) - flow.clades$time.from.root.donor[i]

		phy.merged <- bind.tree(phy.merged, pulled.clade, attachment.crown.node, position=phy.merged$edge.length[which(phy.merged$edge[,2]==attachment.crown.node)] - (flow.clades$time.from.root.donor[i]-nodeheight(phy.merged, attachment.stem.node)))
		if(length(taxa.to.retain)==1) {
			phy.merged <- drop.tip(phy.merged, paste("DUMMY", suffix, sep=""))
		}
	}
	return(phy.merged)
}



#		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#GetMeansModified <- function(x, phy, flow, actual.params) {
#params must be named vector
SimulateTipData <- function(phy, flow, params, suffix="_DUPLICATE") {
#	flow.clades <- LumpIntoClades(phy, flow)
#	recipients <- c()
#	for (i in sequence(dim(flow.clades)[1])) {
#		recipients<-c(recipients, strsplit(flow.clades$recipient.clades, ",")[[1]])
#	}
#	if(length(unique(recipients)) != length(recipients)) {
#		print(flow)
#		print(flow.clades)
#		stop("This function only works if each taxon of hybrid origin only appears once as a recipient")
#	}
	phy.merged <- AttachHybridsToDonor(phy, flow, suffix=suffix)
	phy.merged$edge.length <- phy.merged$edge.length*params["sigma.sq"]
	phy.merged$edge.length[which(phy.merged$edge[,2] <= Ntip(phy.merged))] <- phy.merged$edge.length[which(phy.merged$edge[,2] <= Ntip(phy.merged))]+params["SE"]
	tips <- sim.char(phy.merged, par=matrix(1, nrow=1, ncol=1), nsim=1, model="BM", root=params["mu"])[,,1]	#we have already done the scaling and SE
	hybrid.name.root <- gsub(suffix, "", phy.merged$tip.label[grepl(suffix,  phy.merged$tip.label)])
	for (i in sequence(length(hybrid.name.root))) {
		focal.tips <- tips[c(hybrid.name.root[i], paste(hybrid.name.root[i], suffix, sep=""))]
		focal.gamma <- flow$gamma[which(flow$recipient==hybrid.name.root[i])]
		focal.tips.bt <- focal.tips + log(params['bt'])
		tip.mean <- focal.gamma*focal.tips.bt[2] + (1-focal.gamma) * focal.tips.bt[1]
		tip.final <- rnorm(1, mean=tip.mean, sd=sqrt(params['vh']))
		tips[hybrid.name.root]<-tip.final
	}
	tips<-tips[!grepl(suffix, names(tips))]
	return(tips)
}

#The following short function comes from Ken Takagi at  https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
PlotConvexHull<-function(xcoord, ycoord, lcolor){
	hpts <- chull(x = xcoord, y = ycoord)
	hpts <- c(hpts, hpts[1])
	lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}

ContourFromAdaptiveSampling<-function(sims, params.of.interest=NULL) {
	#sims is a data.frame with the sim results
	if(is.null(params.of.interest)) {
		params.of.interest <- colnames(sims)[-1]
	}
	sims$neglnL <- sims$neglnL - min(sims$neglnL)
	for (param.1 in sequence(length(params.of.interest)-1)) {
		for	(param.2 in c((param.1+1) : length(params.of.interest))) {
			sims.sub <- sims[,c(params.of.interest[param.1], params.of.interest[param.2], 'neglnL')]
			# points.to.fit <- data.frame()
			x.range <- range(sims.sub[which(sims.sub[,3]<10),1])
			y.range <- range(sims.sub[which(sims.sub[,3]<10),2])
			x.grid <- quantile(sims.sub[,1], seq(from=0, to=1, length.out = floor(length(sims.sub[,1])/50)))
			y.grid <- quantile(sims.sub[,2], seq(from=0, to=1, length.out = floor(length(sims.sub[,2])/50)))
			colnames(sims.sub)[1:2] <- c("x", "y")
			# for (x.index in sequence(length(x.grid)-1)) {
				# for (y.index in sequence(length(y.grid)-1)) {
					# relevant.points <- subset(sims.sub, x>=x.grid[x.index] & x<x.grid[x.index+1] & y>=y.grid[y.index] & y<y.grid[y.index+1])
					# points.to.fit <- rbind(points.to.fit, relevant.points[which.max(relevant.points[,3]),])
				# }
			# }
			plot(x=x.range, y=y.range, xlab=params.of.interest[param.1], ylab=params.of.interest[param.2], type="n", bty="n")

			relevant.point.id <- chull(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2])
			#polygon(sims.sub[relevant.point.id,1], sims.sub[relevant.point.id, 2], col="gray", border=NA, fillOddEven=TRUE)
			#polygon(sims.sub[relevant.point.id,1], sims.sub[relevant.point.id, 2], col="gray", border=NA, fillOddEven=FALSE)
			PlotConvexHull(sims.sub[which(sims.sub[,3]<5),1], sims.sub[which(sims.sub[,3]<5),2], "darkgray")
      PlotConvexHull(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2], "black")

			#points(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2], col="green", pch="X")
						#contour(interp(points.to.fit[,1], points.to.fit[,2], points.to.fit[,3]), xlab=params.of.interest[param.1], ylab=params.of.interest[param.2], levels=c(1, 2, 5, 10))
			points(sims.sub[which.min(sims.sub$neglnL),1], sims.sub[which.min(sims.sub$neglnL),2], col="red", pch=20, cex=2)
		}
	}
}

ConvertVectorToMatrix <- function(x) {
  matrix.size <- (-1 + sqrt(1+8*length(x)))/2
  new.mat <- matrix(NA, ncol=matrix.size, nrow=matrix.size)
  new.mat[upper.tri(new.mat, diag=TRUE)] <- x
  new.mat[lower.tri(new.mat, diag=FALSE)] <- new.mat[upper.tri(new.mat, diag=FALSE)]
  return(new.mat)
}

PositiveDefiniteOptimizationFn <- function(x, original) {
  new.mat <- ConvertVectorToMatrix(x)
  distance <- as.numeric(dist(rbind(as.vector(new.mat), as.vector(original))))
  # new.mat.no.diag <- new.mat
  # diag(new.mat.no.diag) <- 0
  # original.no.diag <- original
  # diag(original.no.diag) <- 0
  # distance <- distance + 10 * as.numeric(dist(rbind(as.vector(new.mat.no.diag), as.vector(original.no.diag ))))
  if(min(new.mat)<0) {
    neg.values <- new.mat[which(new.mat<0)]
    distance <- distance * (1+abs(sum(neg.values)))
  }
  if(!IsPositiveDefinite(new.mat)) {
    distance <- distance * max(c(sum(new.mat), 10))
  }
  return(distance)
}

#From Ravi Varadhan, http://r.789695.n4.nabble.com/how-to-randomly-generate-a-n-by-n-positive-definite-matrix-in-R-td846858.html
GenerateRandomPositiveDefiniteMatrix <- function(n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}


# This is inspired by the work of
#Mishra, Sudhanshu K. "The nearest correlation matrix problem: Solution by differential evolution method of global optimization." (2007)
# https://mpra.ub.uni-muenchen.de/44809/9/MPRA_paper_44809.pdf
AlterMatrixUsingDE <- function(V.modified) {
  starting.val.center <- V.modified[upper.tri(V.modified, diag=TRUE)]
  #starting.val.matrix <- matrix(NA, nrow=max(40,round(5*sqrt(length(starting.val.center)))), ncol=length(starting.val.center))
  starting.val.matrix <- matrix(NA, nrow=40, ncol=length(starting.val.center))
  starting.val.matrix[1,] <- starting.val.center
  starting.means <- log(starting.val.center)
  starting.means[!is.finite(starting.means)] <- min(starting.means[is.finite(starting.means)])
  sd.vector <- seq(from=0.001, to=1, length.out=nrow(starting.val.matrix)) #to give some points close to original, some further away
  diag.matrix <- matrix(0, ncol(V.modified), nrow(V.modified))
  diag(diag.matrix) <- max(V.modified)
  starting.val.matrix[2,] <- diag.matrix[upper.tri(diag.matrix, diag=TRUE)] #something that will start positive definite
  diag.matrix.2 <- matrix(0, nrow=nrow(V.modified), ncol=ncol(V.modified))
  diag(diag.matrix.2) <- diag(V.modified)
  starting.val.matrix[3,] <- diag.matrix.2[upper.tri(diag.matrix.2, diag=TRUE)] #something that will start positive definite
  diag.matrix.3 <- V.modified
  diag(diag.matrix.3) <- mean(diag(V.modified))
  starting.val.matrix[4,] <- diag.matrix.3[upper.tri(diag.matrix.3, diag=TRUE)] #something that will start positive definite
  pos.def.candidate <- as.matrix(Matrix::nearPD(V.modified, corr=FALSE)$mat)
  pos.def.eigen <- eigen(pos.def.candidate)$values
  has.neg <- function(x) {
    return(any(x<0))
  }
  negative.taxa <- unique(c(which(apply(pos.def.candidate, 1,has.neg)), which(apply(pos.def.candidate, 2, has.neg))))
  pos.def.values <- pos.def.candidate[upper.tri(pos.def.candidate, diag=TRUE)]
  pos.def.values.abs <- abs(pos.def.values)
  pos.def.values.zeroed <- pos.def.values
  pos.def.values.zeroed <- pos.def.values.zeroed[which(pos.def.values<0)] <- 0
  starting.val.matrix[5,] <- pos.def.values #start with a positive definite matrix (but might not meet the nonnegative constraint)
  starting.val.matrix[6,] <- pos.def.values.abs #start with a potentially positive definite matrix (but might not be, since we've converted neg to positive values)
  starting.val.matrix[7,] <- pos.def.values.zeroed #start with a potentially positive definite matrix (but might not be, since we've converted neg to zero)
  number.filled <- sum(!is.na(starting.val.matrix[,1]))
  for (i in (number.filled+1):(number.filled+10)) {
    local.mat <- GenerateRandomPositiveDefiniteMatrix(ncol(V.modified), ev=pos.def.eigen)
    starting.val.matrix[i,] <- local.mat[upper.tri(local.mat, diag=TRUE)]
  }
  number.filled <- sum(!is.na(starting.val.matrix[,1]))
  proportions <- seq(from=1, to=0, length.out=10)
  for(i in sequence(length(proportions))) {
    V.modified.by.proportions<-(1-proportions[i]) * V.modified + proportions[i] * diag(dim(V.modified)[1]) * diag(V.modified)
    starting.val.matrix[i,] <- V.modified.by.proportions[upper.tri(V.modified.by.proportions, diag=TRUE)]
  }
  number.filled <- sum(!is.na(starting.val.matrix[,1]))
  for (i in (number.filled+1):nrow(starting.val.matrix)) {
    #starting.val.matrix[i,] <- rexp(length(starting.val.center), rate=starting.rates)
    starting.val.matrix[i,] <- rlnorm(length(starting.val.center), meanlog=starting.means, sdlog=sd.vector[i])
  }
  result <- DEoptim::DEoptim(PositiveDefiniteOptimizationFn, lower=rep(0, sum(upper.tri(V.modified, diag=TRUE))), upper=rep(2*max(V.modified), sum(upper.tri(V.modified, diag=TRUE))), control=list(trace=FALSE, initialpop = starting.val.matrix, c=0.1, itermax=20, reltol=1e-1), original=V.modified)
  #result <- DEoptim::DEoptim(PositiveDefiniteOptimizationFn, lower=rep(0, sum(upper.tri(V.modified, diag=TRUE))), upper=rep(2*max(V.modified), sum(upper.tri(V.modified, diag=TRUE))), control=list(trace=FALSE, initialpop = starting.val.matrix, c=0.1, itermax=20, reltol=1e-1, strategy=6, p=0.3), original=V.modified)

  #result <- optim(par=starting.val.center, fn=PositiveDefiniteOptimizationFn, lower=rep(0, sum(upper.tri(V.modified, diag=TRUE))), upper=rep(2*max(V.modified), sum(upper.tri(V.modified, diag=TRUE))), method="L-BFGS-B", original=V.modified)

  #print(paste0("Bestval ",result$optim$bestval, " number of function evals ", result$optim$nfeval, " number of iterations ", result$optim$iter, " smallest value is ", min(ConvertVectorToMatrix(result$optim$bestmem)), " status of being positive definite is ", IsPositiveDefinite(ConvertVectorToMatrix(result$optim$bestmem))))
  #print(V.modified[1:6,1:6])

  #print(ConvertVectorToMatrix(result$par))
  #print(ConvertVectorToMatrix(result$optim$bestmem)[1:6,1:6])
  final.mat <- ConvertVectorToMatrix(result$optim$bestmem)
  rownames(final.mat) <- rownames(V.modified)
  colnames(final.mat) <- colnames(V.modified)
  return(final.mat)
}
