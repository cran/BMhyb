ComputeAICc<-function(n,k,LogLik){
    if((n-k-1) <=0 ) {
        return(NA)
    }
    return(2*n*k/(n-k-1)+2*LogLik)
}

# AICc<-function(n,k,LogLik){
#   return(2*n*k/(n-k-1)+2*LogLik)
# }


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
# BMhyb <- function(data, phy, flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=TRUE, get.se=TRUE, plot.se=TRUE, store.sims=FALSE, precision=2, auto.adjust=FALSE, likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000, measurement.error=NULL, do.kappa.check=FALSE, number.of.proportions=101, number.of.proportions.adaptive=101, allow.restart=TRUE, lower.bounds = c(sigma.sq = 0, mu = -Inf, bt = 1e-06, vh = 0, SE = 0), upper.bounds=c(sigma.sq = 10, mu = Inf,bt = 100,vh = 100,SE = 100), badval.if.not.positive.definite=TRUE, attempt.deletion.fix=FALSE, starting.values=NULL, n.random.start.points=5000, do.Brissette.correction=FALSE, do.Higham.correction=FALSE, do.DE.correction=FALSE) {
#     preset.starting.parameters = NULL
#     flow.problems <- CheckFlow(phy, flow)$problem.taxa
#     if(length(flow.problems)>0) {
#         stop(paste("Sorry, the algorithm cannot work with overlapping hybridization (where any taxon has a history with more than one hybridization event leading to it). In this case, it is multiple events leading to taxon/taxa", paste(flow.problems, collapse=", "), "that are causing the issue. You can edit your flow data.frame manually; you may also use AdjustFlow to randomly delete hybridization events or taxa of hybrid origin."))
#     }
#     if(n.random.start.points>0 & is.null(starting.values)) {
#         grid.results <- BMhybGrid(data=data, phy=phy, flow=flow, verbose=verbose, precision=precision, n.points=n.random.start.points, attempt.deletion.fix=FALSE, measurement.error=measurement.error, get.se=FALSE, plot.se=FALSE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#         starting.values=grid.results$sims[which.min(grid.results$sims$AICc)[1],1:5]
#         starting.names <- colnames(starting.values)
#         starting.values <- as.numeric(starting.values)
#         names(starting.values) <- starting.names
#     }
#     if(min(flow$gamma)<0) {
#         stop("Min value of flow is too low; should be between zero and one")
#     }
#     if(max(flow$gamma)>1) {
#         stop("Max value of flow is too high; should be between zero and one")
#     }
#     results<-list()
#     #hessians <- list()
#     results.summary <-data.frame()
#     if(auto.adjust) {
#         phy <- AdjustForDet(phy)
#     }
#     all.sims<-list()
#
#     starting.from.geiger<-NA
#     #starting.values <- NA
#     if(is.null(starting.values)) {
#         starting.values <- NA
#         if(verbose) {
#             print("Getting starting values from Geiger")
#         }
#         phy.geiger.friendly <- phy #geiger can't handle branch lengths near zero. Let's lengthen them if needed
#         if(min(phy.geiger.friendly$edge.length)<0.00001) {
#             phy.geiger.friendly$edge.length <- phy.geiger.friendly$edge.length + 0.00001
#         }
#         geiger.SE <- data*NA
#         if(!is.null(measurement.error)) {
#             if(length(measurement.error)==1) {
#                 geiger.SE <- rep(measurement.error, length(geiger.SE))
#                 names(geiger.SE) <- phy$tip.label
#             } else {
#                 geiger.SE <- measurement.error
#                 names(geiger.SE) <- phy$tip.label
#
#             }
#             starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
#             starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  0.01*starting.from.geiger$sigsq*max(vcv(phy)), mean(measurement.error)) #sigma.sq, mu, beta, vh, SE
#
#         } else {
#             starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
#             starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  0.01*starting.from.geiger$sigsq*max(vcv(phy)), starting.from.geiger$SE) #sigma.sq, mu, beta, vh, SE
#         }
#         names(starting.values) <- GenerateParamLabels()
#         if(verbose) {
#             print("Done getting starting values")
#         }
#     }
#     if(badval.if.not.positive.definite) {
#         if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, measurement.error=measurement.error))) {
#             if(attempt.deletion.fix) {
#                 phy <- AttemptDeletionFix(phy, flow, starting.values)
#                 tips <- tips[names(tips) %in% phy$tip.label]
#             }
#             #  stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
#         }
#     }
#     for (model.index in sequence(length(models))) {
#         do.run = TRUE
#         preset.starting.parameters = NULL
#         while(do.run) {
#             do.run = FALSE
#             step.count <- 0
#             if(verbose) {
#                 print(paste("Starting model", models[model.index], "of", length(models), "models"))
#             }
#             model <- models[model.index]
#             if(model==1) {
#                 starting.values <- starting.values[-which(names(starting.values=="bt"))]
#             }
#             if(model==2) {
#                 starting.values <- starting.values[-which(names(starting.values=="vh"))]
#             }
#             if(model==3) {
#                 starting.values <- starting.values[-which(names(starting.values=="bt"))]
#                 starting.values <- starting.values[-which(names(starting.values=="vh"))]
#             }
#             if(!is.null(measurement.error)) {
#                 starting.values <- starting.values[-which(names(starting.values=="SE"))]
#             }
#
#             if(is.null(preset.starting.parameters)) {
#                 preset.starting.parameters <- starting.values
#             }
#             # if(badval.if.not.positive.definite) {
#             #   if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, actual.params= free.parameters))) {
#             #     stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
#             #   }
#             # }
#             best.run <- optim(par=preset.starting.parameters, fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#             best.run$par <- best.run$par
#             attempts <- 1
#             while(best.run$convergence!=0){#want to get a convergence code 0
#                 if(verbose) {
#                     print(paste0("Initial search had a convergence code of ", best.run$convergence, ", indicating it did not converge. See ?optim for what the code may mean. Starting again, likely near that point. Negative log likelihood was ", best.run$value))
#                     print("Parameter estimates were")
#                     current.params <- best.run$par
#                     print(current.params)
#                 }
#                 attempts <- attempts+1
#                 if(attempts%%3!=1) {
#                     best.run<-optim(par=GenerateValues(best.run$par, lower=lower.bounds[names(preset.starting.parameters)], upper=upper.bounds[names(preset.starting.parameters)], examined.max=10*best.run$par, examined.min=0.1*best.run$par), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                 } else {
#                     best.run<-optim(par=GenerateRandomValues(data, preset.starting.parameters, lower.bounds[names(preset.starting.parameters)], upper.bounds[names(preset.starting.parameters)]), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                 }
#                 #best.run$par <- ConvertExpm1(best.run$par)
#
#             }
#
#
#
#             if(verbose) {
#                 results.vector<-c(step.count, best.run$value, best.run$par)
#                 names(results.vector) <- c("step","negloglik", names(best.run$par))
#                 print(results.vector)
#             }
#             #this is to continue optimizing; we find that optim is too lenient about when it accepts convergence
#
#             times.without.improvement <- 0
#             while(times.without.improvement<10) {
#                 times.without.improvement <- times.without.improvement+1
#                 new.run <- optim(par=best.run$par, fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                 #new.run$par <- ConvertExpm1(new.run$par)
#                 attempts <- 1
#                 while(new.run$convergence!=0 & attempts < 20){#want to get a convergence code 0
#                     if(verbose) {
#                         print(paste0("This search had a convergence code of ", new.run$convergence, ", indicating it did not converge. See ?optim for what the code may mean (1 = maximum number of steps were hit; 10 means the search optimizer became degenerate). Starting again, likely near that point. Negative log likelihood was ", new.run$value, ". Attempting start ", attempts+1, " of 20 max before we give up"))
#                         print("Parameter estimates were")
#                         current.params <- new.run$par
#                         print(current.params)
#                     }
#                     attempts <- attempts+1
#                     if(attempts%%3!=1) {
#                         new.run<-optim(par=GenerateValues(best.run$par, lower=lower.bounds[names(preset.starting.parameters)], upper=upper.bounds[names(preset.starting.parameters)], examined.max=10*best.run$par, examined.min=0.1*best.run$par), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                     } else {
#                         new.run<-optim(par=GenerateRandomValues(data, preset.starting.parameters, lower.bounds[names(preset.starting.parameters)], upper.bounds[names(preset.starting.parameters)]),  fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, precision=precision, allow.extrapolation=allow.extrapolation, measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                     }
#                     #  new.run$par <- ConvertExpm1(new.run$par)
#
#                 }
#
#                 #print("new.run best.run")
#                 #print(c(new.run$value, best.run$value))
#
#                 if(new.run$value<best.run$value) {
#                     if(best.run$value - new.run$value > likelihood.precision) {
#                         times.without.improvement <- 0
#                         if(verbose) {
#                             print("New improvement found, resetting step counter")
#                         }
#                     } else {
#                         if(verbose) {
#                             print("New improvement found, but slight; taking the best value, but not resetting the step counter")
#                         }
#                     }
#                     best.run <- new.run
#                 }
#                 if(verbose) {
#                     step.count <- step.count+1
#                     results.vector<-c(step.count, times.without.improvement, best.run$value, best.run$par)
#                     names(results.vector) <- c("step", "steps.without.improvement","negloglik", names(best.run$par))
#                     print(results.vector)
#                 }
#             }#end of times.without.improvement<10
#             results[[model.index]] <- best.run
#             #try(hessians[[model.index]] <- hessian(func=CalculateLikelihood, x=new.run$par, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)]))
#             results.vector.full <- c(NA, NA, 1, 0, 0)
#             names(results.vector.full) <- GenerateParamLabels()
#             for (i in sequence(length(best.run$par))) {
#                 results.vector.full[which(names(results.vector.full)==names(best.run$par)[i])] <- best.run$par[i]
#             }
#             #print(hessians[[model.index]])
#             #try(print(solve(hessians[[model.index]])))
#             ci.vector<-rep(NA,10)
#             for(parameter in sequence(5)) {
#                 names(ci.vector)[1+2*(parameter-1)] <- paste(GenerateParamLabels()[parameter],"lower", sep=".")
#                 names(ci.vector)[2+2*(parameter-1)] <- paste(GenerateParamLabels()[parameter],"upper", sep=".")
#             }
#             weird.result <- FALSE
#             if(best.run$value>1e100) {
#                 weird.result = TRUE
#                 warning(paste("It seems your likelihood for this run", best.run$value, "is invalid. This probably reflects a problem with numerical optimization for your tree. Your parameter estimates and AIC weights are likely meaningless. We will not calculate confidence: the confidence intervals for all parameters should be taken as c(-Inf, Inf) [unless the parameter is bounded by zero, in which case it is c(0, Inf)]. You could try a transformation of your branch lengths (and make sure to change the time of the flow matrix, too) -- note that parameter values (rate of evolution) will be based on this new timescale. Alternatively, you might allow SE for tip values to be estimated by setting measurement.error=NULL."))
#                 print(paste("It seems your likelihood for this run", best.run$value, "is invalid. This probably reflects a problem with numerical optimization for your tree. Your parameter estimates and AIC weights are likely meaningless. We will not calculate confidence: the confidence intervals for all parameters should be taken as c(-Inf, Inf) [unless the parameter is bounded by zero, in which case it is c(0, Inf)]. You could try a transformation of your branch lengths (and make sure to change the time of the flow matrix, too) -- note that parameter values (rate of evolution) will be based on this new timescale. Alternatively, you might allow SE for tip values to be estimated by setting measurement.error=NULL."))
#             }
#             if(get.se & !weird.result) {
#                 if(verbose) {
#                     print("Now doing simulation to estimate parameter uncertainty")
#                 }
#                 interval.results <- AdaptiveConfidenceIntervalSampling(best.run$par, fn=CalculateLikelihood, lower=lower.bounds[names(preset.starting.parameters)], upper=upper.bounds[names(preset.starting.parameters)], data=data, phy=phy, flow=flow, allow.extrapolation=allow.extrapolation, n.points=n.points,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions.adaptive, allow.restart=allow.restart, best.lnl = best.run$value, likelihood.precision=likelihood.precision, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#                 interval.results.in <- interval.results[which(interval.results[,1]-min(interval.results[,1])<=2),]
#                 interval.results.out <- interval.results[which(interval.results[,1]-min(interval.results[,1])>2),]
#                 if(best.run$value - min(interval.results[,1]) > likelihood.precision) {
#                     print("The sampling to find confidence in parameters actually found a better part of the likelihood surface. Restarting the run for this model at that point")
#                     best.point <- interval.results[which.min(interval.results[,1]),]
#                     preset.starting.parameters <- best.point[-1]
#                     do.run = TRUE
#                 }
#                 if(plot.se) {
#                     pdf(file=paste("Model",models[model.index], "_uncertainty_plot.pdf", sep=""), height=5, width=5*length(preset.starting.parameters))
#                     par(mfcol=c(1, length(preset.starting.parameters)))
#                     for(parameter in sequence(length(preset.starting.parameters))) {
#                         plot(x=interval.results[,parameter+1], y=interval.results[,1], type="n", xlab=names(present.starting.parameters)[parameter], ylab="NegLnL", bty="n", ylim=c(min(interval.results[,1]), min(interval.results[,1])+10))
#                         points(x=interval.results.in[,parameter+1], y=interval.results.in[,1], pch=16, col="black")
#                         points(x=interval.results.out[,parameter+1], y=interval.results.out[,1], pch=16, col="gray")
#                         points(x= best.run$par[parameter], y= best.run$value, pch=1, col="red", cex=1.5)
#                     }
#                     dev.off()
#                     if(verbose) {
#                         print(paste("Uncertainty plot has been saved in Model",models[model.index], "_uncertainty_plot.pdf in ", getwd(), sep=""))
#                     }
#                 }
#                 if(store.sims) {
#                     colnames(interval.results) <- c("neglnL", names(preset.starting.parameters))
#                     all.sims[[model.index]]<-interval.results
#                 }
#                 free.index=0
#                 for(parameter in sequence(length(GenerateParamLabels()))) {
#
#                     if(grepl(GenerateParamLabels()[parameter], names(preset.starting.parameters))) { #is estimated
#                         free.index <- free.index + 1
#                         ci.vector[1+2*(parameter-1)] <- min(interval.results.in[,free.index+1])
#                         ci.vector[2+2*(parameter-1)] <- max(interval.results.in[,free.index+1])
#                     } else {
#                         ci.vector[1+2*(parameter-1)] <- results.vector.full[parameter]
#                         ci.vector[2+2*(parameter-1)] <- results.vector.full[parameter]
#                     }
#                 }
#             }
#             local.df <- data.frame(matrix(c(models[model.index], results.vector.full, AICc(Ntip(phy),k=length(preset.starting.parameters), best.run$value), best.run$value, length(preset.starting.parameters), ci.vector), nrow=1), stringsAsFactors=FALSE)
#             colnames(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K", names(ci.vector))
#             if(do.Higham.correction) {
#                 local.df$penalty=CalculateLikelihood(best.run$par,data=data, phy=phy, flow=flow,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions.adaptive,  likelihood.precision=likelihood.precision, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction, return.penalty=TRUE)
#             }
#             print(local.df)
#             if(!do.run) { #otherwise, we're going to start it again
#                 results.summary <- rbind(results.summary, local.df)
#             }
#         }
#     }
#     results.summary <- cbind(results.summary, deltaAICc=results.summary$AICc-min(results.summary$AICc))
#     results.summary<-cbind(results.summary, AkaikeWeight = AkaikeWeight(results.summary$deltaAICc))
#     if(store.sims) {
#         return(list(results=results.summary, sims=all.sims))
#     }
#     return(results.summary)
# }

# CheckFlow <- function(phy, flow) {
#     lumped.flow <- LumpIntoClades(phy, flow)
#     problem.events <- c()
#     problem.taxa <- c()
#     lumped.flow$problem <- FALSE
#     for (event.index in sequence(nrow(lumped.flow))) {
#         remaining.flow <- lumped.flow[-event.index,]
#         recipient.taxa <- strsplit(lumped.flow$recipient.clades[event.index], ",")[[1]]
#         for (taxon.index in sequence(length(recipient.taxa))) {
#             if(any(grepl(recipient.taxa[taxon.index], remaining.flow$recipient.clades))) {
#                 problem.taxa <- append(problem.taxa, recipient.taxa[taxon.index])
#                 lumped.flow$problem[event.index] <- TRUE
#             }
#         }
#     }
#     problem.taxa <- unique(problem.taxa)
#     return(list(problem.taxa=problem.taxa, all.events=lumped.flow))
# }

# AdjustFlow <- function(data, phy, flow, remove=c("taxa", "events")) {
#     problems <- CheckFlow(phy, flow)
#     number.of.attempts <- 0
#     while(length(problems$problem.taxa)>0) {
#         number.of.attempts <- number.of.attempts + 1
#         if(remove=="taxa") {
#             taxon.to.delete <- sample(problems$problem.taxa, 1)
#             phy <- drop.tip(phy, taxon.to.delete)
#             pruned <- geiger::treedata(phy, data)
#             phy <- pruned$phy
#             data <- pruned$data
#             flow <- subset(flow, flow$recipient!=taxon.to.delete)
#             flow <- subset(flow, flow$donor!=taxon.to.delete)
#         } else {
#             pruned.flow <- problems$all.events[-sample(which(problems$all.events$problem),1),]
#             flow <- UnlumpIntoTaxa(pruned.flow)
#             flow$problem <- NULL
#             # flow$events_names <- as.numeric(as.factor(paste(flow$time.from.root.donor, flow$time.from.root.recipient))) # finally a reason to use factors.
#             # problem.events <- unique(flow$events_names[which(flow$recipient %in% problems)])
#             # problem.event.to.delete <- sample(problem.events, 1)
#             # flow <- flow[-which(flow$events_names==problem.event.to.delete),]
#             # flow$events_names <- NULL
#         }
#         problems <- CheckFlow(phy, flow)
#     }
#     return(list(data=data, phy=phy, flow=flow, number.of.attempts=number.of.attempts))
# }

# BMhybGrid <- function(data, phy, flow, models=c(1,2,3,4), verbose=TRUE, get.se=TRUE, plot.se=TRUE, store.sims=TRUE, precision=2, auto.adjust=FALSE, likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000, measurement.error=NULL, do.kappa.check=FALSE, number.of.proportions=101, number.of.proportions.adaptive=101, allow.restart=TRUE, lower.bounds = c(sigma.sq = 0, mu = -Inf, bt = 1e-06, vh = 0, SE = 0), upper.bounds = c(sigma.sq = 10, mu = Inf,bt = 100,vh = 100,SE = 100), badval.if.not.positive.definite=TRUE, attempt.deletion.fix=FALSE, starting.values=NULL, do.Brissette.correction=FALSE, do.Higham.correction=FALSE, do.DE.correction=FALSE) {
#     flow.problems <- CheckFlow(phy, flow)$problem.taxa
#     if(length(flow.problems)>0) {
#         stop(paste("Sorry, the algorithm cannot work with overlapping hybridization (where any taxon has a history with more than one hybridization event leading to it). In this case, it is multiple events leading to taxon/taxa", paste(flow.problems, collapse=", "), "that are causing the issue. You can edit your flow data.frame manually; you may also use AdjustFlow to randomly delete hybridization events or taxa of hybrid origin."))
#     }
#     if(min(flow$gamma)<0) {
#         stop("Min value of flow is too low; should be between zero and one")
#     }
#     if(max(flow$gamma)>1) {
#         stop("Max value of flow is too high; should be between zero and one")
#     }
#     results<-list()
#     #hessians <- list()
#     results.summary <-data.frame()
#     phy.geiger.friendly <- phy #geiger can't handle branch lengths near zero. Let's lengthen them if needed
#     if(min(phy.geiger.friendly$edge.length)<0.00001) {
#         phy.geiger.friendly$edge.length <- phy.geiger.friendly$edge.length + 0.00001
#     }
#     if(auto.adjust) {
#         phy <- AdjustForDet(phy)
#     }
#     all.sims<-list()
#     if(is.null(starting.values)) {
#         if(verbose) {
#             print("Getting starting values from Geiger")
#         }
#         starting.from.geiger<-NA
#         starting.values <- NA
#         geiger.SE <- data*NA
#         if(!is.null(measurement.error)) {
#             if(length(measurement.error)==1) {
#                 geiger.SE <- rep(measurement.error, length(geiger.SE))
#                 names(geiger.SE) <- phy$tip.label
#             } else {
#                 geiger.SE <- measurement.error
#                 names(geiger.SE) <- phy$tip.label
#
#             }
#             starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
#             starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  starting.from.geiger$sigsq*max(vcv(phy)), mean(measurement.error)) #sigma.sq, mu, beta, vh, SE
#
#         } else {
#             starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=geiger.SE, ncores=1)$opt
#             starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  starting.from.geiger$sigsq*max(vcv(phy)), starting.from.geiger$SE) #sigma.sq, mu, beta, vh, SE
#         }
#         names(starting.values) <- GenerateParamLabels()
#         if(verbose) {
#             print("Done getting starting values")
#         }
#     }
#     #  if(badval.if.not.positive.definite) {
#     if(!IsPositiveDefinite(GetVModified(starting.values, phy, flow, measurement.error=measurement.error))) {
#         if(attempt.deletion.fix) {
#             phy <- AttemptDeletionFix(phy, flow, starting.values)
#             tips <- tips[names(tips) %in% phy$tip.label]
#         }
#         #    stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
#     }
#     #  }
#     for (model.index in sequence(length(models))) {
#         do.run = TRUE
#         preset.starting.parameters = NULL
#         while(do.run) {
#             do.run = FALSE
#
#
#             step.count <- 0
#             if(verbose) {
#                 print(paste("Starting model", models[model.index], "of", length(models), "models"))
#             }
#             model <- models[model.index]
#             if(model==1) {
#                 starting.values <- starting.values[-which(names(starting.values=="bt"))]
#             }
#             if(model==2) {
#                 starting.values <- starting.values[-which(names(starting.values=="vh"))]
#             }
#             if(model==3) {
#                 starting.values <- starting.values[-which(names(starting.values=="bt"))]
#                 starting.values <- starting.values[-which(names(starting.values=="vh"))]
#             }
#             if(!is.null(measurement.error)) {
#                 starting.values <- starting.values[-which(names(starting.values=="SE"))]
#             }
#
#             starting.mins <- c(0, min(data)-(max(data)-min(data)), 0, 0, 0)
#             for(i in sequence(length(starting.mins))) {
#                 starting.mins[i] <- max(starting.mins[i], lower.bounds[i])
#             }
#             starting.maxes <- c(10*starting.values[1], max(data)+(max(data)-min(data)), 10*starting.values[3], 10*starting.values[4], 10*starting.values[5])
#             for (i in sequence(length(starting.maxes))) {
#                 starting.maxes[i] <- min(starting.maxes[i], upper.bounds[i])
#             }
#             names(starting.mins) <- names(preset.starting.parameters)
#             names(starting.maxes) <- names(preset.starting.parameters)
#             ci.vector<-rep(NA,10)
#             for(parameter in sequence(5)) {
#                 names(ci.vector)[1+2*(parameter-1)] <- paste(GenerateParamLabels()[parameter],"lower", sep=".")
#                 names(ci.vector)[2+2*(parameter-1)] <- paste(GenerateParamLabels()[parameter],"upper", sep=".")
#             }
#
#             # if(badval.if.not.positive.definite) {
#             #   if(!IsPositiveDefinite(GetVModified(preset.starting.parameters, phy, flow, actual.params= free.parameters))) {
#             #     stop("It appears your network is in a part of parameter space where calculating likelihood is numerically impossible under a multivariate normal. The best hope is probably removing taxa.")
#             #   }
#             # }
#
#             if(model==1) {
#                 starting.mins["bt"] <- 1
#                 starting.maxes["bt"] <- 1
#             }
#             if(model==2) {
#                 starting.mins["vh"] <- 0
#                 starting.maxes["vh"] <- 0
#             }
#             if(model==3) {
#                 starting.mins["bt"] <- 1
#                 starting.maxes["bt"] <- 1
#                 starting.mins["vh"] <- 0
#                 starting.maxes["vh"] <- 0
#             }
#             if(!is.null(measurement.error)) {
#                 starting.mins["SE"] <- 0
#                 starting.maxes["SE"] <- 0
#             }
#
#
#             grid.of.points <- lhs::randomLHS(n=n.points, k=length(starting.mins))
#             for(parameter.index in sequence(ncol(grid.of.points))) {
#                 grid.of.points[,parameter.index] <- starting.mins[parameter.index] + grid.of.points[,parameter.index] * (starting.maxes[parameter.index] - starting.mins[parameter.index])
#             }
#             colnames(grid.of.points) <- names(starting.parameters)
#             likelihoods <- rep(NA, n.points)
#
#             for (rep.index in sequence(n.points)) {
#                 local.likelihood <- try(CalculateLikelihood(as.numeric(grid.of.points[rep.index,]), data=data, phy=phy, flow=flow, measurement.error=measurement.error, badval.if.not.positive.definite=badval.if.not.positive.definite, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction))
#                 if(!is.numeric(local.likelihood)) {
#                     local.likelihood <- (0.5)*.Machine$double.xmax
#                 }
#                 likelihoods[rep.index] <- local.likelihood
#                 if(verbose & rep.index%%50==0) {
#                     print(paste0("Now done ", rep.index, " of ", n.points, " to analyze (", round(100*rep.index/n.points, 4), "% done)"))
#                 }
#             }
#
#             best.one <- which.min(likelihoods)[1]
#             best.params <- grid.of.points[best.one,]
#             best.likelihood <- min(likelihoods, na.rm=TRUE)
#
#
#             results.vector.full <- c(NA, NA, 1, 0, 0)
#             names(results.vector.full) <- GenerateParamLabels()
#             for (i in sequence(length(best.params))) {
#                 results.vector.full[which(names(results.vector.full)==names(best.params)[i])] <- best.params[i]
#             }
#
#         }
#
#         if(get.se) {
#             if(verbose) {
#                 print("Now doing simulation to estimate parameter uncertainty")
#             }
#             previous.results <- cbind(likelihoods, grid.of.points)
#             colnames(previous.results)[1] <- "NegLogL"
#             interval.results <- AdaptiveConfidenceIntervalSampling(best.params[names(preset.starting.parameters)], fn=CalculateLikelihood, lower=lower.bounds[names(preset.starting.parameters)], upper=upper.bounds[names(preset.starting.parameters)], data=data, phy=phy, flow=flow, allow.extrapolation=allow.extrapolation, n.points=n.points,  measurement.error=measurement.error, do.kappa.check=do.kappa.check, number.of.proportions=number.of.proportions.adaptive, allow.restart=allow.restart, best.lnl = min(likelihoods), likelihood.precision=likelihood.precision, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction)
#             colnames(interval.results) <- c("NegLogL", GenerateParamLabels())
#             if(!any(grepl("bt", names(interval.results)))) {
#                 interval.results$bt <- 1
#             }
#             if(!any(grepl("vh", names(interval.results)))) {
#                 interval.results$vh <- 0
#             }
#             if(!any(grepl("SE", names(interval.results)))) {
#                 interval.results$SE <- 0
#             }
#             interval.results <- interval.results[, colnames(previous.results)]
#             interval.results.local <- interval.results
#             interval.results <- rbind(previous.results, interval.results)
#             interval.results <- interval.results[is.finite(interval.results[,1]),]
#             interval.results.in <- interval.results[which(interval.results[,1]-min(as.numeric(interval.results[,1]))<=2),]
#             interval.results.out <- interval.results[which(interval.results[,1]-min(as.numeric(interval.results[,1]))>2),]
#
#             for(parameter in sequence(ncol(interval.results)-1)) {
#                 parameter.name <- names(interval.results)[parameter+1]
#                 ci.vector[paste0(parameter.name, ".upper")] <- max(interval.results.in[,parameter+1], na.rm=TRUE)
#                 ci.vector[paste0(parameter.name, ".lower")] <- min(interval.results.in[,parameter+1], na.rm=TRUE)
#             }
#             if(min(interval.results$NegLogL, na.rm=TRUE) <  min(likelihoods)) {
#                 best.params <- interval.results[which.min(interval.results$NegLogL), -1]
#                 results.vector.full <- c(NA, NA, 1, 0, 0)
#                 names(results.vector.full) <- GenerateParamLabels()
#                 for (i in sequence(length(best.params))) {
#                     results.vector.full[which(names(results.vector.full)==names(best.params)[i])] <- best.params[i]
#                 }
#                 best.likelihood <- min(interval.results$NegLogL, na.rm=TRUE)
#             }
#             if(plot.se) {
#                 pdf(file=paste("Model",models[model.index], "_uncertainty_plot.pdf", sep=""), height=5, width=5*length(starting.parameters))
#                 par(mfcol=c(1, length(starting.parameters)))
#                 for(parameter in sequence(length(starting.parameters))) {
#                     plot(x=interval.results[,parameter+1], y=interval.results[,1], type="n", xlab=names(starting.parameters)[parameter], ylab="NegLnL", bty="n", ylim=c(min(interval.results[,1]), min(interval.results[,1])+10))
#                     points(x=interval.results.in[,parameter+1], y=interval.results.in[,1], pch=16, col="black")
#                     points(x=interval.results.out[,parameter+1], y=interval.results.out[,1], pch=16, col="gray")
#                     points(x= best.params[parameter], y= best.likelihood, pch=1, col="red", cex=1.5)
#                 }
#                 dev.off()
#                 if(verbose) {
#                     print(paste("Uncertainty plot has been saved in Model",models[model.index], "_uncertainty_plot.pdf in ", getwd(), sep=""))
#                 }
#             }
#         }
#         local.df <- data.frame(matrix(c(models[model.index], results.vector.full, AICc(Ntip(phy),k=length(starting.parameters), best.likelihood), best.likelihood, length(starting.parameters), ci.vector), nrow=1), stringsAsFactors=FALSE)
#         local.df <- apply(local.df, 2, unlist)
#         names(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K", names(ci.vector))
#
#
#         if(do.Higham.correction) {
#             param.estimates <- unlist(results.vector.full)
#             names(param.estimates) <- names(preset.starting.parameters)
#             local.df$penalty=CalculateLikelihood(param.estimates,data=data, phy=phy, flow=flow,  measurement.error=measurement.error, do.kappa.check=do.kappa.check,  number.of.proportions=number.of.proportions.adaptive,  likelihood.precision=likelihood.precision, lower.b=lower.bounds[names(preset.starting.parameters)], upper.b=upper.bounds[names(preset.starting.parameters)], restart.mode=TRUE, do.Brissette.correction=do.Brissette.correction, do.Higham.correction=do.Higham.correction, do.DE.correction=do.DE.correction, return.penalty=TRUE)
#         }
#         all.points <- data.frame(grid.of.points, stringsAsFactors=FALSE)
#         all.points$NegLogL <- likelihoods
#         all.points$Model <- models[model.index]
#         all.points$AICc <- AICc(Ntip(phy),k=length(starting.parameters), all.points$NegLogL)
#         all.points$K <- length(starting.parameters)
#         results.summary <- rbind(results.summary, data.frame(t(local.df), stringsAsFactors=FALSE))
#         all.sims <- rbind(all.sims, all.points)
#
#     }
#     results.summary <- cbind(results.summary, deltaAICc=as.numeric(results.summary$AICc)-min(as.numeric(results.summary$AICc)))
#     results.summary<-cbind(results.summary, AkaikeWeight = AkaikeWeight(results.summary$deltaAICc))
#     #save(list=ls(), file="~/Desktop/everything.rda")
#     if(store.sims) {
#         all.sims$deltaAICc <- all.sims$AICc - min(as.numeric(all.sims$AICc))
#         all.sims$AkaikeWeight <- AkaikeWeight(all.sims$deltaAICc)
#         return(list(results=results.summary, sims=all.sims))
#     }
#     return(results.summary)
# }
#
# PlotAICRegion <- function(sims, show.finite.only=TRUE, true.params=NULL, ...) {
#     pairs.of.params <- utils::combn(colnames(sims)[1:5],2)
#     sims.to.plot <- sims
#     if(show.finite.only) {
#         sims.to.plot <- sims.to.plot[which(sims.to.plot$NegLogL<1e300),]
#     }
#     best.one <- which.min(sims.to.plot$AICc)[1]
#     par(mfcol=c(2, ceiling(ncol(pairs.of.params)/2)))
#     for(pair.index in sequence(ncol(pairs.of.params))) {
#         plot(sims.to.plot[,pairs.of.params[1,pair.index]], sims.to.plot[,pairs.of.params[2,pair.index]], pch=".", col=rgb(0,0,0,.5), xlim=range(sims[,pairs.of.params[1,pair.index]]), ylim=range(sims[,pairs.of.params[2,pair.index]]), xlab=paste0(pairs.of.params[1,pair.index],ifelse(0==max(sims[,pairs.of.params[1,pair.index]])-min(sims[,pairs.of.params[1,pair.index]]), " FIXED", "")), ylab=paste0(pairs.of.params[2,pair.index],ifelse(0==max(sims[,pairs.of.params[2,pair.index]])-min(sims[,pairs.of.params[2,pair.index]]), " FIXED", "")), ...)
#         if(!is.null(true.params)) {
#             points(x=true.params[pairs.of.params[1,pair.index]], y=true.params[pairs.of.params[2,pair.index]], pch=8, col="blue")
#         }
#         points(x=sims.to.plot[best.one,pairs.of.params[1,pair.index]], y=sims.to.plot[best.one,pairs.of.params[2,pair.index]], pch=5, col="red")
#     }
# }



# DetPass <- function(phy) {
#     det.pass <- TRUE
#     vcv.result <- vcv.phylo(phy)
#     det.tries <- c(det(vcv.result), det(1000*vcv.result), det(0.0001*vcv.result))
#     if(min(det.tries)<0) {
#         det.pass <- FALSE
#     }
#     if(sum(is.finite(det.tries))!=length(det.tries)) {
#         det.pass <- FALSE
#     }
#     return(det.pass)
# }

# AdjustForDet <- function(phy, max.attempts=100) {
#     attempts<-0
#     if(!DetPass(phy)) {
#         warning("Determininant of the phylogeny was difficulty to calculate, so the phylogeny needed to be adjusted. Your results may be approximate as a result")
#         while(!DetPass(phy) && attempts<=max.attempts) {
#             attempts <- attempts+1
#             phy$edge.length <- phy$edge.length+0.00001*attempts
#         }
#     }
#     if(attempts>max.attempts) {
#         stop("Phylogeny could not be adjusted adequately")
#     }
#     return(phy)
# }
#
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

# GetVModified <- function(x,phy,flow,measurement.error){
#     bt <- 1
#     vh <- 0
#     sigma.sq <-x['sigma.sq']
#     mu <- x['mu']
#     SE <- 0
#
#     if(grepl('SE', names(x))) {
#         SE <- x['SE']
#     }
#
#     if(grepl('bt', names(x))) {
#         bt <- x['bt']
#     }
#
#     if(grepl('vh', names(x))) {
#         vh <- x['vh']
#     }
#
#
#     flow <- flow[order(flow$time.from.root.donor),]
#     crossing.flows <- any(!((order(flow$time.from.root.donor) == order(flow$time.from.root.recipient))))
#     times.original<-vcv.phylo(phy,model="Brownian")
#     V.original <- sigma.sq * times.original
#     V.modified <- V.original
#     table.donor <- table(flow$donor)
#     #print(table.donor)
#     for(table.index in 1:length(table.donor)){
#         sub.flow<-flow[flow$donor == names(table.donor)[table.index],]
#         #print(sub.flow)
#         sub.flow<-sub.flow[order(sub.flow$time.from.root.donor),]
#         #print(table.donor[table.index]==1)
#         if(table.donor[table.index]==1){ # one donor to one recipient
#             gamma <- sub.flow$gamma
#             recipient.index <- which(rownames(V.modified) == sub.flow$recipient)
#             donor.index <- which(rownames(V.modified) == sub.flow$donor)
#             V.modified[recipient.index,donor.index] <- gamma*sub.flow$time.from.root.donor*sigma.sq
#             V.modified[donor.index,recipient.index] <- V.modified[recipient.index,donor.index]
#             V.modified[recipient.index, recipient.index] <- (V.original[recipient.index, recipient.index] -  sigma.sq*sub.flow$time.from.root.recipient) + (gamma^2 + (1- gamma)^2) * (sub.flow$time.from.root.recipient)*sigma.sq + 2*gamma*(1-gamma)*V.original[recipient.index, donor.index] + vh
#
#         }else{ # one donor to several recipients
#             #print(table.donor[table.index])
#             gamma <- sub.flow$gamma[1]
#             recipient.index.set <- which(rownames(V.modified) %in% sub.flow$recipient)
#             donor.index <- which(rownames(V.modified) == sub.flow$donor[1]) #one donor here
#             for(index in 1:length(recipient.index.set)){ #so only use recipient index to match the poisition
#                 recipient.index<-recipient.index.set[index]
#                 donor.catch.right.recipient <-which(sub.flow$recipient==rownames(V.modified)[recipient.index])
#                 V.modified[recipient.index,donor.index] <- gamma*sub.flow$time.from.root.donor[donor.catch.right.recipient]*sigma.sq
#                 V.modified[donor.index,recipient.index] <- V.modified[recipient.index,donor.index]
#                 V.modified[recipient.index, recipient.index] <- (V.original[recipient.index, recipient.index] -  sigma.sq*sub.flow$time.from.root.recipient[donor.catch.right.recipient]) + (gamma^2 + (1- gamma)^2) * (sub.flow$time.from.root.recipient[donor.catch.right.recipient])*sigma.sq + 2*gamma*(1-gamma)*V.original[recipient.index, donor.index] + vh
#             }
#
#             for(index.1 in 1:dim(sub.flow)[1]){
#                 recipient.index.1 <- recipient.index.set[index.1]
#                 for(index.2 in 1:dim(sub.flow)[1]){
#                     recipient.index.2 <- recipient.index.set[index.2]
#                     if(recipient.index.1!=recipient.index.2){
#                         V.modified[recipient.index.1,recipient.index.2] = (V.original[recipient.index.1,recipient.index.2]- sub.flow$time.from.root.recipient[1])*sigma.sq + (gamma^2+(1-gamma)^2)*sub.flow$time.from.root.recipient[1]*sigma.sq
#                     }#end of flow.Index.1 != flow.Index.2
#                 }#end flow.Index.2
#             }#end flow.Index.1
#         }#end of else
#     }#end of for table.index
#
#     if(is.null(measurement.error)){
#         diag(V.modified)<-diag(V.modified)+ SE^2
#     }else{
#         diag(V.modified)<-diag(V.modified)+ measurement.error^2
#     }
#     return(V.modified)
# }

# GetMeansModified <- function(x, phy, flow) {
#     badval<-(0.5)*.Machine$double.xmax
#     bt <- 1
#     vh <- 0
#     sigma.sq <-x['sigma.sq']
#     mu <- x['mu']
#     SE <- 0
#
#     if(grepl('SE', names(x))) {
#         SE <- x['SE']
#     }
#
#     if(grepl('bt', names(x))) {
#         bt <- x['bt']
#     }
#
#     if(grepl('vh', names(x))) {
#         vh <- x['vh']
#     }
#
#     times.original <-vcv.phylo(phy, model="Brownian") #is the initial one based on the tree alone, so just time
#     V.original <- sigma.sq * times.original
#
#     means.original <- rep(mu, Ntip(phy))
#     names(means.original) <- rownames(V.original)
#     means.modified <- means.original
#
#     means.original <- rep(mu, Ntip(phy))
#     names(means.original) <- rownames(V.original)
#     means.modified <- means.original
#     for (flow.index in sequence(dim(flow)[1])) {
#         recipient.index <- which(rownames(V.original)==flow$recipient[flow.index])
#         if(length(recipient.index)!=1) {
#             stop(paste("Tried to find ", flow$recipient[flow.index], " but instead found ", paste(rownames(V.original)[recipient.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
#         }
#         donor.index <- which(rownames(V.original)==flow$donor[flow.index])
#         if(length(donor.index)!=1) {
#             stop(paste("Tried to find ", flow$donor[flow.index], " but instead found ", paste(rownames(V.original)[donor.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", sep=""))
#         }
#         means.modified[recipient.index] <- means.original[recipient.index] + log(bt)
#     }
#     return(means.modified)
# }

#precision is the cutoff at which we think the estimates become unreliable due to ill conditioned matrix
# CalculateLikelihood <- function(x, data, phy, flow, precision=2, proportion.mix.with.diag=0, allow.extrapolation=FALSE, measurement.error, do.kappa.check=FALSE, number.of.proportions=101, lower.b=c(0, -Inf, 0.000001, 0, 0), upper.b=c(10,Inf,100,100,100), badval.if.not.positive.definite=TRUE, do.Brissette.correction=FALSE, do.Higham.correction=FALSE, do.DE.correction=FALSE, return.penalty=FALSE, ...) {
#     badval<-(0.5)*.Machine$double.xmax
#     bt <- 1
#     vh <- 0
#     sigma.sq <-x['sigma.sq']
#     mu <- x['mu']
#     SE <- 0
#
#     if(grepl('SE', names(x))) {
#         SE <- x['SE']
#     }
#
#     if(grepl('bt', names(x))) {
#         bt <- x['bt']
#     }
#
#     if(grepl('vh', names(x))) {
#         vh <- x['vh']
#     }
#     if(any(x<lower.b) | any(x>upper.b)) {
#         return(badval)
#     }
#     V.modified <- GetVModified(x, phy, flow, measurement.error=measurement.error)
#     if(do.Brissette.correction) {
#         V.modified <- BrissetteEtAlCorrection(V.modified)
#         if(is.null(V.modified)) {
#             return(badval)
#         }
#     }
#     if(do.Higham.correction & !IsPositiveDefinite(V.modified)) {
#         new.mat <- as.matrix(Matrix::nearPD(V.modified, corr=FALSE, posd.tol = 1e-16, eig.tol = 1e-16, conv.tol = 1e-16)$mat)
#         if(any(new.mat!=V.modified)) {
#             warning("Had to do Higham (2002) correction for not positive definite matrix")
#             likelihood.penalty <- 10+dist(rbind(c(new.mat), c(V.modified)))
#         }
#         V.modified <- new.mat
#         #  if(min(V.modified)<0) {
#         #    warning("Had to remove negative values in VCV after Higham (2002) correction")
#         #    V.modified[(V.modified<0)] <- 0
#         #  }
#     }
#     if(return.penalty) {
#         return(likelihood.penalty)
#     }
#     if(do.DE.correction & !IsPositiveDefinite(V.modified)) {
#         warning("Have to modify variance covariance matrix to make it positive definite, so results are approximate and the analysis will be slow.")
#         V.modified <-  AlterMatrixUsingDE(V.modified)
#     }
#     if(badval.if.not.positive.definite) {
#         if(!IsPositiveDefinite(V.modified)) {
#             return(badval)
#         }
#     }
#     means.modified <- GetMeansModified(x, phy, flow)
#     if(sigma.sq <0 || vh<0 || bt <= 0.0000001 || SE < 0) {
#         return(badval)
#     }
#     #if(min(eigen(V.modified)$values) < 0.00000001) {
#     #  return(0.5*badval) #b/c it might not be very bad, just hard to calculate
#     #}
#     data <- data[match(names(means.modified), names(data))]
#     if(length(data)!=length(means.modified)) {
#         stop("Mismatch between names of taxa in data vector and on phy")
#     }
#     #if(length(which(eigen(V.modified)$values<0))>0) {
#     #	last.bad <- V.modified
#     #	return(badval)
#     #}
#     #NegLogML <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified)))
#     NegLogML <- NULL
#     try(NegLogML <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified) + (1/2)*determinant(V.modified, logarithm=TRUE)$modulus + likelihood.penalty, silent=TRUE)
#     if(is.null(NegLogML)) {
#         NegLogML <- badval
#     }
#     #print(paste0("NegLogML = ", NegLogML-likelihood.penalty, ", penalty=", likelihood.penalty))
#     # NegLogML.dmvnorm <- -dmvnorm(x=data, mean=means.modified, sigma=V.modified, log=TRUE)
#     # if(!is.finite(NegLogML.dmvnorm)) {
#     #
#     #   print("sigma.sq, vh, bt, SE, mu")
#     #   print(c(sigma.sq, vh, bt, SE, mu))
#     #   print("us vs dmvnorm")
#     #   print(c(NegLogML, NegLogML.dmvnorm))
#     # }
#     #if(min(V.modified)<0 || sigma.sq <0 || vh<0 || bt <= 0.0000001 || !is.finite(NegLogML) || SE<0) {
#     if( sigma.sq <0 || vh<0 || bt <= 0.0000001 || !is.finite(NegLogML) || SE<0) {
#         NegLogML<-badval
#         print("badval")
#         print(paste0("min(V.modified) ",min(V.modified)))
#         print(paste0("sigma.sq ", sigma.sq))
#         print(paste0("vh ", vh))
#         print(paste0("bt ", bt))
#     }
#     matrix.condition <- kappa(V.modified, exact=TRUE)
#     #print("condition")
#     #print(kappa(V.modified, exact=TRUE))
#     #print("log(condition)")
#     #print(log(kappa(V.modified, exact=TRUE)))
#
#     #pretty<-c(NegLogML, log(matrix.condition))
#     #names(pretty) <- c("NegLogL", "log(matrix.condition")
#     #print(pretty)
#     #The ratio  of the largest to smallest singular value in the singular value decomposition of a matrix. The base- logarithm of  is an estimate of how many base- digits are lost in solving a linear system with that matrix. In other words, it
#     #estimates worst-case loss of precision. A system is said to be singular if the condition number is infinite, and ill-conditioned if it is too large, where "too large" means roughly  the precision of matrix entries.
#     #if(rcond(V.modified) < .Machine$double.eps^2){
#     if(log(matrix.condition) > precision & do.kappa.check) {
#         proportions <- seq(from=1, to=0, length.out=number.of.proportions)
#         lnl.vector <- rep(NA, length(proportions))
#         max.diff <- 0
#         kappa.vector <- rep(NA, length(proportions))
#         for(i in sequence(length(proportions))) {
#             V.modified.by.proportions<-(1-proportions[i]) * V.modified + proportions[i] * diag(dim(V.modified)[1]) * diag(V.modified)
#             #local.lnl <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified.by.proportions)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified.by.proportions)))
#             local.lnl <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified.by.proportions)%*%(data-means.modified) + (1/2)*determinant(V.modified.by.proportions, logarithm=TRUE)$modulus
#             lnl.vector[i] <- local.lnl
#             kappa.vector[i] <- kappa(V.modified.by.proportions, exact=TRUE)
#             # if(i>6) {
#             # 	very.local.lnl <- lnl.vector[(i-6):(i-1)]
#             # 	max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)])) #looking locally for jumps in the likelihood
#             # 	current.diff <- abs(local.lnl - lnl.vector[i-1])
#             # 	if(current.diff > 2 * max.diff) {
#             # 		#print(paste("breaking after ", i))
#             # 		#break() #the modified matrix is still poorly conditioned, so stop here
#             #     lnl.vector[i] <- NA
#             #     break()
#             # 	}
#             # }
#         }
#         # small.kappas <- which(kappa.vector < 2)
#         # if(length(small.kappas)<10) {
#         #   small.kappas<-order(kappa.vector, decreasing=FALSE)[1:10]
#         # }
#         #smoothing.estimate <- smooth.spline(x=proportions[small.kappas], y=lnl.vector[small.kappas], w=1/kappa.vector[small.kappas], df=4)$spar
#         proportions<-proportions[which(!is.na(lnl.vector))]
#         lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
#         NegLogML.4 <- predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=4), data.frame(proportions =0.000))$y #df of 4: don't want overfitting
#         NegLogML.6 <- predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=6), data.frame(proportions =0.000))$y #df of 4: don't want overfitting
#         if(abs(NegLogML.4 - NegLogML.6)>0.1) { #the result is very sensitive to the degrees of freedom.
#             NegLogML <- badval
#         } else {
#             NegLogML <- NegLogML.4
#         }
#
#         #plot(proportions, lnl.vector)
#         #lines(proportions, predict(smooth.spline(x=proportions, y=lnl.vector, w=1/kappa.vector, df=4), data.frame(proportions))$y[,1])
#         #plot(c(0, proportions), c(NegLogML, lnl.vector), type="n")
#         #points(proportions, lnl.vector, pch=20)
#         #points(0, NegLogML, col="red")
#         if(abs(NegLogML - lnl.vector[length(lnl.vector)]) > 0.001) {	#means this point was extrapolated b/c the likelihood surface got strange
#             if(allow.extrapolation) {
#                 warning("VCV matrix was ill-conditioned, so used splines to estimate its likelihood (allow.extrapolation=TRUE). This could lead to very bad estimates of the likelihood")
#             } else {
#                 NegLogML <- badval
#                 warning("VCV matrix was ill-conditioned at some points in parameter space; treated the likelihood at these points as a bad value rather than estimating it. If you think the MLE is in this region, you could try allow.extrapolation=TRUE")
#             }
#         }
#
#         #print(paste("Did interpolation, got ", NegLogML))
#     }
#
#
#     #print("datadiff")
#     #print(quantile(data-means.modified))
#     #print("middle")
#     #print((1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified))
#     #print("end")
#     #print((1/2)*log(abs(det(V.modified))) )
#     #print(x)
#     #print(V.modified[1:10,1:10])
#     #print(means.modified)
#     # if(NegLogML < -1000 | NegLogML > 100000) {
#     #   stop("bad")
#     # } else {
#     #   stop("good")
#     # }
#     #if(NegLogML< (0)) {
#     #	NegLogML <- badval #since something seems off.
#     #}
#
#     return(NegLogML[1])
#
# }

ConvertLog1P <- function(x) {
    x[which(names(x)=="mu")] <- log1p(x[which(names(x)=="mu")])
}

ConvertExpm1 <- function(x) {
    x[which(names(x)=="mu")] <- expm1(x[which(names(x)=="mu")])
}

# AdaptiveConfidenceIntervalSampling <- function(par, fn, lower=-Inf, upper=Inf, desired.delta = 2, n.points=5000, verbose=TRUE, measurement.error, do.kappa.check=FALSE, allow.restart=TRUE,  best.lnl = -Inf, likelihood.precision=0.001, restart.mode=FALSE, ...) {
#     starting<-fn(par, measurement.error=measurement.error,  ...)
#     if(length(lower) < length(par)) {
#         lower<-rep(lower, length(par))
#     }
#     if(length(upper) < length(par)) {
#         upper<-rep(upper, length(par))
#     }
#     min.multipliers <- rep(1, length(par))
#     max.multipliers <- rep(1, length(par))
#     results<-data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
#     results[1,]<-unname(c(starting, par))
#     for (i in sequence(n.points)) {
#         sim.points<-NA
#         while(is.na(sim.points[1]) | !is.numeric(sim.points[1])) {
#             sim.points<-GenerateValues(par, lower, upper, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
#         }
#         results[i+1,] <- c(fn(sim.points, measurement.error=measurement.error, do.kappa.check=do.kappa.check, ...), sim.points)
#         if(i>5 & restart.mode) {
#             if((best.lnl - min(results[,1], na.rm=TRUE) > likelihood.precision ) & allow.restart) {
#                 results <- results[sequence(i+1),] #stop here and restart
#                 return(results)
#             }
#         }
#         if (i%%20==0) {
#             for (j in sequence(length(par))) {
#                 returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
#                 total.range <- range(results[,j+1], na.rm=TRUE)
#                 width.ratio <- diff(returned.range)/diff(total.range)
#                 if(is.na(width.ratio)) {
#                     width.ratio=1
#                 }
#                 if(width.ratio > 0.5) { #we are not sampling widely enough
#                     min.multipliers[j] <- min.multipliers[j] * 0.9
#                     max.multipliers[j] <- max.multipliers[j] * 1.1 #expand the range
#                 } else {
#                     min.multipliers[j] <- 1
#                     max.multipliers[j] <- 1
#                 }
#             }
#         }
#         if (verbose && i%%100==0) {
#             print(paste(i, "of", n.points, "done"))
#         }
#     }
#     return(results)
# }

ComputeConfidenceIntervals <- function(par, fn, traits, desired.delta = 2, n.points=5000, verbose=TRUE, do.kappa.check=FALSE,  best.lnl = -Inf, confidence.points=5000, likelihood.precision=0.01, ...) {
    starting<-fn(parameters=par, traits=traits, ...)
    lower <- c(sigma.sq=0, mu=min(c(10*min(traits), .1*min(traits))), bt=-10, vh=0,SE=0)[names(par)] #this is why we pass traits explicitly
    upper <- c(sigma.sq=10, mu=max(c(10*max(traits), .1*max(traits))), bt=10, vh=10,SE=10)[names(par)]
    for (i in seq_along(par)) {
        lower[i] <- min(lower[i], par[i])
        upper[i] <- max(upper[i], par[i])
    }
    min.multipliers <- rep(1, length(par))
    max.multipliers <- rep(1, length(par))
    results<-data.frame(data.frame(matrix(nrow=confidence.points+1, ncol=1+length(par))))
    names(results) <- c("negloglik", names(par))
    results[1,]<-(c(starting, par))
    for (i in sequence(confidence.points)) {
        sim.points<-NA
        while(is.na(sim.points[1]) | !is.numeric(sim.points[1])) {
            sim.points<-GenerateValues(par, lower, upper, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
        }
        results[i+1,] <- c(fn(sim.points, traits=traits, ...), sim.points)
        # if(i>5 & restart.mode) {
        #     if((best.lnl - min(results[,1], na.rm=TRUE) > likelihood.precision ) & allow.restart) {
        #         results <- results[sequence(i+1),] #stop here and restart
        #         return(results)
        #     }
        # }
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
            print(paste(i, "of", confidence.points, "done"))
        }
    }
    return(results)
}

# GenerateRandomValues <- function(data, parameters, lower, upper) {
#     new.vals <- c(-Inf, Inf)
#     while(any(new.vals < lower) | any(new.vals>upper)) {
#         new.vals <- c(stats::rexp(1,10), stats::runif(1, min=min(data), max=max(data)), stats::rexp(1,1), stats::rexp(1,10), stats::rexp(1,10))
#         names(new.vals) <- GenerateParamLabels()
#         new.vals <- new.vals[names(parameters)]
#     }
#     return(new.vals)
# }

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
            new.vals[i]<-stats::runif(1, min=ifelse(is.finite(new.vals.bounds[1]),new.vals.bounds[1], 0.000001) , max=ifelse(is.finite(new.vals.bounds[2]), new.vals.bounds[2], 10000))

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
    names(new.vals) <- names(par)
    return(new.vals)
}

# GetClade <- function(phy, clade.size) {
#     nodes <- phy$edge[,1]
#     subtrees <- lapply(nodes, ape::extract.clade, phy=phy)
#     counts <- sapply(subtrees, ape::Ntip)
#     matches<-subtrees[which(counts==clade.size)]
#     if(length(matches)==0) {
#         return(NA)
#     }
#     lucky <- matches[sample.int(length(matches),1)][[1]]
#     return(phytools::findMRCA(phy, tips=lucky$tip.label, type="node"))
# }
#
# GetAncestor <- function(phy, node) {
#     return(phy$edge[which(phy$edge[,2]==node),1])
# }


# AttemptDeletionFix <- function(phy, flow, params=c(1,0,0.1, 0, 0), m.vector = c(1,2)) {
#     if(is.null(names(params))) {
#         names(params) <- GenerateParamLabels()
#     }
#     taxa.to.try.deleting <- phy$tip.label
#     taxa.to.try.deleting <- taxa.to.try.deleting[!(taxa.to.try.deleting %in% flow$recipient)]
#     taxa.to.try.deleting <- taxa.to.try.deleting[!(taxa.to.try.deleting %in% flow$donor)]
#     taxa.to.try.deleting <- taxa.to.try.deleting[sample.int(length(taxa.to.try.deleting), length(taxa.to.try.deleting), replace=FALSE)]
#     if(length(taxa.to.try.deleting)==0) {
#         stop("There are no taxa to delete that aren't involved in hybridization.")
#     }
#     phy.pruned <- phy
#     current.m.index <- 1
#     current.index <- 1
#     combos.to.delete <- utils::combn(taxa.to.try.deleting,m.vector[current.m.index])
#     placeholder.params <- rep(TRUE, 5)
#     names(placeholder.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#     while(!IsPositiveDefinite(GetVModified(params, phy.pruned, flow, measurement.error=0))) {
#         #print(current.index)
#         #print(paste0("Dropping ", paste(combos.to.delete[,current.index], collapse=" ")))
#         phy.pruned <- ape::drop.tip(phy, combos.to.delete[,current.index])
#         current.index <- current.index + 1
#         if(current.index > ncol(combos.to.delete)) {
#             if(current.m.index < length(m.vector)) {
#                 current.m.index <- current.m.index + 1
#                 combos.to.delete <- utils::combn(taxa.to.try.deleting,m.vector[current.m.index])
#                 current.index <- 1
#             } else {
#                 stop(paste0("Correction by removing random taxa up to ", max(m.vector), " taxa at a time failed."))
#             }
#         }
#     }
#     return(phy.pruned)
# }

AddNodeToPhygraph <- function(below.node, depth.below,  phy.graph, tip.label, terminal.length=0) {
    new.phy.graph <- phy.graph
    new.tip.number <- ape::Ntip(phy.graph) + 1 #need to then offset everything else
    new.internal.node.number <- ape::Nnode(phy.graph) + ape::Ntip(phy.graph) + 2
    new.phy.graph$reticulation[which(new.phy.graph$reticulation>ape::Ntip(phy.graph))] <- 1 + new.phy.graph$reticulation[which(new.phy.graph$reticulation>ape::Ntip(phy.graph))]
    new.phy.graph$edge[which(new.phy.graph$edge>ape::Ntip(phy.graph))] <- 1 + new.phy.graph$edge[which(new.phy.graph$edge>ape::Ntip(phy.graph))]
    edge.to.break <- which(phy.graph$edge[,2]==below.node)
    parent.node <- new.phy.graph$edge[edge.to.break,1]
    new.below.node <- new.phy.graph$edge[edge.to.break,2]
    original.edge.length <- phy.graph$edge.length[edge.to.break]
    if(round(depth.below,4) > round(original.edge.length,4)) {
        stop("Node is added below the beginning of the edge (depth.below too large)")
    }

    #delete original edge
    new.phy.graph$edge <- new.phy.graph$edge[-edge.to.break,]
    new.phy.graph$edge.length <- new.phy.graph$edge.length[-edge.to.break]

    #add lower edge
    new.phy.graph$edge <- rbind(new.phy.graph$edge, c(parent.node, new.internal.node.number))
    new.phy.graph$edge.length[1+length(new.phy.graph$edge.length)] <- original.edge.length - depth.below

    #add upper edge
    new.phy.graph$edge <- rbind(new.phy.graph$edge, c(new.internal.node.number, new.below.node))
    new.phy.graph$edge.length[1+length(new.phy.graph$edge.length)] <- depth.below

    #add the new taxon
    new.phy.graph$edge <- rbind(new.phy.graph$edge, c(new.internal.node.number, new.tip.number))
    new.phy.graph$edge.length[1+length(new.phy.graph$edge.length)] <- terminal.length
    new.phy.graph$tip.label <- c(new.phy.graph$tip.label,tip.label)

    new.phy.graph$Nnode <- 1 + new.phy.graph$Nnode



    return(new.phy.graph)
}

RenumberPhygraph <- function(phy.graph) {
    new.phy.graph <- phy.graph
    heightnode <- function(node, phy) {
        return(phytools::nodeheight(phy, node))
    }
    heights <- sapply(sequence(ape::Nnode(phy.graph) + ape::Ntip(phy.graph)), heightnode, phy=phy.graph) # does not line up with edge matrix or edge lengths, remember
    old2new <- data.frame(old=sequence(ape::Nnode(phy.graph) + ape::Ntip(phy.graph)), new=NA, height=heights)
    old2new <- subset(old2new, old2new$old > ape::Ntip(phy.graph))
    old2new <- old2new[order(old2new$height),]
    old2new$new <- seq(from=ape::Ntip(phy.graph)+1, to=ape::Nnode(phy.graph) + ape::Ntip(phy.graph), by=1)
    for (i in sequence(nrow(old2new))) {
        new.phy.graph$edge[which(phy.graph$edge==old2new$old[i])] <- old2new$new[i]
        new.phy.graph$reticulation[which(phy.graph$reticulation==old2new$old[i])] <- old2new$new[i]
    }
    return(new.phy.graph)
}

PruneDonorsRecipientsFromVCV <- function(VCV) {
  VCV <- VCV[!grepl("donor_", rownames(VCV)), !grepl("donor_", colnames(VCV))]
  VCV <- VCV[!grepl("recipient_", rownames(VCV)), !grepl("recipient_", colnames(VCV))]
  return(VCV)
}

PruneDonorsRecipientsFromMeans <- function(x) {
  x <- x[!grepl("donor_", names(x))]
  x <- x[!grepl("recipient_", names(x))]
  return(x)
}

PruneRecipientsFromPhyGraph <- function(phy.graph) {
  new.phy.graph <- phy.graph
  recipient.ids <- which(grepl("recipient_", phy.graph$tip.label))
  if(length(recipient.ids)>0) {
    for(r.index in seq_along(recipient.ids)) {
      recipient.id <- recipient.ids[r.index]
      recipient.ancestor.row <- which(phy.graph$edge[,2]==recipient.id)
      recipient.ancestor.id <- phy.graph$edge[recipient.ancestor.row, 1]
      new.phy.graph$reticulation[which(phy.graph$reticulation[,2]==recipient.id),2] <- recipient.ancestor.id
      new.phy.graph$edge[recipient.ancestor.row,] <- c(NA, NA) # to delete later
      new.phy.graph$edge.length[recipient.ancestor.row] <- NA
    }
    recipient.ids <- sort(recipient.ids, decreasing=TRUE)
    for (r.index in seq_along(recipient.ids)) {
      new.phy.graph$edge[which(new.phy.graph$edge>recipient.ids[r.index])] <- new.phy.graph$edge[which(new.phy.graph$edge>recipient.ids[r.index])] - 1
      new.phy.graph$reticulation[which(new.phy.graph$reticulation>recipient.ids[r.index])] <- new.phy.graph$reticulation[which(new.phy.graph$reticulation>recipient.ids[r.index])] - 1
    }
    new.phy.graph$edge <- new.phy.graph$edge[which(!is.na(new.phy.graph$edge[,1])),]
    new.phy.graph$edge.length <- new.phy.graph$edge.length[which(!is.na(new.phy.graph$edge.length))]
    new.phy.graph$tip.label <- new.phy.graph$tip.label[-recipient.ids]
  }
  return(new.phy.graph)
}
#
# RemoveZeroTerminalsPhygraph <- function(phy.graph) {
#     new.phy.graph <- phy.graph
#     terminal.nodes <- sequence(ape::Ntip(phy.graph)) #this is so as we delete taxa we don't lose track of which ones we've examined
#     for (terminal.index in sequence(ape::Ntip(phy.graph))) {
#         if(new.phy.graph$edge.length[which(new.phy.graph$edge[,2]==terminal.nodes[terminal.index])]==0) {
#             ancestral.node <- new.phy.graph$edge[which(new.phy.graph$edge[,2]==terminal.nodes[terminal.index]),1]
#             new.phy.graph$reticulation[which(new.phy.graph$reticulation==terminal.nodes[terminal.index])] <- ancestral.node
#             new.phy.graph$tip.label <- new.phy.graph$tip.label[-terminal.nodes[terminal.index]] #delete from tips
#             new.phy.graph$edge.length <- new.phy.graph$edge.length[-which(new.phy.graph$edge[,2]==terminal.nodes[terminal.index])] # and from brlen
#             new.phy.graph$edge <- new.phy.graph$edge[-which(new.phy.graph$edge[,2]==terminal.nodes[terminal.index]),] #delete from edges
#
#             #renumber
#             new.phy.graph$edge[which(new.phy.graph$edge>terminal.nodes[terminal.index])] <- new.phy.graph$edge[which(new.phy.graph$edge>terminal.nodes[terminal.index])] - 1
#             new.phy.graph$reticulation[which(new.phy.graph$reticulation>terminal.nodes[terminal.index])] <- new.phy.graph$reticulation[which(new.phy.graph$reticulation>terminal.nodes[terminal.index])] - 1
#
#             terminal.nodes[which(terminal.nodes>terminal.nodes[terminal.index])] <- terminal.nodes[which(terminal.nodes>terminal.nodes[terminal.index])] - 1
#         }
#     }
#     return(new.phy.graph)
# }
#
# ReorderPhygraph <- function(phy.graph, order="cladewise") {
#     new.order <- ape::reorder.phylo(ape::as.phylo(phy.graph), order="pruning",index.only=TRUE) #doing this b/c if just do default cladewise ape doesn't actually reorder
#     phy.graph$edge <- phy.graph$edge[new.order,]
#     phy.graph$edge.length <- phy.graph$edge[new.order]
#
#     new.order <- ape::reorder.phylo(ape::as.phylo(phy.graph), order="post",index.only=TRUE)
#     phy.graph$edge <- phy.graph$edge[new.order,]
#     phy.graph$edge.length <- phy.graph$edge[new.order]
#
#     new.order <- ape::reorder.phylo(ape::as.phylo(phy.graph), order=order,index.only=TRUE)
#     phy.graph$edge <- phy.graph$edge[new.order,]
#     phy.graph$edge.length <- phy.graph$edge[new.order]
#
#     return(phy.graph)
# }

#' Convert phylo object to evonet
#'
#' ape::evonet converts a phylo object to evonet, but requires having at least one hybridization event. This lets you convert to evonet without having a hybridization event.
#'
#' @param phy A phylo object (ape's basic tree format)
#' @return An evonet object, suitable for passing as phy.graph into many of BMhyb's functions
#' @export
CreateHybridlessEvonet <- function(phy) {
  phy.graph <- ape::evonet(phy, from=0, to=0) #0 is just a placeholder here; NA not allowed
  phy.graph$reticulation <- phy.graph$reticulation[-1,]
  return(phy.graph)
}

ConvertPhyAndFlowToPhygraph <- function(phy, flow) {
  flow.aggregate <- LumpIntoClades(phy, flow)
  phy.graph <- CreateHybridlessEvonet(phy)
  for (i in sequence(nrow(flow.aggregate))) {
    phy.graph <- AddHybridization(phy.graph, from.clade=strsplit(flow.aggregate$donor.clades[i], ",")[[1]], to.clade=strsplit(flow.aggregate$recipient.clades[i], ",")[[1]],  time.from.root=flow.aggregate$time.from.root.donor[i], ghost.length=flow.aggregate$time.from.root.recipient[i] - flow.aggregate$time.from.root.donor[i])
  }
  attr(phy.graph, "order")<- NULL
  phy.graph <- RenumberPhygraph(phy.graph)
  return(phy.graph)
}

#' Add hybrid events to a phy.graph
#'
#' Given an evonet object, and info on where the gene flow is from and to, and when this occurs, add a hybridization event. The edges things move from and to are specified by the list of descendant taxa of those edges (basically the edge is the subtending branch for the clade). You do not have to list all taxa, only those spanning the node at the end of the edge. You can enter a single taxon to have gene flow to or from a terminal branch. You also ideally will specify when the gene flow happens. This can be given as time from the root of the tree to when the event starts or time from the tip of the tree back to when the gene flow starts (but you must give one of these). If gene flow goes through an unsampled ghost intermediate, you can enter the length of time it spends there. If you do not specify any of these, flow is assumed to directly from the source to the recipient, with the time set at the start of the recent of the two branches (i.e., if flow goes from taxon A to taxon D, if D is younger the flow is assumed to happen partway up the terminal branch of A to directly connect to the start of D.
#'
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param from.clade A vector of names specifying taxa spanning the node descended from the focal edge for the start of the hybridization event
#' @param to.clade A vector of names specifying taxa spanning the node descended from the focal edge for the start of the hybridization event
#' @param time.from.root When the hybridization event starts, as measured from the root of the tree
#' @param time.from.tip When the hybridization event starts, as measured from the tips of the tree (assumed to be coeval)
#' @param ghost.length How long the hybrid genes spend in an unsampled species before arriving in their recipient
#' @return An evonet object with the new hybridization event
#' @export
AddHybridization <- function(phy.graph, from.clade, to.clade, time.from.root=NULL, time.from.tip=NULL, ghost.length=0) {
  heightnode <- function(node, phy) {
      return(phytools::nodeheight(phy, node))
  }
  heights <- sapply(sequence(ape::Nnode(phy.graph) + ape::Ntip(phy.graph)), heightnode, phy=phy.graph) # does not line up with edge matrix or edge lengths, remember
  donor.node <- ape::getMRCA(phy.graph, from.clade) #returns NULL if only one taxon given
  if(length(from.clade)==1) {
    donor.node <- which(phy.graph$tip.label==from.clade)
  }
  recipient.node <- ape::getMRCA(phy.graph, to.clade)
  if(length(to.clade)==1) {
    recipient.node <- which(phy.graph$tip.label==to.clade)
  }
  donor.height.from.root <- NA

  donor.node.anc <- phy.graph$edge[which(phy.graph$edge[,2]==donor.node),1]
  recipient.node.anc <- phy.graph$edge[which(phy.graph$edge[,2]==recipient.node),1]


  if(!is.null(time.from.root)) {
    donor.height.from.root <- time.from.root
  } else if (!is.null(time.from.tip)) {
    donor.height.from.root <- max(ape::vcv(phy.graph)) - time.from.tip
  } else {
    donor.height.from.root <- max(c(heights[donor.node.anc], heights[recipient.node.anc]))
  }
  recipient.height.from.root <- donor.height.from.root + ghost.length
  new.donor <- ape::Ntip(phy.graph) + 1
  new.recipient <- new.donor + 1
  phy.graph <- AddNodeToPhygraph(below.node=donor.node, depth.below=heights[donor.node] - donor.height.from.root,  phy.graph=phy.graph, tip.label=paste0("donor_", nrow(phy.graph$reticulation)+1), terminal.length=recipient.height.from.root - donor.height.from.root)
  # if matching to internal node, have to renumber since its number was updated by earlier step
  phy.graph <- AddNodeToPhygraph(below.node=ifelse(recipient.node >= (ape::Ntip(phy.graph)+1), recipient.node+1, recipient.node ), depth.below=heights[recipient.node] - recipient.height.from.root,  phy.graph=phy.graph, tip.label=paste0("recipient_", nrow(phy.graph$reticulation)+1), terminal.length=0)
  phy.graph$reticulation <- rbind(phy.graph$reticulation, c(new.donor, new.recipient))
  attr(phy.graph, "order")<- NULL
  phy.graph <- RenumberPhygraph(phy.graph)
  phy.graph$edge <- matrix(as.integer(phy.graph$edge), ncol=2)
  phy.graph$Nnode <- as.integer(max(phy.graph$edge[,1]-ape::Ntip(phy.graph)))
  return(phy.graph)
}


#' Simulate trait data
#'
#' For a given phylogenetic network generate tip data. Any values not specified use default values
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param sigma.sq The Brownian motion wiggle rate
#' @param mu The population mean (in the absence of hybridization)
#' @param bt The beta value (multiplier on expected value for each hybridization event)
#' @param vh The burst of variance that comes from a hybridization event
#' @param SE Uniform uncertainty at the tips
#' @param measurement.error Uncertainty at the tips, especially if it varies between species
#' @param gamma In a hybridization event, what proportion of the trait comes from the donating parent. 0.5 means half comes from each parent
#' @param exclude.donors.recipients If TRUE, do not generate for any donors or recipient placeholder taxa
#' @return A vector of trait values
#' @export
#' @examples
#' network <- SimulateNetwork(ntax=5, nhybridizations=2)
#' tips <- SimulateTips(network, mu=1.1, bt=3, vh=1.1, SE=1)
SimulateTips <- function(phy.graph, sigma.sq=1, mu=0, bt=1, vh=0, SE=0, measurement.error=0, gamma=0.5, exclude.donors.recipients=TRUE) {
  phy.graph <- PruneRecipientsFromPhyGraph(phy.graph)
  means.modified <- ComputeMeans(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma)
  V.modified <- ComputeVCV(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma)
  tips <- MASS::mvrnorm(1, mu=means.modified, Sigma=V.modified, tol=1e-100)
  if(exclude.donors.recipients) {
    tips <- tips[!grepl("donor_", names(tips))]
    tips <- tips[!grepl("recipient_", names(tips))]
  }
  return(tips)
}

#' Simulate a phylogenetic network
#'
#' This uses a birth death process (TreeSim::sim.bd.taxa.age) to make a tree, then randomly adds hybridization events. The events are placed uniformly with time (not with numbers of taxa). If you use the phy.graph argument, you can pass in an existing phylogenetic network and it will add hybridization events to that; if you use a phy argument, it will add hybridization events to that. Note that currently there is no checking for multiple events between the same two branches. While hybridization events happen between taxa alive at the same instant of time, it is possible that the donor taxon later goes extinct with no descendants (other than the taxa of hybrid origin). These are basically ghost lineages, and this process (which then looks like gene flow going forward in time) is permitted if allow.ghost is TRUE.
#'
#' @param ntax How many surviving taxa to have on the tree (extinct taxa are pruned, with the exception of donors if allow.ghost=TRUE)
#' @param nhybridizations How many hybridization events to have
#' @param birth Birth rate (instantaneous rate)
#' @param death Extinction rate (instantaneous)
#' @param sample.f What fraction of taxa alive at the present to sample (resulting in the final ntax)
#' @param tree.height Root to tip height of the final tree
#' @param allow.ghost Allow a hybridization events from an unsampled ancestor
#' @param phy.graph If not NULL, uses this network and adds hybridizations to it to reach nhybridizations in total
#' @param phy If not NULL, uses this tree to create hybridization events on
#' @return A phy.graph object with hybridizations
#' @export
#'
#' @examples
#' p <- SimulateNetwork(ntax=10 ,nhybridizations=2)
#' plot(p)
SimulateNetwork <- function(ntax=100, nhybridizations=10, birth = 1, death = 1, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE, phy.graph=NULL, phy=NULL) {
    done = FALSE
    if(is.null(phy.graph) & is.null(phy)) {
        phy <-  TreeSim::sim.bd.taxa.age(n=ntax, numbsim=1, lambda=birth, mu=death, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
        phy.graph <- ape::evonet(phy, from=0, to=0) #0 is just a placeholder here; NA not allowed
        phy.graph$reticulation <- phy.graph$reticulation[-1,]
    } else if (!is.null(phy)) {
        phy.graph <- ape::evonet(phy, from=0, to=0) #0 is just a placeholder here; NA not allowed
        phy.graph$reticulation <- phy.graph$reticulation[-1,]
    }



    #now use phytools::bind.tip() to add a taxon to the tree for receiving and donating gene flow (contemporaneously if all.ghost=FALSE). Will then delete this tip but keep the internal node (unless we want a ghost lineage).

    while(!done) {
        donor.height.from.root <- stats::runif(1, 0, max(ape::vcv(phy.graph)))
        recipient.height.from.root <- donor.height.from.root
        if(allow.ghost) {
            recipient.height.from.root <- stats::runif(1, donor.height.from.root, max(ape::vcv(phy.graph)))
        }
        heightnode <- function(node, phy) {
            return(phytools::nodeheight(phy, node))
        }
        GetEdgeLengthByNodeNumber <- function(node, phy) {
            final.length <- phy$edge.length[which(phy$edge[,2]==node)]
            if(length(final.length)==0) { # we're at the root, with no root.edge
                final.length <- 0
            }
            return(final.length)
        }
        heights <- sapply(sequence(ape::Nnode(phy.graph) + ape::Ntip(phy.graph)), heightnode, phy=phy.graph) # does not line up with edge matrix or edge lengths, remember
        lengths <- sapply(sequence(ape::Nnode(phy.graph) + ape::Ntip(phy.graph)), GetEdgeLengthByNodeNumber, phy=phy.graph)
        subtending.node.heights <- heights - lengths
        eligible.starts.donors <- which(subtending.node.heights<donor.height.from.root)
        eligible.ends.donors <- which(heights>donor.height.from.root)
        eligible.donors <- intersect(eligible.starts.donors, eligible.ends.donors)
        eligible.starts.recipients <- which(subtending.node.heights<recipient.height.from.root)
        eligible.ends.recipients <- which(heights>recipient.height.from.root)
        eligible.recipients <- intersect(eligible.starts.recipients, eligible.ends.recipients)
        recipient.node <- 0
        donor.node <- 0
        while(recipient.node == donor.node) {
            recipient.node <- sample(eligible.recipients, 1)
            donor.node <- sample(eligible.donors, 1)
            # if((recipient.node %in% phy.graph$reticulation) | (donor.node %in% phy.graph$reticulation) | (phytools::getParent(phy.graph, donor.node) %in% phy.graph$reticulation) | (phytools::getParent(phy.graph, recipient.node) %in% phy.graph$reticulation)) {
            #   recipient.node <- 0
            #   donor.node <- 0
            # }
        }

        new.donor <- ape::Ntip(phy.graph) + 1
        new.recipient <- new.donor + 1
        phy.graph <- AddNodeToPhygraph(below.node=donor.node, depth.below=heights[donor.node] - donor.height.from.root,  phy.graph=phy.graph, tip.label=paste0("donor_", nrow(phy.graph$reticulation)+1), terminal.length=recipient.height.from.root - donor.height.from.root)
        # if matching to internal node, have to renumber since its number was updated by earlier step
        phy.graph <- AddNodeToPhygraph(below.node=ifelse(recipient.node >= (ape::Ntip(phy.graph)+1), recipient.node+1, recipient.node ), depth.below=heights[recipient.node] - recipient.height.from.root,  phy.graph=phy.graph, tip.label=paste0("recipient_", nrow(phy.graph$reticulation)+1), terminal.length=0)
        phy.graph$reticulation <- rbind(phy.graph$reticulation, c(new.donor, new.recipient))
        if(nrow(phy.graph$reticulation)==nhybridizations) {
            done=TRUE
        }
    }
    attr(phy.graph, "order")<- NULL
    phy.graph <- PruneRecipientsFromPhyGraph(RenumberPhygraph(phy.graph))
    return(phy.graph)
}

#allow.ghost allows ghost lineage: something that persists for awhile, hybridizes, goes extinct. Otherwise, hybridization events must between coeval edges with extant descendants
# SimulateNetwork <- function(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 1, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE) {
#     done = FALSE
#     used.recipients <- c()
#     available.recipient.ids <- sequence(ntax.nonhybrid + ntax.hybrid)
#     flow <- data.frame()
#     phy<-NA
#     phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=death, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
#     if(origin.type=="clade" && ntax.hybrid==1) {
#         warning("For ntax.hybrid = 1 and clade sampling, this will do individual sampling instead (which is equivalent in this case)")
#         origin.type<-"individual"
#     }
#     if(origin.type=="clade") {
#         while(is.na(GetClade(phy, ntax.hybrid))) { #not all trees of a given size have a clade of a given size, so may need to resimulate it
#             phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=0.5, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
#         }
#     }
#     while(!done) {
#         donors <- c()
#         recipients <- c()
#         recipient.ids <- c()
#         focal.node <- c()
#         if (origin.type=="clade") {
#             focal.node <- GetClade(phy, ntax.hybrid)
#             if(is.na(focal.node)) {
#                 done=FALSE
#                 break()
#             }
#             recipients <- phy$tip.label[getDescendants(phy, node=focal.node)]
#             recipients <- recipients[!is.na(recipients)] #since we want just the tips
#             recipient.ids <- which(phy$tip.label %in% recipients)
#             used.recipients <- append(used.recipients, recipients)
#         } else {
#             focal.node<-sample(available.recipient.ids, 1, replace=FALSE)
#             recipient.ids <- focal.node
#             recipients <- phy$tip.label[focal.node]
#             used.recipients <- append(used.recipients, recipients)
#         }
#         available.recipient.ids <- available.recipient.ids[!available.recipient.ids %in% recipient.ids]
#         longest.from.root <- nodeheight(phy, node=focal.node)
#         shortest.from.root <- nodeheight(phy, node=GetAncestor(phy, focal.node))
#         all.heights <- nodeHeights(phy)
#         #idea here: take a recipient clade. The gene flow must happen on its stem edge, which starts at shortest.from.root and goes up to longest.from.root. Gene flow can't go back in time
#         qualifying.lower <- which(all.heights[,1]<longest.from.root) #if this is false, gene flow goes back in time
#         qualifying.upper <- sequence(dim(all.heights)[1]) #in general, gene flow can go forward in time via ghost lineages
#         if(!allow.ghost) {
#             qualifying.upper <- which(all.heights[,2]>shortest.from.root) #if no ghost lineages, then there must be temporal overlap between the donor and recipient lineages. So the tipward end of the donor edge must be later than the rootward end of the recipient edge
#         }
#         qualifying.upper <- qualifying.upper[which(phy$edge[qualifying.upper,2]!=focal.node)] #let's not hybridize with ourselves
#         qualifying.all <- qualifying.upper[qualifying.upper %in% qualifying.lower]
#         if(length(qualifying.all)==0) {
#             break()
#         }
#         donor.edge <- sample(qualifying.all, 1)
#         donors <- phy$tip.label[getDescendants(phy, phy$edge[donor.edge,2])]
#         donors <- donors[!is.na(donors)] #getDescendants includes all descendant nodes, including internal ones. We just want the terminal taxa
#         time.in <- runif(1, min=max(all.heights[donor.edge,1],shortest.from.root), max=longest.from.root)
#         time.out <- runif(1, min=all.heights[donor.edge,1], max=min(time.in, all.heights[donor.edge,2]))
#         if (!allow.ghost) {
#             time.in <- runif(1, min=max(shortest.from.root, all.heights[donor.edge,1]), max=min(longest.from.root, all.heights[donor.edge,2])) #if no ghost lineages, must move from the overlapping interval
#             time.out <- time.in
#         }
#         pairs <- expand.grid(donors, recipients)
#         for (pairs.index in sequence(dim(pairs)[1])) {
#             flow <- rbind(flow, data.frame(donor=pairs[pairs.index,1], recipient=pairs[pairs.index,2], gamma=flow.proportion, time.from.root.donor=time.out, time.from.root.recipient=time.in, stringsAsFactors=FALSE))
#         }
#         if(length(used.recipients)==ntax.hybrid) {
#             done=TRUE
#         }
#         if(length(CheckFlow(phy, flow)$problem.taxa)>0) {
#             done=FALSE
#         }
#     }
#     flow$donor <- as.character(flow$donor)
#     flow$recipient <- as.character(flow$recipient)
#     flow$gamma <- as.numeric(as.character(flow$gamma))
#     flow$time.from.root.donor <-as.numeric(as.character(flow$time.from.root.donor))
#     flow$time.from.root.recipient <-as.numeric(as.character(flow$time.from.root.recipient))
#     return(list(phy=phy, flow=flow))
# }

# PlotNetwork <- function(phy, flow, col.non="black", col.hybrid="red", col.donor="blue", name.padding=1.5, cex=1, xlab="", bty="n", head.length=0.2, edge.width=2, col.tree="darkgray", col.arrow="red", arrow.width=1, try.rotations=FALSE, ...) {
#     phy<-reorder(phy, "pruningwise")
#     if(try.rotations) {
#         phy <- phytools::rotateNodes(phy, nodes="all")
#     }
#     phy4 <- as(phy, "phylo4")
#     xxyy <- phyloXXYY(phy4)
#     #plot(phy4)
#     plot(x=c(min(xxyy$xx), name.padding*max(xxyy$xx)), y=range(xxyy$yy), type="n", xaxt="n", xlab=xlab, yaxt="n", ylab="", bty=bty, ...)
#     arrows(x0=xxyy$segs$v0x, x1=xxyy$segs$v1x, y0=xxyy$segs$v0y, y1=xxyy$segs$v1y, length=0, lwd=edge.width, col=col.tree)
#     arrows(x0=xxyy$segs$h0x, x1=xxyy$segs$h1x, y0=xxyy$segs$h0y, y1=xxyy$segs$h1y, length=0, lwd=edge.width, col=col.tree)
#     label.colors <- rep(col.non, Ntip(phy))
#     for (i in sequence(Ntip(phy))) {
#         if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$donor) {
#             label.colors[i]<-col.donor
#         }
#         if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$recipient) {
#             label.colors[i]<-col.hybrid
#         }
#
#     }
#     text(x=rep(max(xxyy$xx), Ntip(phy)), y=xxyy$yy[which(edges(phy4)[xxyy$eorder,2] %in% sequence(Ntip(phy)))], names(getNode(phy4, xxyy$torder)), col=label.colors, pos=4, cex=cex)
#     for (i in sequence(dim(flow)[1])) {
#         recipient.node <- getNode(phy4, flow$recipient[i])
#         recipient.path <- c(recipient.node, ancestors(phy4, recipient.node))
#         recipient.path.heights <- nodeHeight(phy4, recipient.path, from="root")
#         valid.recipients <- recipient.path[which(recipient.path.heights > flow$time.from.root.recipient[i])]
#         recipient <- valid.recipients[length(valid.recipients)] #do it from the earliest qualifying tipward node
#         if(length(valid.recipients)>1 && length(which(recipient.path.heights==flow$time.from.root.recipient[i]))>0) { #the latter condition means we're moving to an existing node
#             recipient <- valid.recipients[length(valid.recipients)-1]
#         }
#         y1 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == recipient)]
#         donor.node <- getNode(phy4, flow$donor[i])
#         donor.path <- c(donor.node, ancestors(phy4, donor.node))
#         donor.path.heights <- nodeHeight(phy4, donor.path, from="root")
#         valid.donors <- donor.path[which(donor.path.heights > flow$time.from.root.donor[i])]
#         donor <- valid.donors[length(valid.donors)] #do it from the earliest qualifying tipward node
#         y0 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == donor)]
#         arrows(x0=flow$time.from.root.donor[i]/max(vcv(phy)), x1=flow$time.from.root.recipient[i]/max(vcv(phy)), y1=y1, y0=y0, col=col.arrow, lwd=arrow.width, length=head.length) #rescale since it goes from zero to 1 in height
#         #grid.arrows(x=c(flow$time.from.root[i],flow$time.from.root[i]), y=c(y0, y1))
#     }
# }
#
# PlotNetwork <- function(phy.graph, col.non="black", col.hybrid="red", col.donor="blue", name.padding=1.5, cex=1, xlab="", bty="n", head.length=0.2, edge.width=2, col.tree="darkgray", col.arrow="red", arrow.width=1, try.rotations=FALSE, ...) {
#     #phy<-reorder(phy, "pruningwise")
#     #if(try.rotations) {
#     #    phy <- phytools::rotateNodes(phy, nodes="all")
#     #}
#     phy.graph$edge.length[phy.graph$edge.length<0] <- 0
#     phy4 <- as(ape::as.phylo(phy.graph), "phylo4")
#     xxyy <- phylobase::phyloXXYY(phy4)
#     #plot(phy4)
#     plot(x=c(min(xxyy$xx), name.padding*max(xxyy$xx)), y=range(xxyy$yy), type="n", xaxt="n", xlab=xlab, yaxt="n", ylab="", bty=bty, ...)
#     arrows(x0=xxyy$segs$v0x, x1=xxyy$segs$v1x, y0=xxyy$segs$v0y, y1=xxyy$segs$v1y, length=0, lwd=edge.width, col=col.tree)
#     arrows(x0=xxyy$segs$h0x, x1=xxyy$segs$h1x, y0=xxyy$segs$h0y, y1=xxyy$segs$h1y, length=0, lwd=edge.width, col=col.tree)
#     label.colors <- rep(col.non, Ntip(phy))
#     for (i in sequence(Ntip(phy))) {
#         if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$donor) {
#             label.colors[i]<-col.donor
#         }
#         if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$recipient) {
#             label.colors[i]<-col.hybrid
#         }
#
#     }
#     text(x=rep(max(xxyy$xx), Ntip(phy)), y=xxyy$yy[which(edges(phy4)[xxyy$eorder,2] %in% sequence(Ntip(phy)))], names(getNode(phy4, xxyy$torder)), col=label.colors, pos=4, cex=cex)
#     # for (i in sequence(dim(flow)[1])) {
#     #     recipient.node <- getNode(phy4, flow$recipient[i])
#     #     recipient.path <- c(recipient.node, ancestors(phy4, recipient.node))
#     #     recipient.path.heights <- nodeHeight(phy4, recipient.path, from="root")
#     #     valid.recipients <- recipient.path[which(recipient.path.heights > flow$time.from.root.recipient[i])]
#     #     recipient <- valid.recipients[length(valid.recipients)] #do it from the earliest qualifying tipward node
#     #     if(length(valid.recipients)>1 && length(which(recipient.path.heights==flow$time.from.root.recipient[i]))>0) { #the latter condition means we're moving to an existing node
#     #         recipient <- valid.recipients[length(valid.recipients)-1]
#     #     }
#     #     y1 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == recipient)]
#     #     donor.node <- getNode(phy4, flow$donor[i])
#     #     donor.path <- c(donor.node, ancestors(phy4, donor.node))
#     #     donor.path.heights <- nodeHeight(phy4, donor.path, from="root")
#     #     valid.donors <- donor.path[which(donor.path.heights > flow$time.from.root.donor[i])]
#     #     donor <- valid.donors[length(valid.donors)] #do it from the earliest qualifying tipward node
#     #     y0 <- xxyy$yy[which(edges(phy4)[xxyy$eorder,2] == donor)]
#     #     arrows(x0=flow$time.from.root.donor[i]/max(vcv(phy)), x1=flow$time.from.root.recipient[i]/max(vcv(phy)), y1=y1, y0=y0, col=col.arrow, lwd=arrow.width, length=head.length) #rescale since it goes from zero to 1 in height
#     #     #grid.arrows(x=c(flow$time.from.root[i],flow$time.from.root[i]), y=c(y0, y1))
#     # }
# }



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

# UnlumpIntoTaxa <- function(lumped.flow) {
#     flow <- data.frame()
#     for (event.index in sequence(nrow(lumped.flow))) {
#         donor.taxa <- strsplit(lumped.flow$donor.clades[event.index], ",")[[1]]
#         recipient.taxa <- strsplit(lumped.flow$recipient.clades[event.index], ",")[[1]]
#         flow <- rbind(flow, expand.grid(donor=donor.taxa, recipient=recipient.taxa, gamma=lumped.flow$gamma[event.index], time.from.root.donor = lumped.flow$time.from.root.donor[event.index], time.from.root.recipient =  lumped.flow$time.from.root.recipient[event.index], event=event.index))
#     }
#     return(flow)
# }

# AttachHybridsToDonor <- function(phy, flow, suffix="_DUPLICATE") {
#     flow.clades <- LumpIntoClades(phy, flow)
#     phy.merged <- phy
#     for (i in sequence(dim(flow.clades)[1])) {
#         taxa.to.retain <- strsplit(flow.clades$recipient.clades[i], ",")[[1]]
#         donor.taxa <- strsplit(flow.clades$donor.clades[i], ",")[[1]]
#         pulled.clade <- c()
#         if(length(taxa.to.retain)>1) {
#             pulled.clade <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% taxa.to.retain])
#         } else {
#             pulled.clade <- structure(list(edge = structure(c(3L, 3L, 1L, 2L), .Dim = c(2L,
#                                                                                         2L)), edge.length = rep(flow.clades$time.from.root.recipient[i],2),
#                                            tip.label = c(taxa.to.retain, "DUMMY"), Nnode = 1L), .Names = c("edge",
#                                                                                                            "edge.length", "tip.label", "Nnode"), class = "phylo")
#         }
#         pulled.clade$tip.label <- paste(pulled.clade$tip.label, suffix, sep="")
#         attachment.crown.node <- which(phy.merged$tip.label==donor.taxa[1])
#         if(length(donor.taxa)>1) {
#             attachment.crown.node <- findMRCA(phy.merged, tips=donor.taxa, type=c("node"))
#         }
#         attachment.stem.node <- GetAncestor(phy.merged, attachment.crown.node)
#
#         pulled.clade$root.edge<-max(vcv(phy)) - max(vcv(pulled.clade)) - flow.clades$time.from.root.donor[i]
#
#         phy.merged <- bind.tree(phy.merged, pulled.clade, attachment.crown.node, position=phy.merged$edge.length[which(phy.merged$edge[,2]==attachment.crown.node)] - (flow.clades$time.from.root.donor[i]-nodeheight(phy.merged, attachment.stem.node)))
#         if(length(taxa.to.retain)==1) {
#             phy.merged <- drop.tip(phy.merged, paste("DUMMY", suffix, sep=""))
#         }
#     }
#     return(phy.merged)
# }
#

# SimulateTipData <- function(phy, flow, params, measurement.error, suffix="_DUPLICATE") {
#     #	flow.clades <- LumpIntoClades(phy, flow)
#     #	recipients <- c()
#     #	for (i in sequence(dim(flow.clades)[1])) {
#     #		recipients<-c(recipients, strsplit(flow.clades$recipient.clades, ",")[[1]])
#     #	}
#     #	if(length(unique(recipients)) != length(recipients)) {
#     #		print(flow)
#     #		print(flow.clades)
#     #		stop("This function only works if each taxon of hybrid origin only appears once as a recipient")
#     #	}
#     phy.merged <- AttachHybridsToDonor(phy, flow, suffix=suffix)
#     phy.merged$edge.length <- phy.merged$edge.length*params["sigma.sq"]
#     phy.merged$edge.length[which(phy.merged$edge[,2] <= Ntip(phy.merged))] <- phy.merged$edge.length[which(phy.merged$edge[,2] <= Ntip(phy.merged))]+params["SE"]
#     tips <- sim.char(phy.merged, par=matrix(1, nrow=1, ncol=1), nsim=1, model="BM", root=params["mu"])[,,1]	#we have already done the scaling and SE
#     hybrid.name.root <- gsub(suffix, "", phy.merged$tip.label[grepl(suffix,  phy.merged$tip.label)])
#     for (i in sequence(length(hybrid.name.root))) {
#         focal.tips <- tips[c(hybrid.name.root[i], paste(hybrid.name.root[i], suffix, sep=""))]
#         focal.gamma <- flow$gamma[which(flow$recipient==hybrid.name.root[i])]
#         focal.tips.bt <- focal.tips + log(params['bt'])
#         tip.mean <- focal.gamma*focal.tips.bt[2] + (1-focal.gamma) * focal.tips.bt[1]
#         tip.final <- rnorm(1, mean=tip.mean, sd=sqrt(params['vh']))
#         tips[hybrid.name.root[i]]<-tip.final
#     }
#     tips<-tips[!grepl(suffix, names(tips))]
#     if(!is.null(measurement.error)) {
#         tips <- tips+measurement.error
#     }
#     return(tips)
# }

#The following short function comes from Ken Takagi at  https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
PlotConvexHull<-function(xcoord, ycoord, lcolor){
    hpts <- grDevices::chull(x = xcoord, y = ycoord)
    hpts <- c(hpts, hpts[1])
    graphics::lines(xcoord[hpts], ycoord[hpts], col = lcolor)
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
            x.grid <- stats::quantile(sims.sub[,1], seq(from=0, to=1, length.out = floor(length(sims.sub[,1])/50)))
            y.grid <- stats::quantile(sims.sub[,2], seq(from=0, to=1, length.out = floor(length(sims.sub[,2])/50)))
            colnames(sims.sub)[1:2] <- c("x", "y")
            # for (x.index in sequence(length(x.grid)-1)) {
            # for (y.index in sequence(length(y.grid)-1)) {
            # relevant.points <- subset(sims.sub, x>=x.grid[x.index] & x<x.grid[x.index+1] & y>=y.grid[y.index] & y<y.grid[y.index+1])
            # points.to.fit <- rbind(points.to.fit, relevant.points[which.max(relevant.points[,3]),])
            # }
            # }
            graphics::plot(x=x.range, y=y.range, xlab=params.of.interest[param.1], ylab=params.of.interest[param.2], type="n", bty="n")

            relevant.point.id <- grDevices::chull(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2])
            #polygon(sims.sub[relevant.point.id,1], sims.sub[relevant.point.id, 2], col="gray", border=NA, fillOddEven=TRUE)
            #polygon(sims.sub[relevant.point.id,1], sims.sub[relevant.point.id, 2], col="gray", border=NA, fillOddEven=FALSE)
            PlotConvexHull(sims.sub[which(sims.sub[,3]<5),1], sims.sub[which(sims.sub[,3]<5),2], "darkgray")
            PlotConvexHull(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2], "black")

            #points(sims.sub[which(sims.sub[,3]<2),1], sims.sub[which(sims.sub[,3]<2),2], col="green", pch="X")
            #contour(interp(points.to.fit[,1], points.to.fit[,2], points.to.fit[,3]), xlab=params.of.interest[param.1], ylab=params.of.interest[param.2], levels=c(1, 2, 5, 10))
            graphics::points(sims.sub[which.min(sims.sub$neglnL),1], sims.sub[which.min(sims.sub$neglnL),2], col="red", pch=20, cex=2)
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
    distance <- as.numeric(stats::dist(rbind(as.vector(new.mat), as.vector(original))))
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
GenerateRandomPositiveDefiniteMatrix <- function(n, ev = stats::runif(n, 0, 10)) {
    Z <- matrix(ncol=n, stats::rnorm(n^2))
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
        starting.val.matrix[i,] <- stats::rlnorm(length(starting.val.center), meanlog=starting.means, sdlog=sd.vector[i])
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

# VerifyActualParams <- function(x) {
#     if(class(x)!="logical" | is.null(names(x))) {
#         stop("This requires a free.params or actual.params argument that is a boolean vector with names")
#     }
# }

GenerateParamLabels <- function() {
    return(c("sigma.sq", "mu", "bt", "vh", "SE"))
}


# New attempt to do network

# get enewick using ape, make an evonet object. Call it phy.net. For igraph object, call phy.i

#' Convert an evonet object into igraph
#'
#' ape can already convert from evonet to igraph; the advantage of this function is that it uses the node ids from the evonet object for labels in igraph.
#'
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @return An igraph network
#'
#' @examples
#' phy <- ape::rcoal(5)
#' phy.evo <- ape::evonet(phy, from=1, to=2)
#' plot(phy.evo) # this is the ape plot
#' phy.igraph <- ConvertEvonetToIgraphWithNodeNumbers(phy.evo)
#' plot(phy.igraph)
#' @export
ConvertEvonetToIgraphWithNodeNumbers <- function(phy.graph) {
    #phy.i <- ape::as.igraph(phy.graph)
    phy.graph$tip.label <- sequence(ape::Ntip(phy.graph)) #converting to ape node numbers
    phy.graph$node.label <- seq(from=(1+ape::Ntip(phy.graph)), to=ape::Ntip(phy.graph)+ape::Nnode(phy.graph), by=1)
    return(ape::as.igraph.evonet(phy.graph))
}




GetAllPathTopologies <- function(phy.graph) {
    phy.i <- ConvertEvonetToIgraphWithNodeNumbers(phy.graph)
    root <- ape::Ntip(phy.graph)+1
    #all.paths <- list(rep(NA, ape::Ntip(phy.graph)))
    all.paths <- data.frame()
    for (taxon in sequence(ape::Ntip(phy.graph))) {
        simple.paths <-  igraph::all_simple_paths(phy.i, from=root, to=taxon)
        local.paths <- unlist(lapply(lapply(simple.paths, as.numeric), paste, collapse="_"))
        local.df <- data.frame(final.taxon=taxon, paths=local.paths, stringsAsFactors=FALSE)
        all.paths <- rbind(all.paths, local.df)
    }

    return(all.paths) #one element per taxon; each element is itself a list of multiple possible paths from the root to that taxon
}

ComputeAllEdges <- function(phy.graph, gamma=0.5) {
    normal.edges <- data.frame(node.from=phy.graph$edge[,1], node.to=phy.graph$edge[,2], length=phy.graph$edge.length, flow.prob=1, type="normal", stringsAsFactors=FALSE)
    #we actually don't need this: flow.prob is prob of a gene being passed on from that node. For regular speciation nodes, it's 1
  #  for (edge.index in sequence(nrow(normal.edges))) {
        #normal.edges$flow.prob[edge.index] <- 1/length(which(normal.edges$node.from == normal.edges$node.from[edge.index])) #so can deal with nodes with outdegree of 1, 3+, etc.
  #  }
    all.edges <- normal.edges
    if(nrow(phy.graph$reticulation)>0) {
        hybrid.edges <- data.frame(node.from=phy.graph$reticulation[,1], node.to=phy.graph$reticulation[,2], length=0, flow.prob=2*gamma, type="hybridflow", stringsAsFactors=FALSE) # 2 x gamma here so that it balances with other paths
        all.edges <- rbind(normal.edges, hybrid.edges)
    }
    return(all.edges)
}


GetHybridNodes <- function(phy.graph, gamma=0.5) {
    all.edges <- ComputeAllEdges(phy.graph, gamma)
    immediate.hybrids <- subset(all.edges, all.edges$type=="hybridflow")$node.to
    hybrid.descendants <- c()
    for (i in sequence(length(immediate.hybrids))) {
        hybrid.descendants <- append(hybrid.descendants, phytools::getDescendants(phy.graph, node=immediate.hybrids[i]))
    }
    hybrid.counts <- table(c(hybrid.descendants, immediate.hybrids))
    hybrid.descendants <- unique(hybrid.descendants)
    return(list(hybrid.descendants=hybrid.descendants, immediate.hybrids=as.numeric(immediate.hybrids), hybrid.counts=hybrid.counts))
}

ScaleAllEdges <- function(phy.graph, sigma.sq=1, mu=0, bt=1, vh=0, SE=0, measurement.error=0, gamma=0.5) {
    all.edges <- ComputeAllEdges(phy.graph, gamma)
    all.edges$length <- all.edges$length*sigma.sq
    if(any(all.edges$length<0)) {
      warning(paste0("Some of the edges were negative in length. If this is tiny, don't worry about it; if not, do. The negative edge with largest magnitude was ", min(all.edges$length), ". We are converting all of the negative edges to zero length."))
      all.edges$length[all.edges$length<0] <- 0
    }
    immediates <- GetHybridNodes(phy.graph, gamma)$immediate.hybrids
    for (immediate.node.index in sequence(length(immediates))) {
        all.edges[which(all.edges$node.from==immediates[immediate.node.index]),]$length <- all.edges[which(all.edges$node.from==immediates[immediate.node.index]),]$length + vh
    }
    for (terminal.taxon in sequence(ape::Ntip(phy.graph))) { # b/c with ape, terminal taxa are 1:Ntip
        all.edges[which(all.edges$node.to==terminal.taxon),]$length <- all.edges[which(all.edges$node.to==terminal.taxon),]$length + SE
    }
    #TODO: add measurement error
    if(length(measurement.error)==1) {
        all.edges[which(all.edges$node.to<=ape::Ntip(phy.graph)),]$length <- all.edges[which(all.edges$node.to<=ape::Ntip(phy.graph)),]$length + measurement.error
    }
    if(length(measurement.error)>1) {
      if(is.null(names(measurement.error))) {
        stop("measurement.error, if it is more than one value, must have names in the vector corresponding to names of the taxa. We do not assume the order of names is the same, so they must be done explicitly")
      } else {
        for(taxon.index in seq_along(measurement.error)) {
          tip.id <- which(phy.graph$tip.label==names(measurement.error)[taxon.index])
          all.edges[which(all.edges$node.to==tip.id),]$length <- all.edges[which(all.edges$node.to==tip.id),]$length + measurement.error[taxon.index]
        }
      }
    }
    return(all.edges)
}

GetProbabilityOfIndividualPath <- function(path, all.edges) {
    path<-strsplit(path, "_")[[1]]
    total.prob <- 1
    for (node.start in sequence(length(path)-1)) {
        total.prob <- total.prob * subset(all.edges, all.edges$node.from==path[node.start] & all.edges$node.to==path[node.start+1])$flow.prob
    }
    return(total.prob)
}

GetProbabilityOfAllPaths <- function(phy.graph, all.edges) {
    all.paths <- GetAllPathTopologies(phy.graph)
    all.paths$probabilities <- sapply(all.paths$paths,GetProbabilityOfIndividualPath, all.edges=all.edges)
    all.paths$normalized.probabilities <- NA
    for (taxon in sequence(length(unique(all.paths$final.taxon)))) {
        matching.rows <- which(all.paths$final.taxon==taxon)
        all.paths$normalized.probabilities[matching.rows] <-   all.paths$probabilities[matching.rows] / sum(all.paths$probabilities[matching.rows])
    }
    return(all.paths)
}

ComputePathPairs <- function(path) {
    path<-strsplit(path, "_")[[1]]
    GetPair <- function(node.start, path) {
      return(c(node.from=path[node.start], node.to=path[node.start+1], node.fromto = paste(path[c(node.start, node.start+1)], collapse="_")))
    }
    result <- sapply(sequence(length(path)-1), GetPair, path=path)
    #path.pairs <- data.frame()
    # for (node.start in sequence(length(path)-1)) {
    #     path.pairs <- rbind(path.pairs, data.frame(node.from=path[node.start], node.to=path[node.start+1], node.fromto = paste(path[c(node.start, node.start+1)], collapse="_"), stringsAsFactors=FALSE))
    # }
    # return(path.pairs)
    #return(data.frame(node.from=as.numeric(result[1,]), node.to=as.numeric(result[2,]), node.fromto=as.character(result[3,]), stringsAsFactors=FALSE))
    return(t(result))
}

#' Compute the variance-covariance matrix
#'
#' Creates a variance-covariance matrix for a network and parameters.
#'
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param sigma.sq Value for sigma squared
#' @param mu Value for state at the root
#' @param bt Value for beta parameter
#' @param vh Value for Vh, the variance that comes from a hybridization event
#' @param SE Standard error
#' @param measurement.error How much uncertainty there is in tip values; a single number is applied to all taxa, a vector is applied to the corresponding taxa
#' @param gamma In a hybridization event, what proportion of the trait comes from the donating parent. 0.5 means half comes from each parent
#'
#' @return Returns the variance-covariance matrix
#' @export
ComputeVCV <- function(phy.graph, sigma.sq=1, mu=0, bt=1, vh=0, SE=0, measurement.error=0, gamma=0.5) {
    phy.graph <- PruneRecipientsFromPhyGraph(phy.graph)
    all.edges <- ScaleAllEdges(phy.graph=phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma)
    all.paths <- GetProbabilityOfAllPaths(phy.graph=phy.graph, all.edges=all.edges)
    VCV <- matrix(0, nrow=ape::Ntip(phy.graph), ncol=ape::Ntip(phy.graph))
    rownames(VCV) <- phy.graph$tip.label
    colnames(VCV) <- phy.graph$tip.label
    for (row.index in sequence(nrow(VCV))) {
        for (col.index in seq(from=row.index, to=ncol(VCV), by=1)) {
            paths.left <- subset(all.paths, all.paths$final.taxon==row.index)
            paths.top <- subset(all.paths, all.paths$final.taxon==col.index)
            value <- 0
            for(left.path.index in sequence(nrow(paths.left))) {
                for(top.path.index in sequence(nrow(paths.top))) {
                    this.path.left <- ComputePathPairs(paths.left$paths[left.path.index])
                    weight.left <- paths.left$normalized.probabilities[left.path.index]
                    this.path.top <- ComputePathPairs(paths.top$paths[top.path.index])
                    weight.top <- paths.top$normalized.probabilities[top.path.index]
                    #matches <- this.path.left[this.path.left$node.fromto %in% this.path.top$node.fromto,]
                    matching.indices <- which(this.path.left[,'node.fromto'] %in% this.path.top[,'node.fromto'])
                    for (match.index in seq_along(matching.indices)) {
                        value <- value + (weight.left * weight.top) * subset(all.edges, all.edges$node.from==this.path.left[matching.indices[match.index], 'node.from'] & all.edges$node.to==this.path.left[matching.indices[match.index], 'node.to'])$length
                    }
                }
            }
            VCV[row.index, col.index] <- value
            VCV[col.index, row.index] <- value

        }
    }
    return(VCV)
}

ComputeMeans <- function(phy.graph, sigma.sq=1, mu=0, bt=1, vh=0, SE=0, measurement.error=0, gamma=0.5) {
    phy.graph <- PruneRecipientsFromPhyGraph(phy.graph)
    means.vector <- rep(mu, ape::Ntip(phy.graph))
    hybrids <- GetHybridNodes(phy.graph, gamma)$hybrid.counts #hybrid.descendants=hybrid.descendants, immediate.hybrids, hybrid.counts
    #hybrids <- unique(c(hybrids$hybrid.descendants, hybrids$immediate.hybrids))
    hybrids <- hybrids[which(as.numeric(names(hybrids))<=ape::Ntip(phy.graph))]
    for(i in sequence(length(hybrids))) {
        means.vector[as.numeric(names(hybrids)[i])] <- means.vector[as.numeric(names(hybrids)[i])] + hybrids[i]*log(bt)
    }

    names(means.vector) <- phy.graph$tip.label
    return(means.vector)
}

# ComputeLikelihoodUsingChol <- function(parameters, phy.graph, traits, measurement.error=0, gamma=0.5, do.Higham.correction=FALSE, do.Brissette.correction=FALSE, do.DE.correction=FALSE) {
#     badval<-(0.5)*.Machine$double.xmax
#     sigma.sq=1
#     mu=0
#     bt=1
#     vh=0
#     SE=0
#     for(i in seq_along(parameters)){
#         assign(names(parameters)[i],parameters[i])
#     }
#     phy.graph <- PruneRecipientsFromPhyGraph(phy.graph)
#     means.modified <- ComputeMeans(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma)
#     V.modified <- ComputeVCV(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma)
#
#     if(sigma.sq <0 || vh<0 || bt <= 0.0000001 || SE < 0) {
#         return(badval)
#     }
#
#
#     V.chol <- chol(PruneDonorsRecipientsFromVCV(V.modified))
#
#     traits <- traits[match(rownames(V.chol), names(traits))] #reorder so traits in same order as VCV and means
#   #  traits <-
# #    try(NegLogML_raw <- ((ncol(V.modified)/2)*log(2*pi)+(1/2)*t(traits[!is.na(names(traits))]-(means.modified))%*%corpcor::pseudoinverse(V.modified)%*%(traits[!is.na(names(traits))]-(means.modified)) + (1/2)*determinant(V.modified, logarithm=TRUE)$modulus)[1,1])
# #    print(NegLogML_raw)
#
#
#     traits.modified <- solve((V.chol)) %*% (traits-means.modified)
#
#     #PruneDonorsRecipientsFromMeans(means.modified)
#     #means.modified <- solve((V.chol)) %*% PruneDonorsRecipientsFromMeans(means.modified)
#
#    #print(solve(t(V.chol))%*%t(V.chol) )
#
#
#     NegLogML <- NULL
#     #try(NegLogML <- (ncol(V.chol)/2)*log(2*pi)+(1/2)*t(traits.modified-means.modified)%*%diag(ncol(V.chol))%*%(traits.modified-means.modified), silent=TRUE)
#     try(NegLogML <- -(sum(dnorm(traits.modified, log=TRUE))))
#     if(is.null(NegLogML)) {
#         NegLogML <- badval
#     }
#     return(NegLogML[1])
# }

TryComputeLikelihood <- function(...) {
  NegLogL<-(0.5)*.Machine$double.xmax
  try(NegLogL <- ComputeLikelihood(...))
  return(NegLogL)
}

#' Compute the likelihood for a set of parameters
#'
#' Computes likelihood for a given network, set of traits, and parameters.
#'
#' This takes an ape::evonet object. If all you have is a tree (an ape::phylo object), you can use CreateHybridlessEvonet() to convert the tree to an evonet object. You can then use the AddHybridization() function to add hybrid events to this object. Note that networks created in this way can, by chance, result in orders of nodes in the internal edge matrix that cause ape's reorder.phylo function to crash, which is called in many of the plot and write functions. You can still use the plot functions in this package, however.
#'
#' @param parameters Named vector of parameter values; expected names are sigma.sq, mu, SE, bt, and vh
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param traits A vector of trait values, with names equal to the names of taxa on the phylogeny
#' @param measurement.error How much uncertainty there is in tip values; a single number is applied to all taxa, a vector is applied to the corresponding taxa
#' @param gamma In a hybridization event, what proportion of the trait comes from the donating parent. 0.5 means half comes from each parent
#' @param do.Higham.correction Variance-covariance matrices for this model are sometimes poorly conditioned; this is a hack to reduce the impact of that
#' @param do.Brissette.correction Applies method of Brissette et al. 2007 to also try to fix matrix condition
#' @param do.DE.correction Inspired by Mishra, Sudhanshu K. "The nearest correlation matrix problem: Solution by differential evolution method of global optimization." (2007)
#'
#' @return Returns the negative log likelihood
#' @export
ComputeLikelihood <- function(parameters, phy.graph, traits, measurement.error=0, gamma=0.5, do.Higham.correction=FALSE, do.Brissette.correction=FALSE, do.DE.correction=FALSE) {
    badval<-(0.5)*.Machine$double.xmax
    sigma.sq=1
    mu=0
    bt=1
    vh=0
    SE=0
    for(i in seq_along(parameters)){
        assign(names(parameters)[i],parameters[i])
    }
    phy.graph <- PruneRecipientsFromPhyGraph(phy.graph)
    means.modified <- PruneDonorsRecipientsFromMeans(ComputeMeans(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma))
    V.modified <- PruneDonorsRecipientsFromVCV(ComputeVCV(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error, gamma=gamma))

    if(sigma.sq <0 || vh<0 || bt <= 0.0000001 || SE < 0) {
        return(badval)
    }
    likelihood.penalty <- 0
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
            likelihood.penalty <- 10+stats::dist(rbind(c(new.mat), c(V.modified)))
        }
        V.modified <- new.mat
    }

    if(do.DE.correction & !IsPositiveDefinite(V.modified)) {
        warning("Have to modify variance covariance matrix to make it positive definite, so results are approximate and the analysis will be slow.")
        new.mat <- AlterMatrixUsingDE(V.modified)
        likelihood.penalty <- 10+stats::dist(rbind(c(new.mat), c(V.modified)))
        V.modified <- new.mat
    }

    traits <- traits[match(names(means.modified), names(traits))] #reorder so traits in same order as VCV and means
    traits <- traits[!is.na(names(traits))]

    prune.taxa <- names(means.modified)[!names(means.modified)%in%names(traits)] # We might not have trait data for all tips. This is especially true if we want forward in time hybridization events, which are really a lineage branching off, living for a while, then moving genes to another lineage and then going extinct (or at least unsampled). We can allow a tree like that with just no trait values for the ghost taxon.

    if(nrow(V.modified) - length(prune.taxa) < 3) {
        stop(paste("You have pruned ", paste(prune.taxa, collapse=", "), " from your tree because they did not match taxa in your traits vector (where some of the names are ", paste(names(traits), collapse=", "), ") and this leaves too few taxa to do an analysis"))
    }
    if(length(prune.taxa)>0) {
        V.modified <- V.modified[-which(rownames(V.modified) %in% prune.taxa), -which(colnames(V.modified) %in% prune.taxa)]
        means.modified <- means.modified[-which(names(means.modified) %in% prune.taxa)]
    }



    NegLogML <- NULL
    try(NegLogML <- ((ncol(V.modified)/2)*log(2*pi)+(1/2)*t(traits-means.modified)%*%corpcor::pseudoinverse(V.modified)%*%(traits-means.modified) + (1/2)*determinant(V.modified, logarithm=TRUE)$modulus + likelihood.penalty)[1,1], silent=TRUE)
    if(is.null(NegLogML)) {
        NegLogML <- badval
    }
    return(NegLogML[1])
}

#' Optimize model
#'
#' Fits a BMhyb model to your data.
#'
#' This takes an ape::evonet object. If all you have is a tree (an ape::phylo object), you can use CreateHybridlessEvonet() to convert the tree to an evonet object. You can then use the AddHybridization() function to add hybrid events to this object. Note that networks created in this way can, by chance, result in orders of nodes in the internal edge matrix that cause ape's reorder.phylo function to crash, which is called in many of the plot and write functions. You can still use the plot functions in this package, however.
#'
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param traits A vector of trait values, with names equal to the names of taxa on the phylogeny
#' @param free.parameter.names What parameters you want to optimize rather than use defaults; options are sigma.sq, mu, SE, bt, and vh
#' @param confidence.points How many points to use to estimate parameter uncertainty
#' @param measurement.error How much uncertainty there is in tip values; a single number is applied to all taxa, a vector is applied to the corresponding taxa
#' @param gamma In a hybridization event, what proportion of the trait comes from the donating parent. 0.5 means half comes from each parent
#' @param do.Higham.correction Variance-covariance matrices for this model are sometimes poorly conditioned; this is a hack to reduce the impact of that
#' @param do.Brissette.correction Applies method of Brissette et al. 2007 to also try to fix matrix condition
#' @param verbose If TRUE, BMhyb will chat about its progress
#' @param likelihood.precision When optimizing, how much of a lnL improvement is required to restart optimization between starts
#' @param max.steps The number of restarts without improvement it will attempt
#' @param confidence.lnl For figuring out the confidence interval, how wide you want the confidence region to be in lnL space
#' @param control List of options to pass to optim. ?optim for help.
#'
#' @return Returns an object of class BMhybResult which contains best (a data.frame of the solution), good.region (data.frame of the points making up those in the confidence.lnl region), bad.region (all the other points sampled), phy.graph (same as what you put in), traits (same as what you put in), and free.parameter.names.
#'
#' @examples
#' \dontrun{
#' utils::data("cichlid")
#' result <- BMhyb(phy.graph=cichlid$phy.graph, traits=cichlid$trait,
#'   free.parameter.names=c("sigma.sq", "mu"))
#' }
#' @export
BMhyb <- function(phy.graph, traits, free.parameter.names=c("sigma.sq", "mu", "SE", "bt", "vh"), confidence.points = 5000, measurement.error=0, gamma=0.5, do.Higham.correction=FALSE, do.Brissette.correction=FALSE, verbose=TRUE, likelihood.precision=0.01, max.steps=10, confidence.lnl = 2, control=list(reltol=1e-3)) {
    if(verbose) {
      print("Now starting analysis")
    }
    best.results <- OptimizeThoroughly(phy.graph=phy.graph, traits=traits, free.parameter.names=free.parameter.names, measurement.error=measurement.error, gamma=gamma, do.Higham.correction=do.Higham.correction, do.Brissette.correction=do.Brissette.correction, max.steps=max.steps, verbose=verbose, control=control)

    local.df <- data.frame(t(best.results$par), AICc=ComputeAICc(n=ape::Ntip(phy.graph),k=length(best.results$par), LogLik=best.results$value),  NegLogLik=best.results$value, K=length(best.results$par))

    interval.results <- ComputeConfidenceIntervals(best.results$par, fn=TryComputeLikelihood, phy.graph=phy.graph, traits=traits, confidence.points=confidence.points,  best.lnl = best.results$value, likelihood.precision=likelihood.precision,  measurement.error=measurement.error, gamma=gamma, do.Higham.correction=do.Higham.correction, do.Brissette.correction=do.Brissette.correction )

    interval.results.in <- interval.results[which(interval.results[,1]-min(interval.results[,1])<=confidence.lnl),]
    interval.results.out <- interval.results[which(interval.results[,1]-min(interval.results[,1])>confidence.lnl),]
    result.object <- list(best=local.df, good.region=interval.results.in, bad.region=interval.results.out, phy.graph=phy.graph, traits=traits, free.parameter.names=free.parameter.names)
    class(result.object) <- "BMhybResult"
    return(result.object)
}

#' Exhaustively evaluate models
#'
#' Fits all possible BMhyb models to your data.
#'
#' This takes an ape::evonet object. If all you have is a tree (an ape::phylo object), you can use CreateHybridlessEvonet() to convert the tree to an evonet object. You can then use the AddHybridization() function to add hybrid events to this object. Note that networks created in this way can, by chance, result in orders of nodes in the internal edge matrix that cause ape's reorder.phylo function to crash, which is called in many of the plot and write functions. You can still use the plot functions in this package, however.
#'
#' This will return a list with one model result per element: you can plot these individually (see ?hybResult). By default, these results will include the information about uncertainty. We also compute a summary table so you can see the point estimates for each model and the likelihoods. It is often advisable to average across models, weighting each by its AICc weight, so this is also done automatically. We also return the single best model as an object for convenience, though for most users, we would suggest using the model average and looking at a set of fairly good models rather than look only at the single best one: there are often others that are nearly as good.
#'
#' We do not expect large AIC difference between models unless you have a really large tree, and so you may get a warning if this happens. It is likely something has gone wrong with optimization. Look at all the models and examine for outliers. This issue can come up with certain combinations of networks and parameters (even, very rarely, in Brownian motion with no hybridization), where a step in the likelihood (inverting a matrix) does not yield a numerically stable result (the matrix is poorly conditioned). The 'likelihoods' in such cases are wrong, and they can look too good or too bad. Neither is ideal, but you should especially beware cases where the 'best' model has likelihoods much below some of the other models -- you will often see bad parameter estimates, too. If you get this, do not believe the results -- perhaps look at models with better condition.
#'
#' To try to help with this, if one or more of the models has poor condition at the maximum likelihood estimate, we report this as it having an obvious problem. It is still returned in the results and the original.summary.df objects, but it is excluded from model averaging, the summary.df, and the best.model return (though note the ModelNumber column in summary.df allowing you to get the matching model in the results list). A model not having an obvious problem does *not* mean it worked well, just that it does not exhibit one particular problematic issue. Essentially we're saying, "This model does not have a lion eating its foot" -- which suggest it's not unhealthy in that way, but doesn't mean there's not a crocodile eating its hand. User beware. Plotting the confidence using the plot functions can help.
#'
#' @param phy.graph An ape::evonet object (a phylogeny stored in phylo format that also includes a reticulation matrix)
#' @param traits A vector of trait values, with names equal to the names of taxa on the phylogeny
#' @param measurement.error How much uncertainty there is in tip values; a single number is applied to all taxa, a vector is applied to the corresponding taxa
#' @param ncores Number of cores to use. By default, uses parallel package to detect what's available and uses all but one.
#' @param ... All other parameters to pass to BMhyb (see ?BMhyb)
#'
#' @return Returns a list of objects of class BMhybResult (results), a summary data frame (summary.df) with parameter estimates and weights for all models where we do not see obvious problems, a summary data frame of all the models, whether or no they seemed to fail (original.summary.df), the model averaged result weighted by AICc weights of the unproblematic models (model.average), and the best unproblematic model (best.model).
#'
#' @examples
#' \dontrun{
#' utils::data("cichlid")
#' traits.only <- cichlid$traits_and_SE$trait
#' names(traits.only) <- rownames(cichlid$traits_and_SE)
#' all.models <- BMhybExhaustive(phy.graph=cichlid$phy.graph, traits=traits.only)
#' print(all.models$summary.df)
#' }
#' @export
BMhybExhaustive <- function(phy.graph, traits, measurement.error=0, ncores=max(c(1, parallel::detectCores()-1), na.rm=TRUE),...) {
  results <- list()
  summary.df <- data.frame()
  free.parameter.matrix <- expand.grid(mu=TRUE, sigma.sq=TRUE, SE=c(TRUE, FALSE), bt=c(TRUE, FALSE), vh=c(TRUE, FALSE))

  matrix.conditions <- c()

  if(ncores>1) {
    print("Note that BMhyb will chattily print results as they go, but with parallel runs, the results from different models running in parallel will all be mixed together. You could pass verbose=FALSE to suppress much of the output")
  }

  do_single_model <- function(model.index, free.parameter.matrix, phy.graph, traits, measurement.error, ...) {
    free.parameter.row <- free.parameter.matrix[model.index,]
    free.parameters <- colnames(free.parameter.row)[unlist(free.parameter.row)]

    result <- BMhyb::BMhyb(phy.graph=phy.graph, traits=traits, free.parameter.names=free.parameters, measurement.error=measurement.error, ...)
    return(result)
  }

  results <- parallel::mclapply(sequence(nrow(free.parameter.matrix)), do_single_model, free.parameter.matrix=free.parameter.matrix, phy.graph=phy.graph, traits=traits, measurement.error=measurement.error, ..., mc.cores=ncores)

  # for (model.index in sequence(nrow(free.parameter.matrix))) {
  #   results[[model.index]] <- result
  # }

  for (model.index in sequence(nrow(free.parameter.matrix))) {

    summary.df <- plyr::rbind.fill(summary.df, results[[model.index]]$best)

    sigma.sq=1
    mu=0
    bt=1
    vh=0
    SE=0
    best.values <- unlist(results[[model.index]]$best[1,1:(length(results[[model.index]]$best)-3)])
    for(i in seq_along(best.values)){
        assign(names(best.values)[i],best.values[i])
    }

    final.VCV <- PruneDonorsRecipientsFromVCV(ComputeVCV(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, measurement.error=measurement.error))
    matrix.conditions <- c(matrix.conditions, kappa(final.VCV))
  }

  summary.df$MatrixCondition <- matrix.conditions

  summary.df$ObviousProblem <- (log(summary.df$MatrixCondition)>10)

  summary.df$ModelType <- paste0(ifelse(is.na(summary.df$sigma.sq), "_", "S"), ifelse(is.na(summary.df$mu), "_", "M"), ifelse(is.na(summary.df$vh), "_", "V"), ifelse(is.na(summary.df$bt), "_", "B"), ifelse(is.na(summary.df$SE), "_", "E"))


  summary.df$ModelNumber <- sequence(nrow(summary.df))

  original.summary.df <- summary.df

  summary.df <- summary.df[!summary.df$ObviousProblem,]


  model.average <- NULL
  best.model <- NULL

  if(nrow(summary.df)>0) {

    summary.df$deltaAICc <- summary.df$AICc - min(summary.df$AICc)
    rel.lik <- exp(-0.5* summary.df$deltaAICc)
    summary.df$AkaikeWeight <- rel.lik / sum(rel.lik)


    input.for.average <- summary.df[,-which(colnames(summary.df) %in% c("ModelNumber", "ModelType"))]
    input.for.average$vh[which(is.na(input.for.average$vh))] <- 0
    input.for.average$bt[which(is.na(input.for.average$bt))] <- 1
    input.for.average$SE[which(is.na(input.for.average$SE))] <- 0
    model.average <- data.frame(t(apply(input.for.average, 2, stats::weighted.mean, w=input.for.average$AkaikeWeight)))

    if(max(summary.df$deltaAICc)>100) {
      warning("Large differences in AICc values suggest something is wrong. See the documentation for this function and look at the MatrixConditions in the summary.df", immediate.=TRUE)
    }
    best.model <- results[[summary.df$ModelNumber[which(summary.df$deltaAICc==0)]]]

  } else {
    warning("no models seemed good (well conditioned) so the likelihoods and thus parameter estimates are unreliable. Try enabling some of the corrections (see ?BMhyb) to try to improve this", immediate.=TRUE)
  }




  return.object <- list(results=results, summary.df=summary.df, original.summary.df=original.summary.df, model.average=model.average, best.model=best.model)
  class(return.object) <- "BMhybExhaustiveResult"
  return(return.object)
}

#' Function to merge adaptive sampling sims for plotting
#'
#' @param exhaustive.object Return of BMhybExhaustive
#' @return Returns a single BMhyb object with results from all models merged (use for plotting)
#' @export
MergeExhaustiveForPlotting <- function(exhaustive.object) {
  #result.object <- list(best=local.df, good.region=interval.results.in, bad.region=interval.results.out, phy.graph=phy.graph, traits=traits, free.parameter.names=free.parameter.names)
  all.sims <- data.frame()
  best <- plyr::rbind.fill(lapply(exhaustive.object$results, '[[', 'best'))
  best$SE[which(is.na(best$SE))] <- 0
  best$bt[which(is.na(best$bt))] <- 1
  best$vh[which(is.na(best$vh))] <- 0
  best <- best[which.min(best$NegLogLik),]
  for (i in seq_along(exhaustive.object$results)) {
    all.sims <- plyr::rbind.fill(all.sims, exhaustive.object$results[[i]]$good.region)
    all.sims <- plyr::rbind.fill(all.sims, exhaustive.object$results[[i]]$bad.region)
  }
  all.sims$SE[which(is.na(all.sims$SE))] <- 0
  all.sims$bt[which(is.na(all.sims$bt))] <- 1
  all.sims$vh[which(is.na(all.sims$vh))] <- 0
  result.object <- list(
    best=unlist(best),
    good.region = all.sims[which(all.sims$negloglik<=(2+min(all.sims$negloglik))),],
    bad.region = all.sims[which(all.sims$negloglik>(2+min(all.sims$negloglik))),],
    phy.graph = exhaustive.object[[1]]$phy.graph,
    traits = exhaustive.object[[1]]$traits,
    free.parameter.names = colnames(all.sims)[-1]
  )
  class(result.object) <- "BMhybResult"
  return(result.object)
}

#' Get convex hull at a given threshold
#'
#' For a given delta lnL, get the convex hull (blob encircling the points) for two dimensions
#'
#' @param threshold What value to exclude numbers worse than
#' @param df The data.frame
#' @param height The variable name to use for the height threshold
#' @param x The first variable to look at for the hull
#' @param y The second variable to look at for the hull
#'
GetConvexHull <- function(threshold=2, df, height, x, y) {
  df <- df[df[,height]<= threshold,]
  df <- df[grDevices::chull(df[,x], df[,y]),]
  df$x <- df[,x]
  df$y <- df[,y]
  df$z <- df[,height]
  return(df)
}

#' Plot BMhyb result
#'
#' Shows the plot of confidence regions with MLEs indicated (red dots) or a plot of pairs of traits together. Note that for the latter plot, it converts the sampled points to an even grid with interpolation; it sets any points with likelihood worse than ten units to just ten units worse so that you can see the colors near the area of the optimum.
#'
#' @param x A BMhyb object (result of a BMhyb() call)
#' @param style Either univariate or contour
#' @param focal.color Color for the point showing the maximum likelihood estimate
#' @param inregion.color Color for univariate plot, points in the good region
#' @param outregion.color Color for univariate plot, points in the bad region
#' @param gradientworst.color Color for contour plot, color of the worst contour region
#' @param gradientbest.color Color for contour plot, color of the best contour region
#' @param contour.color Color showing the contour line for the best threshold
#' @param contour.threshold What delta log likelihood to use for the best/worst threshold for the contour plot
#' @param nrow The number of rows to plot in the grid for contour (will set it automatically if NULL)
#' @param ... Other arguments to pass to plot (for univariate only; the contour plot uses ggplot2)
#' @export
#' @rawNamespace S3method(plot, BMhybResult)
plot.BMhybResult <- function(x,style="univariate", focal.color="red", inregion.color = "black", outregion.color="gray", gradientworst.color="black", gradientbest.color="white", contour.color="red", contour.threshold=2, nrow=NULL, ...) {
    if(style=="univariate") {
      x$par <- x$best[1:(length(x$best)-3)]
      graphics::par(mfcol=c(1, length(x$par)))
      all.results <- rbind(x$good.region, x$bad.region)
      for(parameter in sequence(length(x$par))) {
          graphics::plot(x=all.results[,parameter+1], y=all.results[,1], type="n", xlab=names(x$par)[parameter], ylab="NegLnL", bty="n", ylim=c(min(all.results[,1]), min(all.results[,1])+10),...)
          graphics::points(x=x$good.region[,parameter+1], y=x$good.region[,1], pch=16, col=inregion.color)
          graphics::points(x=x$bad.region[,parameter+1], y=x$bad.region[,1], pch=16, col=outregion.color)
          graphics::points(x= x$best[parameter], y= x$best['NegLogLik'], pch=1, col=focal.color, cex=1.5)
      }
    } else {
      contour_plot <- function(data, x, y, x.best, y.best, breaks= c(1,2,5,10)) {
        data$delta_likelihood <- data$negloglik-min(data$negloglik)
        p <- ggplot2::ggplot(data, ggplot2::aes_string(x=x, y=y, z="delta_likelihood", fill="delta_likelihood")) + ggplot2::scale_fill_gradient(low=gradientbest.color, high=gradientworst.color, breaks=breaks, limits=c(min(breaks)-1, max(breaks)))
        #print(p)
        breaks <- sort(breaks, decreasing=TRUE)
        for (break_index in seq_along(breaks)) {
          polygon_points <- GetConvexHull(threshold=breaks[break_index], df=data, height="delta_likelihood", x=x, y=y)
          p <- p + ggplot2::geom_polygon(data=polygon_points, ggplot2::aes(x=x,y=y))
          #print(p)
        }
        polygon_points <- GetConvexHull(threshold=contour.threshold, df=data, height="delta_likelihood", x=x, y=y)
        p <- p + ggplot2::geom_polygon(data=polygon_points, ggplot2::aes(x=x,y=y), color=contour.color, fill=NA)
        p <- p + ggplot2::theme(legend.position = "none")
        p <- p+ ggplot2::geom_point(data=data, ggplot2::aes_string(x=x.best, y=y.best), colour=focal.color)
      #  print(p)
        return(p)
      }
      breaks <- sequence(10)
      x$par <- x$best[1:(length(x$best)-3)]
      all.results <- rbind(x$good.region, x$bad.region)
      #all.results <- x$good.region
      plotlist <- list()
      for(i in sequence(length(x$par))) {
        for (j in sequence(length(x$par))) {
          if (i < j) {
            plotlist[[length(plotlist)+1]] <- contour_plot(all.results, names(x$par)[i], names(x$par)[j], unname(unlist(x$par[i])[1]), unname(unlist(x$par[j])[1]), breaks=breaks)
            #print(contour_plot(all.results, names(x$par)[i], names(x$par)[j], unname(x$par[i][1,1]), unname(x$par[j][1,1])))
          }
        }
      }
      cowplot::plot_grid(plotlist=plotlist, nrow=nrow)






      # contour_plot <- function(data, x, y, x.best, y.best) {
      #   # using interpolation advice from https://stackoverflow.com/questions/35018971/3d-data-with-ggplot
      #   data$delta_likelihood <- data$negloglik-min(data$negloglik)
      #   im <- akima::interp(x=data[,x], y=data[,y], z=data[,"delta_likelihood"], nx = 500, ny = 500, extrap=TRUE, linear=FALSE)
      #   df2 <- data.frame(expand.grid(x = im$x, y = im$y), delta_likelihood = c(im$z))
      #   df2$delta_likelihood[which(df2$delta_likelihood>10)] <- 10
      #   colnames(df2) <- c(x, y, "delta_likelihood")
      #   return(ggplot2::ggplot(df2, ggplot2::aes_string(x=x, y=y, z="delta_likelihood", fill="delta_likelihood")) + metR::geom_contour_fill(ggplot2::aes(fill = ..level..)) +
      #   ggplot2::geom_contour(color = "red", size = .5, breaks=c(2)) + ggplot2::scale_fill_gradient(low="black", high="white", breaks=sequence(12)) + ggplot2::geom_point(ggplot2::aes_string(x=x.best, y=y.best), colour="red") + ggplot2::theme(legend.position = "none"))
      # }
      # x$par <- x$best[1:(length(x$best)-3)]
      # all.results <- rbind(x$good.region, x$bad.region)
      # #all.results <- x$good.region
      # plotlist <- list()
      # for(i in sequence(length(x$par))) {
      #   for (j in sequence(length(x$par))) {
      #     if (i < j) {
      #       plotlist[[length(plotlist)+1]] <- contour_plot(all.results, names(x$par)[i], names(x$par)[j], unname(unlist(x$par[i])[1]), unname(unlist(x$par[j])[1]))
      #       #print(contour_plot(all.results, names(x$par)[i], names(x$par)[j], unname(x$par[i][1,1]), unname(x$par[j][1,1])))
      #     }
      #   }
      # }
      # cowplot::plot_grid(plotlist=plotlist)
      # #return(cowplot::plot_grid(plotlist=plotlist))

    }
}

#' Print BMhyb result
#'
#' @param x A BMhyb object (result of a BMhyb() call)
#' @param ... Other arguments to pass to this function
#' @export
print.BMhybResult <- function(x, ...) {
  print(summary.BMhybResult(x,...))
}

#' Summarize BMhyb result
#'
#' @param object A BMhyb object (result of a BMhyb() call)
#' @param ... Other arguments to pass to this function
#' @return A data.frame with summarized results
#' @export
summary.BMhybResult <- function(object, ...) {
  best <- data.frame(matrix(object$best, nrow=1))
  colnames(best) <- names(object$best)
  min.good <- apply(object$good.region, 2, min)[-1]
  max.good <- apply(object$good.region, 2, max)[-1]
  results <- plyr::rbind.fill(best=best, lower=as.data.frame(t(min.good)), upper=as.data.frame(t(max.good)))
  rownames(results) <- c("MLE", "lower", "upper")
  return(results)
}

#' Print BMhybExhaustive result
#'
#' @param x A BMhybExhaustive object (result of a BMhybExhaustive() call)
#' @param ... Other arguments to pass to this function
#' @export
print.BMhybExhaustiveResult <- function(x, ...) {
  print(summary.BMhybExhaustiveResult(x, ...))
}

#' Summarize BMhybExhaustive result
#'
#' @param object A BMhybExhaustive object (result of a BMhybExhaustive() call)
#' @param ... Other arguments to pass to this function
#' @return A data.frame with summarized results
#' @export
summary.BMhybExhaustiveResult <- function(object, ...) {
  object$model.average$MatrixCondition <- NULL
  object$model.average$ObviousProblem <- NULL
  object$model.average$deltaAICc <- NULL
  object$model.average$AICc <- NULL
  object$model.average$NegLogLik <- NULL

  object$model.average$AkaikeWeight <- NULL
  merged <- MergeExhaustiveForPlotting(object)
  summary.merged <-  summary.BMhybResult(merged)
  results <- plyr::rbind.fill(as.data.frame((object$model.average)),summary.merged['lower',], summary.merged['upper',], as.data.frame((object$best.model$best)), object$summary.df)
  rownames(results) <- c("model.average", "lower", "upper", "best.model", paste0("model.", sequence(nrow(object$summary.df))))
  return(results)
}


#' Plot BMhybExhaustive result
#'
#' Note this aggregates the info from all the model runs and plots the MLE across all of them and the contours from all the sims
#'
#' @param x A BMhybExhaustive object (result of a BMhybExhaustive() call)
#' @param ... Other parameters to pass to hybResult
#' @export
plot.BMhybExhaustiveResult <- function(x,...) {
  x <- MergeExhaustiveForPlotting(x)
  plot.BMhybResult(x, ...)
}

# #' Plot evonet object
# #'
# #' ape's reorder.evonet function can case C crashes if the edge numbering is not in the manner ape expects (even if obvious requirements like terminal nodes being 1:Ntip and root node being Ntip+1 are met). Instead, we use igraph plotting.
# #' @param x An evonet object
# #' @param ... Other arguments to pass to plot
# #' @export
# #' @rawNamespace S3method(plot, evonet)
# plot.evonet <- function(x, size=0, shape="none", color=NA, frame.color=NA, vertex.size=0, vertex.color=NA, arrow.size=0.01, arrow.width=0.0, edge.arrow.size=0.1, ...) {
# #  phy.igraph <- BMhyb:::ConvertEvonetToIgraphWithNodeNumbers(x)
#   x$node.label <- seq(from=(1+ape::Ntip(x)), to=ape::Ntip(x)+ape::Nnode(x), by=1)
#   phy.igraph <- ape::as.igraph.evonet(x)
#   plot(phy.igraph, layout=igraph::layout_as_tree(phy.igraph, circular=TRUE), size=size, shape=shape, color=color, frame.color=frame.color, vertex.size=vertex.size, vertex.color=vertex.color, arrow.size=arrow.size, arrow.width=arrow.width, edge.arrow.size=edge.arrow.size, ...)
# }


OptimizeThoroughly <- function(phy.graph, traits, free.parameter.names=c("sigma.sq", "mu", "SE", "bt", "vh"), measurement.error=0, gamma=0.5, do.Higham.correction=FALSE, do.Brissette.correction=FALSE, do.DE.correction=FALSE, verbose=TRUE, likelihood.precision=0.01, max.steps=10, control=list()) {
    simple.phy <- ape::collapse.singles(ape::as.phylo(phy.graph))
    cleaned <- geiger::treedata(simple.phy, traits, warnings=FALSE, sort=TRUE)
    starting.from.geiger <- geiger::fitContinuous(cleaned$phy, cleaned$data, model="BM", SE=0, ncores=1)$opt
    starting.values <- c(sigma.sq=starting.from.geiger$sigsq, mu=starting.from.geiger$z0, bt=1,  vh=0.01*starting.from.geiger$sigsq*max(ape::vcv(simple.phy)), SE=0.01*starting.from.geiger$sigsq*max(ape::vcv(simple.phy))) #sigma.sq, mu, beta, vh, SE
    starting.values <- starting.values[free.parameter.names]

    # now test initial values to make sure not in a bad part of matrix condition
    sigma.sq=1
    mu=0
    bt=1
    vh=0
    SE=0
    for(i in seq_along(starting.values)){
        assign(names(starting.values)[i],starting.values[i])
    }

    starting.VCV <- PruneDonorsRecipientsFromVCV(ComputeVCV(phy.graph, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, measurement.error=measurement.error))
    original.starting.values <- starting.values
    while(log(kappa(starting.VCV))>15) { # really hard to get accurate likelihood scores here
      if(verbose) {
        print("Initial starting values were in an area where the precision of the likelihood would be low (due to poor condition of the VCV matrix), so we're trying in a new area")
      }
      starting.values <- stats::rnorm(length(starting.values), mean=starting.values, sd=abs(starting.values)/10)
      if(stats::runif(1) < 0.1) { # lest we wander too far from good values, let's restart periodically near the original estimates
        starting.values <- stats::rnorm(length(original.starting.values), mean=starting.values, sd=abs(starting.values)/10)
      }
      if(stats::runif(1) < 0.1) { # maybe try some more normal values, too.
        starting.values <- original.starting.values
        if("vh" %in% names(starting.values)) {
          starting.values$vh <- 0
        }
        if("bt" %in% names(starting.values)) {
          starting.values$bt <- 1
        }
        if("SE" %in% names(starting.values)) {
          starting.values$SE <- 0
        }
      }
      for (i in seq_along(starting.values)){
          if(any(grepl(names(starting.values)[i], c("sigma.sq", "vh", "bt", "SE")))) {
            starting.values[i] <- abs(starting.values[i])
          }
      }
    }



    if(verbose) {
      print("Starting with initial values")
      print(starting.values)
    }
    best.run <- stats::optim(par=starting.values, fn=TryComputeLikelihood, traits=traits, phy.graph=phy.graph, measurement.error=measurement.error, gamma=gamma, do.Higham.correction=do.Higham.correction, do.Brissette.correction=do.Brissette.correction, control=control)
    attempts <- 1
    step.count <- 1
    while(best.run$convergence!=0 && attempts < 10){#want to get a convergence code 0
        if(verbose) {
            print(paste0("Initial search had a convergence code of ", best.run$convergence, ", indicating it did not converge. See ?optim for what the code may mean. Starting again, likely near that point. Negative log likelihood was ", best.run$value))
            print("Parameter estimates were")
            print(best.run$par)
        }
        new.starting.values <- stats::runif(length(starting.values), min=starting.values - attempts*.1*starting.values, max=starting.values + attempts*.1*starting.values)
        names(new.starting.values) <- names(starting.values)
        best.run <- stats::optim(par=new.starting.values, fn=TryComputeLikelihood, traits=traits, phy.graph=phy.graph, measurement.error=measurement.error, gamma=gamma, do.Higham.correction=do.Higham.correction, do.Brissette.correction=do.Brissette.correction, control=control)
        attempts <- attempts+1
    }
    if(verbose) {
        results.vector<-c(step.count, best.run$value, best.run$par)
        names(results.vector) <- c("step","negloglik", names(best.run$par))
        print(results.vector)
    }

    #this is to continue optimizing; we find that optim is too lenient about when it accepts convergence

    times.without.improvement <- 0
    starting.values <- best.run$par
    while(times.without.improvement<max.steps) {
        times.without.improvement <- times.without.improvement+1
        step.count <- step.count + 1
        new.run <- stats::optim(par=starting.values, fn=TryComputeLikelihood, traits=traits, phy.graph=phy.graph, measurement.error=measurement.error, gamma=gamma, do.Higham.correction=do.Higham.correction, do.Brissette.correction=do.Brissette.correction, control=control)
        new.starting.values <- stats::runif(length(starting.values), min=starting.values - attempts*.1*starting.values, max=starting.values + attempts*.1*starting.values)
        names(new.starting.values) <- names(starting.values)
        starting.values <- new.starting.values
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
            results.vector<-c(step.count, best.run$value, best.run$par)
            names(results.vector) <- c("step","negloglik", names(best.run$par))
            print(results.vector)
        }
    }
    return(best.run)
}

## Stuff below is from aborted attempt to use network. Trying again.

# PlotPhyGraph <- function(phy.graph) {
#     # TODO: get brlen, suppress internal node names
#     plot(phy.graph, layout = layout_as_tree, vertex.size=5)
# }

# ConvertPhyWithFlowIntoMultiphy <- function(phy, flow) {
#     phy.graph <- as.igraph(multiphy[[1]], directed=TRUE)
#     phy.graph <- set_edge_attr(phy.graph, "length", value=multiphy[[1]]$edge.length)
#     phy.graph <- set_edge_attr(phy.graph, "weight", value=1)
#     return(phy.graph)
# }

# AddNodeLabels <- function(phy) {
#     nodes <- seq(from=ape::Ntip(phy)+1, length.out=ape::Nnode(phy))
#     for (node.index in sequence(length(nodes))) {
#         phy$node.label[node.index] <- paste0("InternalNode_",paste(sort(ape::extract.clade(phy, nodes[node.index])$tip.label), collapse="_"))
#     }
#     return(phy)
# }

# MergeTreesIntoPhyGraph<- function(multiphy) {
#     stop("This function is not good. It creates too many input edges to non-hybrid nodes")
#     phy.graph <- as.igraph(AddNodeLabels(multiphy[[1]]), directed=TRUE)
#     phy.graph <- set_edge_attr(phy.graph, "length", value=multiphy[[1]]$edge.length)
#     phy.graph <- set_edge_attr(phy.graph, "weight", value=1)
#     if(length(multiphy)>1) {
#         for (i in 2:length(multiphy)) {
#             phy.graph2 <- as.igraph(AddNodeLabels(multiphy[[i]]), directed=TRUE)
#             phy.graph2 <- set_edge_attr(phy.graph2, "length", value=multiphy[[i]]$edge.length)
#             phy.graph2 <- set_edge_attr(phy.graph2, "weight", value=1)
#
#             phy.graph <- phy.graph  %u% phy.graph2
#         }
#     }
#     edge.properties <- data.frame(edge_attr(phy.graph))
#     all.weights <- edge.properties[,grepl("weight", colnames(edge.properties))]
#     sum.weights <- rowSums(all.weights, na.rm=TRUE)
#     phy.graph <- set_edge_attr(phy.graph, "overall_weight", value=sum.weights)
#     return(phy.graph)
# }

# ExhaustiveV <- function(x, phylogenies, phylogeny.weights, actual.params, measurement.error=NULL, drop.internal=TRUE) {
#     bt <- 1
#     vh <- 0
#     sigma.sq <- x[1]
#     mu <- x[2]
#     SE <- 0
#     if(is.null(measurement.error)) {
#         SE <- x[length(x)]
#     }
#     bt.location <- which(names(actual.params)=="bt")
#     if(length(bt.location)==1) {
#         bt<-x[bt.location]
#     }
#     vh.location <- which(names(actual.params)=="vh")
#     if(length(vh.location)==1) {
#         vh<-x[vh.location]
#     }
#
#     all.V.matrices <- sigma.sq*lapply(phylogenies, ape::vcv)
#     V.matrix <- 0*all.V.matrices[[1]]
#     mean.vector <- rep(0, length(igraph::V(phy.graph)))
#     for (tree.rep in sequence(length(phylogenies))) {
#         for (i in sequence(nrow(V.matrix))) {
#             row.taxon <- rownames(V.matrix)[i]
#             for(j in sequence(ncol(V.matrix))) {
#                 col.taxon <- colnames(V.matrix)[j]
#                 V.matrix[i,j] <- V.matrix[i,j] + (phylogeny.weights[tree.rep]^2)*all.V.matrices[[tree.rep]][row.taxon, col.taxon]
#             }
#         }
#     }
# }

# GetVandMFromIgraph <- function(x, phy.graph, actual.params, measurement.error=NULL, drop.internal=TRUE) {
#     bt <- 1
#     vh <- 0
#     sigma.sq <- x[1]
#     mu <- x[2]
#     SE <- 0
#     if(is.null(measurement.error)) {
#         SE <- x[length(x)]
#     }
#     bt.location <- which(names(actual.params)=="bt")
#     if(length(bt.location)==1) {
#         bt<-x[bt.location]
#     }
#     vh.location <- which(names(actual.params)=="vh")
#     if(length(vh.location)==1) {
#         vh<-x[vh.location]
#     }
#
#     V.matrix <- matrix(0, nrow=length(igraph::V(phy.graph)), ncol=length(igraph::V(phy.graph)))
#     mean.vector <- rep(0, length(igraph::V(phy.graph)))
#     postorder.traversal <- igraph::dfs(phy.graph, names(V(phy.graph))[which.max(nchar(names(V(phy.graph))))])$order
#
#     rownames(V.matrix) <- names(postorder.traversal)
#     colnames(V.matrix) <- names(postorder.traversal)
#     names(mean.vector) <- names(postorder.traversal)
#     parents <- igraph::adjacent_vertices(phy.graph, names(postorder.traversal), mode="in")
#     all.attribute.names <- names(igraph::get.edge.attribute(phy.graph))
#     matching.attribute.names <- all.attribute.names[grepl("length_",all.attribute.names )]
#     all.lengths <- matrix(nrow=length(matching.attribute.names), ncol=length(igraph::E(phy.graph)))
#
#     matching.weight.names <- all.attribute.names[grepl("weight_",all.attribute.names )]
#
#     all.weights <- matrix(0, nrow=length(matching.weight.names), ncol=length(igraph::E(phy.graph)))
#     for (length.index in sequence(length(matching.attribute.names))) {
#         all.lengths[length.index,] <- get.edge.attribute(phy.graph, name=matching.attribute.names[length.index])
#     }
#
#     for (weight.index in sequence(length(matching.weight.names))) {
#         all.weights[weight.index,] <- get.edge.attribute(phy.graph, name=matching.weight.names[weight.index])
#     }
#
#     all.weights <- sweep(all.weights,MARGIN=2,FUN="/",STATS=colSums(all.weights, na.rm=TRUE))
#     all.weights[is.na(all.weights)] <- 0
#     all.lengths[is.na(all.lengths)] <- 0
#
#
#     # Algorithm from Bastide et al., Syst Biol. 2018 in press
#     for (focal.index in sequence(length(names(postorder.traversal)))) {
#         if(focal.index==1) {
#             mean.vector[1] <- mu
#         } else {
#             if(length(parents[[focal.index]])==1) { #tree node
#                 focal.node <- names(postorder.traversal)[focal.index]
#                 parent.node <- names(parents[[focal.index]])
#                 for (other.node.index in sequence(focal.index-1)) { #Bastide et al. eq 3
#                     other.node <- names(postorder.traversal)[other.node.index]
#                     V.matrix[focal.node, other.node] <- V.matrix[parent.node, other.node]
#                     V.matrix[other.node, focal.node] <- V.matrix[focal.node, other.node] #do the upper and lower tri
#                 }
#                 focal.edge <- get.edge.ids(phy.graph, c(parent.node, focal.node))
#                 if(any(is.na(V.matrix))) {
#                     print("841")
#                     stop()
#                 }
#                 V.matrix[focal.node, focal.node] <- V.matrix[parent.node, parent.node] + sigma.sq * weighted.mean(all.lengths[,focal.edge], all.weights[,focal.edge], na.rm=TRUE)
#                 mean.vector[focal.node] <- mean.vector[parent.node]
#                 if(any(is.na(V.matrix))) {
#                     print("847")
#                     stop()
#                 }
#             } else { #hybrid node
#                 focal.node <- names(postorder.traversal)[focal.index]
#                 parent.nodes <- names(parents[[focal.index]])
#                 for (other.node.index in sequence(focal.index-1)) { #Bastide et al. eq 3
#                     other.node <- names(postorder.traversal)[other.node.index]
#                     V.matrix[focal.node, other.node] <- 0
#                     for (parent.index in sequence(length(parent.nodes))) {
#                         focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[parent.index], focal.node))
#                         V.matrix[focal.node, other.node] <- V.matrix[focal.node, other.node] + sum(0, all.weights[parent.index, focal.edge] * V.matrix[parent.nodes[parent.index], other.node], na.rm=TRUE)
#                         if(any(is.na(V.matrix))) {
#                             print("860")
#                             stop()
#                         }
#                     }
#                     V.matrix[other.node, focal.node] <- V.matrix[focal.node, other.node] #do the upper and lower tri
#                 }
#
#                 #now for Vii
#                 V.matrix[focal.node, focal.node] <- 0
#                 focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[1], focal.node))
#
#                 for (parent.index in sequence(length(parent.nodes))) {
#                     focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[parent.index], focal.node))
#                     V.matrix[focal.node, focal.node] <- V.matrix[focal.node, focal.node] + (all.weights[parent.index, focal.edge]^2) * (sum(0,V.matrix[parent.nodes[parent.index], parent.nodes[parent.index]], na.rm=TRUE) + sigma.sq * all.lengths[parent.index,focal.edge])
#                     if(any(is.na(V.matrix))) {
#                         print("876")
#                         stop()
#                     }
#                 }
#
#                 if(length(parent.nodes)>2) {
#                     stop("This code, from Bastide et al. eq 4, only envisions two parents for a node")
#                 }
#                 V.matrix[focal.node, focal.node] <- sum(vh, V.matrix[focal.node, focal.node], 2*all.weights[1, focal.edge]*all.weights[2, focal.edge]*V.matrix[parent.nodes[1], parent.nodes[2]], na.rm=TRUE)
#                 mean.vector[focal.node] <- sum(log(bt),  mean.vector[parent.nodes[1]]*all.weights[1, focal.edge], mean.vector[parent.nodes[2]]*all.weights[2, focal.edge], na.rm=TRUE)
#                 if(any(is.na(V.matrix))) {
#                     print("886")
#                     stop()
#                 }
#             }
#         }
#
#     }
#     if(drop.internal) {
#         elements.to.keep <- !grepl("Node_", names(mean.vector))
#         mean.vector <- mean.vector[elements.to.keep]
#         elements.to.keep <- !grepl("Node_", rownames(V.matrix))
#         V.matrix <- V.matrix[elements.to.keep,elements.to.keep]
#     }
#     diag(V.matrix) <- diag(V.matrix) + SE^2
#     return(list(V.modified=V.matrix, means.modified=mean.vector))
# }
#
# ExpandEvonet <- function(phy.graph) {
#     phy.graph <- reorder(phy.graph, "cladewise")
#     phy.graph$edge.network <- rbind(phy.graph$edge, phy.graph$reticulation)
#     reticulation.lengths <- rep(NA, nrow(phy.graph$reticulation))
#     for (i in sequence(nrow(phy.graph$reticulation))) {
#         reticulation.lengths[i] <- abs(phytools::nodeheight(phy.graph, phy.graph$reticulation[i,1])-phytools::nodeheight(phy.graph, phy.graph$reticulation[i,2]))
#     }
#     phy.graph$edge.length.network <- c(phy.graph$edge.length, reticulation.lengths)
#     return(phy.graph)
# }

# GetVandMFromExpandedEvonet <- function(x, phy.graph, actual.params, measurement.error=NULL, drop.internal=TRUE) {
#     bt <- 1
#     vh <- 0
#     sigma.sq <- x[1]
#     mu <- x[2]
#     SE <- 0
#     if(!("edge.length.network" %in% names(phy.graph))) {
#         phy.graph <- ExpandEvonet(phy.graph)
#     }
#     if(is.null(measurement.error)) {
#         SE <- x[length(x)]
#     }
#     bt.location <- which(names(actual.params)=="bt")
#     if(length(bt.location)==1) {
#         bt<-x[bt.location]
#     }
#     vh.location <- which(names(actual.params)=="vh")
#     if(length(vh.location)==1) {
#         vh<-x[vh.location]
#     }
#     N.all.nodes <- ape::Nnode(phy.graph)+ape::Ntip(phy.graph)
#     V.matrix <- matrix(0, nrow=N.all.nodes, ncol=N.all.nodes)
#     mean.vector <- rep(0, N.all.nodes)
#     rownames(V.matrix) <- colnames(V.matrix) <- names(mean.vector) <- sequence(N.all.nodes)
#
#     ##################
#     # Start adding here. Go up the tree using the $edge.network, top to bottom, using $edge.length.network for edge lengths
#
#     # Algorithm from Bastide et al., Syst Biol. 2018 in press
#     for (focal.index in sequence(length(names(postorder.traversal)))) {
#         if(focal.index==1) {
#             mean.vector[1] <- mu
#         } else {
#             if(length(parents[[focal.index]])==1) { #tree node
#                 focal.node <- names(postorder.traversal)[focal.index]
#                 parent.node <- names(parents[[focal.index]])
#                 for (other.node.index in sequence(focal.index-1)) { #Bastide et al. eq 3
#                     other.node <- names(postorder.traversal)[other.node.index]
#                     V.matrix[focal.node, other.node] <- V.matrix[parent.node, other.node]
#                     V.matrix[other.node, focal.node] <- V.matrix[focal.node, other.node] #do the upper and lower tri
#                 }
#                 focal.edge <- get.edge.ids(phy.graph, c(parent.node, focal.node))
#                 if(any(is.na(V.matrix))) {
#                     print("841")
#                     stop()
#                 }
#                 V.matrix[focal.node, focal.node] <- V.matrix[parent.node, parent.node] + sigma.sq * weighted.mean(all.lengths[,focal.edge], all.weights[,focal.edge], na.rm=TRUE)
#                 mean.vector[focal.node] <- mean.vector[parent.node]
#                 if(any(is.na(V.matrix))) {
#                     print("847")
#                     stop()
#                 }
#             } else { #hybrid node
#                 focal.node <- names(postorder.traversal)[focal.index]
#                 parent.nodes <- names(parents[[focal.index]])
#                 for (other.node.index in sequence(focal.index-1)) { #Bastide et al. eq 3
#                     other.node <- names(postorder.traversal)[other.node.index]
#                     V.matrix[focal.node, other.node] <- 0
#                     for (parent.index in sequence(length(parent.nodes))) {
#                         focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[parent.index], focal.node))
#                         V.matrix[focal.node, other.node] <- V.matrix[focal.node, other.node] + sum(0, all.weights[parent.index, focal.edge] * V.matrix[parent.nodes[parent.index], other.node], na.rm=TRUE)
#                         if(any(is.na(V.matrix))) {
#                             print("860")
#                             stop()
#                         }
#                     }
#                     V.matrix[other.node, focal.node] <- V.matrix[focal.node, other.node] #do the upper and lower tri
#                 }
#
#                 #now for Vii
#                 V.matrix[focal.node, focal.node] <- 0
#                 focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[1], focal.node))
#
#                 for (parent.index in sequence(length(parent.nodes))) {
#                     focal.edge <- get.edge.ids(phy.graph, c(parent.nodes[parent.index], focal.node))
#                     V.matrix[focal.node, focal.node] <- V.matrix[focal.node, focal.node] + (all.weights[parent.index, focal.edge]^2) * (sum(0,V.matrix[parent.nodes[parent.index], parent.nodes[parent.index]], na.rm=TRUE) + sigma.sq * all.lengths[parent.index,focal.edge])
#                     if(any(is.na(V.matrix))) {
#                         print("876")
#                         stop()
#                     }
#                 }
#
#                 if(length(parent.nodes)>2) {
#                     stop("This code, from Bastide et al. eq 4, only envisions two parents for a node")
#                 }
#                 V.matrix[focal.node, focal.node] <- sum(vh, V.matrix[focal.node, focal.node], 2*all.weights[1, focal.edge]*all.weights[2, focal.edge]*V.matrix[parent.nodes[1], parent.nodes[2]], na.rm=TRUE)
#                 mean.vector[focal.node] <- sum(log(bt),  mean.vector[parent.nodes[1]]*all.weights[1, focal.edge], mean.vector[parent.nodes[2]]*all.weights[2, focal.edge], na.rm=TRUE)
#                 if(any(is.na(V.matrix))) {
#                     print("886")
#                     stop()
#                 }
#             }
#         }
#
#     }
#     if(drop.internal) {
#         elements.to.keep <- !grepl("Node_", names(mean.vector))
#         mean.vector <- mean.vector[elements.to.keep]
#         elements.to.keep <- !grepl("Node_", rownames(V.matrix))
#         V.matrix <- V.matrix[elements.to.keep,elements.to.keep]
#     }
#     diag(V.matrix) <- diag(V.matrix) + SE^2
#     return(list(V.modified=V.matrix, means.modified=mean.vector))
# }
