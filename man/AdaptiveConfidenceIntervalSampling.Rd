\name{AdaptiveConfidenceIntervalSampling}
\alias{AdaptiveConfidenceIntervalSampling}
\title{
Confidence interval under adaptive cluster sampling technique.
}
\description{
This function uses an adaptive cluster sampling technique to generate confidence interval for parameters of interest.
}
\usage{
AdaptiveConfidenceIntervalSampling(par, fn, lower=-Inf, upper=Inf, desired.delta = 2,
n.points=5000, verbose=TRUE,  measurement.error=NULL, do.kappa.check=FALSE,
allow.restart=TRUE,  best.lnl = -Inf, likelihood.precision=0.001, restart.mode=FALSE, ...)
}

\arguments{
  \item{par}{parameter of interest for sampling in the model.}
  \item{fn}{the negative log-likelihood function.}
  \item{lower}{the lower bound for the values allowed for sampling.}
  \item{upper}{the upper bound for the values allowed for sampling.}
  \item{desired.delta}{a numeric value with a default value of 2 for the criteria that the desired log-likelihood shall no more than 2 unit away from the maximum.}
  \item{n.points}{number of required points for calculating likelihood value.}
  \item{verbose}{whether to print detailed information during the run.}
  \item{measurement.error}{single value or vector of fixed measurement errors.}
  \item{do.kappa.check}{if TRUE, check for matrix condition when calculating likelihood.}
  \item{allow.restart}{if TRUE, stop and go back to calling function if this finds a better value}
  \item{best.lnl}{the value of the best likelihood from the original run. Only relevant if allow.restart=TRUE.}
  \item{likelihood.precision}{a numerical value used for verifying the convergence of the estimation.}
  \item{restart.mode}{restart if better values found}
  \item{...}{further arguments passed(see details).}
}
\details{
This function starts with a set of parameter values, generates new points using the function \code{GeneratedValues}, and the likelihood value is calculated. Following an idea by Edwards (1992), a support region is based on points within a certain likelihood score of the best one (by default, 2 log likelihood units). The function tries to look over a region wide enough to encompass parameter values that are in this range. Note this is done for all traits at once. This generates wider (more conservative) confidence regions than looking at each parameter value separately. For example, if there were a ridge such that two values can covary without changing the likelihood, a univariate procedure that changes one but holds the others constant would give a narrow region, while our approach of trying many points would give a much wider region.
}

\value{
 a data frame where the first column contains the calculated likelihood values and the rest of columns are the grid points of the parameters generated under the adaptive sampling technique.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
  Edwards, A.W.F. 1993. Likelihood. Johns Hopkins University Press, Baltimore

  Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted.

	Seber G.A.F., Salehi M.M.  2013. Adaptive sampling designs: Inference for spatial and clustered population.  Springer.
  }

\examples{
	library(corpcor)
	#assign the number of non hybrid taxa
	ntax.nonhybrid<-8
	#assign the number of hybrid
	ntax.hybrid<-3
	#simulate the network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 3, allow.ghost=FALSE)

	#generate the tip data
	parameters <- c(0.05,1,1.7,0.005,0.5)
	names(parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	tips <- SimulateTipData(network$phy, network$flow, parameters)

	#set free parameters
	free.parameters<-rep(TRUE, 5)
	names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")


	#Simulate 100 samples
	interval.results <- AdaptiveConfidenceIntervalSampling(parameters, fn=CalculateLikelihood,
	lower=c(0, -Inf, 0, 0, 0)[which(free.parameters)], n.points=100,data=tips,
	phy=network$phy, flow=network$flow, actual.params =
	free.parameters[which(free.parameters)], allow.extrapolation=TRUE)

  #show the results
	interval.results.in <- interval.results[which(interval.results[,1] -
	min(interval.results[,1])<=2),]
	interval.results.out <- interval.results[which(interval.results[,1] -
	min(interval.results[,1])>2),]
	interval.results.in
	interval.results.out
	}
