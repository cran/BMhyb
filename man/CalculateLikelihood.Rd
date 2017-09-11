\name{CalculateLikelihood}
\alias{CalculateLikelihood}
\title{
	Calculate the likelihood value for the model
	}

\description{
	This function calculates the likelihood value for the model.
	}
\usage{
	CalculateLikelihood(x, data, phy, flow, actual.params, precision=2,
	proportion.mix.with.diag=0, allow.extrapolation=FALSE,
	measurement.error = NULL, do.kappa.check = FALSE, number.of.proportions=101,
	lower.b=c(0, -Inf, 0.000001, 0, 0), upper.b=c(10,Inf,100,100,100),
	badval.if.not.positive.definite=TRUE, do.Brissette.correction=FALSE, 
	do.Higham.correction=TRUE, do.DE.correction=FALSE, return.penalty=FALSE, ...)
	}

\arguments{
  \item{x}{parameter of interest.}
  \item{data}{the trait values.}
  \item{phy}{a tree of phylo class.}
  \item{flow}{strcuture of gene flow.}
  \item{actual.params}{the free parameters for hybridization.}
  \item{precision}{a number to verify the condition of the variance covariation for the network model.}
  \item{proportion.mix.with.diag}{the proportion value that applies to the diagonal of the vcv matrix.}
	\item{allow.extrapolation}{a TRUE/FALSE argument.}
	\item{measurement.error}{estimate or fixed (see details).}
	\item{do.kappa.check}{check for VCV matrix condition, approximate if poor. Leave TRUE unless you're cleverer than us.}
  \item{number.of.proportions}{control how many points to use for extrapolation for ill-conditioned matrices.}
	\item{lower.b}{minimum values of all possible parameters (non-free ones will be deleted).}
	\item{upper.b}{maximum values of all possible parameters (non-free ones will be deleted).}
	\item{badval.if.not.positive.definite}{If the VCV is not positive definite, return a large bad value for likelihood.}
	\item{do.Brissette.correction}{If VCV is not positive definite, do Brissette et al. (2007) correction.}
	\item{do.Higham.correction}{If VCV is not positive definite, do Higham (2002).}
	\item{do.DE.correction}{If VCV is not positive definite, find a nearby matrix that is.}
	\item{return.penalty}{If TRUE, return a penalty for approximating not the likelihood itself.}
	\item{...}{place to absorb other arguments.}
}

\details{
	This function calculates the likelihood value. As described in the argument, the \eqn{x} contains the parameter of interests including over all mean \eqn{\mu}, the rate of evolution \eqn{\sigma^2}, the measurement error SE and the hybrid vigor \eqn{\beta} and the variation at the burst of hybridization \eqn{v_H}. Names of the taxa in data vector and the phylogeny must be match for further analysis, otherwise it would terminate immediately. The negative log likelihood function is calculated and a reasonable likelihood value is returned by checking the variance covariance matrix where a precision number is used for verifying the log conditioned number of the variance covariance of the network model. When the matrix is ill conditioned, we modestly adjust the matrix by the shrink the value of the off diagonal matrix using a proportion array. Then the likelihood is calculated using a spline approximation.
	}

\value{
	negative log likelihood value.
	 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#set the number of hybrid
	ntax.nonhybrid<-5
	#set the number of hybrid
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
sigma.sq <- 0.01;mu <- 1;SE <- 0
	#simulatedata
	data<-rnorm(ntax.nonhybrid+ntax.hybrid)
	names(data)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")
	#calculate the likelihood value
	CalculateLikelihood(c(sigma.sq,mu,SE), data, network$phy, network$flow, "vh",
	precision=2, proportion.mix.with.diag=0, allow.extrapolation=TRUE)
	}
