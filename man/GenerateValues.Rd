\name{GenerateValues}
\alias{GenerateValues}
\title{
	Generate parameter values
	}

\description{
	This function generates new parameter values in a range.
	}

\usage{
	GenerateValues(par, lower, upper, max.tries=100, expand.prob=0, examined.max,
	examined.min)
	}

\arguments{
  \item{par}{parameters of interest.}
  \item{lower}{the actual lower bound for the parameters.}
  \item{upper}{the actual upper bound for the parameters.}
  \item{max.tries}{maximum number of attempts to generate the parameter value. The default is set to 100.}
  \item{expand.prob}{a probability value to expand.}
  \item{examined.max}{the allowed maximum for the parameters values.}
  \item{examined.min}{the allowed minimum for the parameters values}
  }

\details{
	This function generates new parameter values using the uniform distribution with the interval\eqn{ (L, U)}  where for each parameter of interest, the lower bound of the interval is set to  \eqn{ L=\max (lower, 0.9*examined.min )} and the upper bound of the interval is set to \eqn{ U =\min (upper, 1.1*examined.max )}.
	}

\value{
	The generated parameter values.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted.

	Seber G.A.F., Salehi M.M.  2013. Adaptive sampling designs: Inference for spatial and clustered population.  Springer.
	}

\examples{
	#parameter of interests
	mu<-1;sigma.sq<-0.1;bt<-1;vH<-0;SE<-0
	#parameters of interest
	par<-c(mu,sigma.sq,bt,vH,SE)
	#the lower bound
	lower=c(-2,0,-10,0,0)
	#the upper bound
	upper<-c(2,1,10,100,100)
	#the examined minimum value
	examined.min<-c(-1.8,0.001,-9.8,0,0)
	#the examined maximum value
	examined.max<-c(1.90,1.1,10.2,100,100)
	#simulate points
	sim.points<-GenerateValues(par, lower, upper, max.tries=100, expand.prob=0,
	examined.max, examined.min)
	names(sim.points)<-c("mu","sigma.sq","bt","vH","SE")
	#show result
	sim.points
	}
