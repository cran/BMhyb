\name{AICc}
\alias{AICc}
\title{
	AICc
	}

\description{
	Calculates the second order Akaike's information criterion score for models of interest.
	}

\usage{
	AICc(n, k, LogLik)
	}

\arguments{
  \item{n}{
		number of taxa for the given phylogenetic tree. It represents the sample size(the number of species on the tip of phylogeny).
	}
  \item{k}{
	number of free parameters in the model.
	}
  \item{LogLik}{
   the minimum of the negative log-likelihood value obtained by optimizing the likelihood function.
	}
	}

\details{
	'AICc' is a function to compute the AICc values and is valid to select among different models. \eqn{AICc = 2*n*k/(n-k-1) -2 log L} where \eqn{L} is the maximum likelihood for the model.
	}

\value{
	The AICc values.
 	}

\author{
	Dwueng-Chwuan Jhwueng
	}

\references{
	Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.
	}
\examples{
  #assign the size 
	n<-5
	#assign the number of parameter
	k<-3
	#assign the negative log likelihood value.
	LogLik<- -2
	#compute the AICc score
	AICc(n,k,LogLik)
	# result AICc value of 26.
	}
