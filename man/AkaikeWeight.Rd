\name{AkaikeWeight}
\alias{AkaikeWeight}
\title{
	Akaike weight
	}

\description{
	Calculate Akaike weight for the models.
	}

\usage{
	AkaikeWeight(Delta.AICc.Array)
	}

\arguments{
  \item{Delta.AICc.Array}{
	Delta.AICc.Array is defined as the difference between the AICc value and the minimum AICc value among candidate models.
 	}
	}

\details{
	For \eqn{n} models of interest, the Akaikei weight for the \eqn{i} th model is defined as \eqn{w = exp(-0.5* \Delta AICc_i)/ sum_i dAICc_i} where \eqn{\Delta AICc_i} is the difference of the AICc value between the the \eqn{i}th model and the best model. The weights can be used in model averaging in advanced.
	}

\value{
 	Akaike weights.
 	}

\author{
	Dwueng-Chwuan Jhwueng
	}

\references{
	Burham, K.P., Anderson, D.R. (2002) Model selection and multimodel inference: a practical information-theoretic approach. Second edition. Springer. New York.
	}
\examples{
	#simulate 4 AICc values for 4 models.
	AICc_Array<-rnorm(4, mean=10,sd=1)
	#calculate the delta AICc
	Delta.AICc.Array<-AICc_Array-min(AICc_Array)
	#calculate the Akaike weight
	AkaikeWeight(Delta.AICc.Array)
	}
