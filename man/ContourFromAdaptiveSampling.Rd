\name{ContourFromAdaptiveSampling}
\alias{ContourFromAdaptiveSampling}
\title{
	Contour maps
	}

\description{
	This function generates the contour maps for a pair of parameters of interest.
	}

\usage{
	ContourFromAdaptiveSampling(sims, params.of.interest=NULL)
    }

\arguments{
  \item{sims}{a data frame with the simulation results.}
  \item{params.of.interest}{name of parameter that is of interest to be plotted.}
   }

\details{
	This function generates the contour plots for the parameter of interest. It uses the simulated data from generating the confidence interval under adaptive sampling method where the likelihood value and parameters values are stored. Then it uses the function \code{PlotConvexHull} to generate the contour maps.
	}

\value{
	It returns the plot of the contour map for a pair of parameters.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#simulate a data frame for two parameters
	sims<-matrix(rnorm(300),ncol=3)
	sims<-as.data.frame(sims)
	colnames(sims)<-c("neglnL","param1","param2")
	#plot the contour map for the parameters
	ContourFromAdaptiveSampling(sims, params.of.interest=NULL)
	}
