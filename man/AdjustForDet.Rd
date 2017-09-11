\name{AdjustForDet}
\alias{AdjustForDet}
\title{
 Adjust the phylogeny.
 }
\description{
	Utility function to adjust the phylogeny to deal with numerical issues.
	}

\usage{
	AdjustForDet(phy, max.attempts=100)
	}

\arguments{
  \item{phy}{
	 	an object of class 'phylo'.
 		}

 	\item{max.attempts}{
 		the maximum number of attempts to adjust the phylogeny. The default is 100.
 		}
	}

\details{
	This function first calculates the determinant of the phylogeny using the function \code{DetPass}. When the determinant of the variance covariance matrix for the phylogeny is negative, this function will adjust the phylogeny by slightly lengthening the edge lengths of the tree. This procedure repeats until a well adjusted phylogeny is returned.}

\value{
 an object of class 'phylo'.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted.
	}

\examples{
	library(ape)
	#simulate a tree of 5 taxa
	phy<-rtree(5)
	#check and adjust the phylogeny
	AdjustForDet(phy,max.attempts=100)
	}
