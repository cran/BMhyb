\name{DetPass}
\alias{DetPass}
\title{
	determinant of the matrix
	}

\description{
	Calculate the determinant of the matrix and check whether it is positive.
	}

\usage{
	DetPass(phy)
	}

\arguments{
  \item{phy}{
 an object of class 'phylo'.
 }
}

\details{
	This function first calculates variance covariance matrix C for the tree. It then calculates the determinant of the C, 0.0001*C, and 1000*C and verifies whether all the three determinants are positive and finite. If yes, a TRUE value is return which indicates that the C matrix is good for further use.
	}

\value{
 a TRUE/FALSE value.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	library(ape)
	##simulate a tree of 5 taxa
	phy<-rtree(5)
	#check whether the determiant is positive
	DetPass(phy)
	}
