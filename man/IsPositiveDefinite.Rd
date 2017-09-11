\name{IsPositiveDefinite}
\alias{IsPositiveDefinite}
\title{
	Check to make sure a variance covariance network (from a tree or network) is positive definite
	}

\description{
	Return TRUE if it is positive definite, FALSE otherwise
	}

\usage{
	IsPositiveDefinite(V.modified)
	}

\arguments{
  \item{V.modified}{
 a matrix
 }
}

\details{
	Calculates the eigen values from a matrix, sees if the minimum is greater than zero
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
