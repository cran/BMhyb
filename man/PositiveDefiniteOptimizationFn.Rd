\name{PositiveDefiniteOptimizationFn}
\alias{PositiveDefiniteOptimizationFn}
\title{
  Distance function for optimization
	}

\description{
Distance Function between two matrices for optimization searching of positive definite matrix.
	}

\usage{
	PositiveDefiniteOptimizationFn(x, original)
	}

\arguments{
  \item{x}{
 x: a vector that contains elements in the upper triangular part of the adjusted variance covariance matrix.
 }
 \item{original}{
  original:  a vector that contains elements in the upper triangular part of the original variance covariance matrix.
 }
}

\details{
It returns a distance between the new matrix and the original matrix.
If the new matrix contains negative element or is not positive definite, then it will return more penalty for the distance.
	}

\value{
a numerica value of distance
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
