\name{AlterMatrixUsingDE}
\alias{AlterMatrixUsingDE}
\title{
	Transforming current variance covariance matrix to a positive definite matrix
	}

\description{
	Use differntial evolution method to search for a positive definite matrix.
	}

\usage{
	AlterMatrixUsingDE(V.modified)
	}

\arguments{
  \item{V.modified}{
 a matrix
 }
}

\details{
  The function is developed based on the theoretical method in Misha (2007), we apply it to search the positive definite matrix.
	It first starts search with the element in upper triangular matrix of the variance covariance matrix, then tries to search the
	positive definite matrix with all nonnegative elements.
	}

\value{
 a postive definite matrix.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
\references{
   Mishra, S.K. 2007. \emph{The nearest correlation matrix problem: Solution by differential evolution method of global optimization}. dx.doi.org/10.2139/ssrn.980403.

		}
