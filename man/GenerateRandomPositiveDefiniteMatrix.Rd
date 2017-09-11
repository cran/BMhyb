\name{GenerateRandomPositiveDefiniteMatrix}
\alias{GenerateRandomPositiveDefiniteMatrix}
\title{
	Generate positive definite matrix
	}

\description{
	Generate postive definite matrix using QR decomposition
	}

\usage{
	GenerateRandomPositiveDefiniteMatrix(n, ev = runif(n, 0, 10))
	}

\arguments{
  \item{n}{
   size of matrix
   }
	\item{ev}{
	positive number
	}
}

\details{
  Use QR decomposition to generate positive definite matrix. See reference.

	}

\value{
 A positive definite matrix
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
\references{
  Ravi Varadhan, http://r.789695.n4.nabble.com/how-to-randomly-generate-a-n-by-n-positive-definite-matrix-in-R-td846858.html
	}
