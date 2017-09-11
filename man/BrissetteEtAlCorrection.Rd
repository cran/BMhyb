\name{BrissetteEtAlCorrection}
\alias{BrissetteEtAlCorrection}
\title{
	Transforming current variance covariance matrix to a positive definite matrix
	}

\description{
	Use  method to search for a positive definite matrix.
	}

\usage{
	BrissetteEtAlCorrection(V.modified, min.eigenvalue=1e-6, max.attempts=10)
	}

\arguments{
  \item{V.modified}{
   a matrix
   }
  \item{min.eigenvalue}{
   minimum of eigenvalues of V.modified
   }
  \item{max.attempts}{
   maimum number of attempt for transformation
  }
}

\details{
The function replaces negative eigenvalues of the variance covariance matrix to a small positive value,
and then transform it back to the matrix.
	}

\value{
 a postive definite matrix.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
\references{
  Brissette F., Khalili M. Leconte R. 2007. \emph{Efficient stochastic generation of mult-site synthetic precipitation data}. Journal of Hydrology 345:121-133.
	}
