\name{ConvertLog1P}
\alias{ConvertLog1P}
\title{
	Logarithms
      }

\description{
	ConvertLog1P computes log(1+x) accurately also for |x| << 1.
	}

\usage{
	ConvertLog1P(x)
	}

\arguments{
  \item{x}{
         x: a numeric or complex vector
         }
         }

\details{
  see function log1p(x) in for more detail description.
	}

\value{
 A vector of the same length as 'x' containing the transformed value.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
\references{
  Becker, R.A., Chambers, J.M. and Wilks, A.R. 1988. \emph{The New S Language}.  Wadsworth and Brooks Cole.
	}
