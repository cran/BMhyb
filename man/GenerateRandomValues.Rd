\name{GenerateRandomValues}
\alias{GenerateRandomValues}
\title{
  Generate random values
	}

\description{
 Generate random value for parameters
	}

\usage{
	GenerateRandomValues(data, free.parameters, lower, upper)
	}

\arguments{
  \item{data}{
  numerical vector that contains trait data
  	}
  \item{free.parameters}{
		TRUE/FALSE argument for free parameters
		}
	\item{lower}{
		lower bound for the new values
		}
  \item{upper}{
	  upper bound for the new values
		}
}

\details{
The function generates random values for parameter estimate. The bound for each paramter is set up using exponential and uniform distribution accordingly.
	}

\value{
 A vector contains the generate values for parameters.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
