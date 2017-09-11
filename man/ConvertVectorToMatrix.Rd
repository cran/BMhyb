\name{ConvertVectorToMatrix}
\alias{ConvertVectorToMatrix}
\title{
	Convert a vector into a matrix
	}

\description{
  The function outputs a matrix given a vector as its elements in the upper triangular part.
	}

\usage{
	ConvertVectorToMatrix(x)
	}

\arguments{
  \item{x}{
 x: a postive integer that indicates the number of elements in the upper triangular matrix.
 }
}

\details{
 The length of vector equal the number of elements in a squared upper triangular matrix.
 It assigns value in the vector to the matrix accordingly. Eventually, the output of the matrix
 will be the matrix with lower triangular part equal to the upper triangular part.
 }

\value{
A square matrix
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
