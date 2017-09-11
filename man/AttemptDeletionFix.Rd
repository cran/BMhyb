\name{AttemptDeletionFix}
\alias{AttemptDeletionFix}
\title{
 Delete taxa to improve matrix condition
 }
\description{
	This will delete a taxon or taxa not involved as a source of hybrid gene flow or recipient of gene flow to see if this fixes a network that is not positive definite. It's drastic action.
	}

\usage{
	AttemptDeletionFix(phy, flow, params=c(1,0,0.1, 0, 0), m.vector = c(1,2))
	}

\arguments{
  \item{phy}{
	 	an object of class 'phylo'.
 		}

 	\item{flow}{
 		the flow matrix
 		}

    \item{params}{
      parameters to use for creating the modified variance covariance matrix
      }

      \item{m.vector}{
        how many taxa to try deleting.
	}
  }

\details{
m.vector=c(1,2) means it will try deleting a single taxon at random; if it works, stop there, if not, replace that one and delete another, until all single taxon deletions have been tried. Then it'll repeat for all possible pairs of taxa. You can add on to this (all groups of three taxa, etc.) though the space of things to try gets very large.
}

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
