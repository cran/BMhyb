\name{GetClade}
\alias{GetClade}
\title{
	Get clade from the tree
	}

\description{
	Search and report the nodes of the tree that have the desired clade size.
	}

\usage{
	GetClade(phy, clade.size)
	}

\arguments{
  \item{phy}{an object of class 'phylo'.}
  \item{clade.size}{size of the clade.}
	}

\details{
	This function uses the 'phylo' class where the edges(ancestor-descendant relationship) of the tree are used to identify the interior node with desired number of descendent (clade.size) on the tips. It applies the function \code{findMRCA} in \emph{phytools} to search the ancestor.
	}

\value{
	the interior node that has descendants of size clade.size on the tip.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted.
	}

\examples{
	library(ape)
	library(phytools)
	#simulate a tree
 	phy<-rtree(3)
 	#plot the tree
 	plot(phy)
 	#set up the clade size
 	clade.size<-2
 	#search the nodes that has desired clade size
 	GetClade(phy,clade.size) #node 5 will be return
 	}
