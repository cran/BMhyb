\name{GetAncestor}
\alias{GetAncestor}
\title{
	Identify the ancestor
	}

\description{
	Given the descendant node for the tree, the function will return its ancestor node.
	}

\usage{
	GetAncestor(phy, node)
	}

\arguments{
  \item{phy}{an object of class 'phylo'.}
  \item{node}{descedant node.}
	}

\value{
	the ancestor node.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	library(ape)
	#simulate a tree of 3 taxa
 	phy<-rtree(3)
 	#plot the tree
 	plot(phy)
 	#descedant node
 	node<-1
 	#get the ancestor node, it will return 5.
 	GetAncestor(phy,node)
 	}
