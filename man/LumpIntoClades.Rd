\name{LumpIntoClades}
\alias{LumpIntoClades}
\title{
Lump into clades
}
\description{
This function Lump the hybrid species into a clades.
}

\usage{
LumpIntoClades(phy,flow)
}
\arguments{
  \item{phy}{a phylogeny of phylo class.}
  \item{flow}{structure of gene flow.}
   }

\details{
 The function summarizes the flow structure into the clade where the donor clades contribute more than one recipient clades.
    }

\value{
	a flow table that contains the donor clade, recipient clade, proportion of the gene flow,time from root donor, and time from root recipient.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#set the number of non hybrid species
	ntax.nonhybrid<-5
	#set the number of hybrid species
	ntax.hybrid<-2
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
	#show the flow
	network$flow
	#show the new flow table
	LumpIntoClades(network$phy, network$flow)
	}
