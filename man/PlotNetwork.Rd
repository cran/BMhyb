\name{PlotNetwork}
\alias{PlotNetwork}
\title{
	Phylogenetic Network Plot
	}

\description{
	Plots the phylogenetic network: including the phylogenetic tree and arrows that indicate the direction of gene flows.
	}

\usage{
	PlotNetwork(phy, flow, col.non="black", col.hybrid="red", col.donor="blue",
	name.padding=1.5, cex=1, xlab="", bty="n", head.length=0.2, edge.width=2,
	col.tree="darkgray", col.arrow="red", arrow.width=1, try.rotations=FALSE, ...)
 }

\arguments{
  \item{phy}{a class phylo tree object.}
  \item{flow}{a flow structure of gene flow.}
  \item{col.non}{color of non hybrid taxa.}
  \item{col.hybrid}{color of hybrid taxa.}
  \item{col.donor}{color of the donor.}
  \item{name.padding}{the size that allowed for the plot in the x axis.}
  \item{cex}{number indicating the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50 percent larger, 0.5 is 50 percent smaller, etc.}
  \item{xlab}{X axis label using font and character expansion par("font.lab") and color par("col.lab").}
  \item{bty}{the type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".}
  \item{head.length}{length of the head.}
  \item{edge.width}{width of the edges.}
  \item{col.tree}{color of the tree.}
  \item{col.arrow}{color of gene flow arrow.}
  \item{arrow.width}{width of the arrow.}
	\item{try.rotations}{Each time try randomly rotating nodes of the tree.}
  \item{...}{other arguments}
  }

\details{
	This function generates plots for network models. It first draws the phylogenetic tree for the given species. Then uses the gene flow data to draw the arrow from the donor species to recipient species. The hybrid taxa are shown in \eqn{col.hybrid} with the gene flow direction from the donor species (colored in \eqn{col.donor}).
  }

\value{
	plot of the network that involves the phylogeny and the gene flow direction.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#set the number of non hybrid species
	ntax.nonhybrid<-10
	#set the number of hybrid species
	ntax.hybrid<-3
	#simulate the network with desired species
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
	#plot the network
	PlotNetwork(network$phy,network$flow)
	}
