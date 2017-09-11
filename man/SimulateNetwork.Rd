\name{SimulateNetwork}
\alias{SimulateNetwork}
\title{
	Simulate phylogenetic network
	}

\description{
	Simulate network with desired number of taxa and number of hybrids.
	}

\usage{
	SimulateNetwork(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5,
	origin.type=c("clade", "individual"), birth = 1, death = 1, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
	}

\arguments{
  \item{ntax.nonhybrid}{
	number of non hybrid taxa.
	}

  \item{ntax.hybrid}{
	number of hybrid taxa.
	}

  \item{flow.proportion}{
  	the gene flow proportion from the parent.
  	}

  \item{origin.type}{
  	the type where the hybrids were formed. clade: hybrid formed and then speciate; individual: hybrid species formed and evolved without further speciation.
  	}

  \item{birth}{
  	birth rate for the tree speciation.
  	}

  \item{death}{
  	death rate (extinction) of the tree.
  	}

  \item{sample.f}{
  	sampling frequency.
  	}

  \item{tree.height}{
  	the height of the tree.
  	}

  \item{allow.ghost}{
  	allows ghost lineage that persists for a while, donates genes to a hybrid, and eventually goes extinct (or at least unsampled). Otherwise, hybridization event
  	will be between coeval species with extant descendants.
  	}
	}

\details{
	This function generates a tree with relevant hybridization information (the structure of gene flow information).	Gene flow cannot go back in time. In general, gene flow can go forward in time via ghost lineages. If no ghost lineages are allowed, then there must be temporal overlap between the donor and recipient lineages.
	}

\value{
	a list that contains
	 \item{phy}{a birth death tree with number of taxa of ntax.nonhybrid + ntax.hybrid}
 	\item{flow}{The flow structure of the gene flow}
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	library(TreeSim)
	#set up the number for non hybrid species
	non.hybrid <- 3
	#set up the number for hybrid species
	hybrid <- 1
	#set up the gene flow proportion
	flow<- 0.5
	#set up the hybridization type to original
	origins<- "individual"
	#start to simulate the network
	network<-SimulateNetwork(ntax.nonhybrid=non.hybrid, ntax.hybrid=hybrid,
	flow.proportion=flow, origin.type="individual", birth = 1, death = 0.5,
	sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	#print out the result
	(network)
	}
