\name{SimulateTipData}
\alias{SimulateTipData}
\title{
	Simulate data on a network
	}

\description{
	Simulate comparative data under the given a network and relevant parameters.
	}

\usage{
	SimulateTipData(phy, flow, params, suffix="_DUPLICATE")
	}

\arguments{
  \item{phy}{
 a	phylo object.
	}

  \item{flow}{
	data.frame of the flow parameters.
	}

  \item{params}{
  	named vector of BMhyb parameter values.
  	}

  \item{suffix}{
	taxa are duplicated on the network; this suffix distinguishes them.
	 }
}
 \details{
	This function simulates tip data on a network, under a BMhyb model.
	}

\value{
	a named vector of tip values.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted.
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
	parameters <- c(0.01,1,1,0,0)
	names(parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	tip.data <- SimulateTipData(network$phy, network$flow, parameters)
	print(tip.data)
	}
