\name{GetVModified}
\alias{GetVModified}
\title{
	 Variance covariance matrix for the network model
	}

\description{
	A function returns the variance covariance function for the network model under the Brownian motion for trait evolution.
	}

\usage{
	GetVModified(x, phy, flow, actual.params, measurement.error = NULL)
	}

\arguments{
  \item{x}{
	  vector of the parameters: including rate of evolution \eqn{\sigma^2}, overall mean \eqn{\mu}, measurement error SE, hybrid vigor \eqn{\beta}, and variation at the burst of hybridization \eqn{v_H}.
 	}
 \item{phy}{
 	a class of phylo tree.
 	}
 \item{flow}{
 	a struture of gene flow.
 	}
 \item{actual.params}{
 	The parameters that are related to the hybridization including the hybrid vigor \eqn{\beta} and variation at the burst of hybridization \eqn{v_H}.
 	}
	\item{measurement.error}{estimate or fixed (see details).}
}

\details{
	This function calculates the variance covariance matrix under the network model. The function first read tree in \emph{phylo} class and convert the variance covariance matrix use the Brownian motion model. It then uses the structure of the gene flow to identify the hybrid species, parents and the direction of the flow. The covariance of between the hybrid and non hybrid species is the weighted sum of the covariance from evolution along the tree plus evolution along the migration path. For more detail about the model description, please see Jhwueng and O'Meara 2016.
	}

\value{
 	the variance covariance matrix of size equal to the number of species on the tip of the phylogeny.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Study trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#set the number of non hybrid species
	ntax.nonhybrid<-4
	#set the number of hybrid species
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
	#assign the parameter values
	sigma.sq <- 0.01;mu <- 1;SE <- 0
	#get the variance covariance matrix
	GetVModified(c(sigma.sq,mu,SE), network$phy, network$flow, "bt")
}
