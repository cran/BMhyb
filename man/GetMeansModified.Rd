\name{GetMeansModified}
\alias{GetMeansModified}
\title{
 	The species means for the network model
	}

\description{
	This function returns the species mean for the network model. For the non hybrid species, it returns as a parameter \eqn{\mu} and for hybrid species it returns
  \eqn{\mu + \log \beta}.
	}

\usage{
	GetMeansModified(x, phy, flow, actual.params)
	}

\arguments{
  \item{x}{
   vector of the parameters: including the rate of evolution \eqn{\sigma^2}, overall mean \eqn{\mu}, measurement error SE, hybrid vigor \eqn{\beta}, variation at the burst of hybridization \eqn{v_H}.
 	}
 \item{phy}{
 	a class of phylo tree.
 	}
 \item{flow}{
 	a flow struture of gene flow.
 	}
 \item{actual.params}{
 	The parameters that are related to the hybridization including the hybrid vigor \eqn{\beta} and variation at the burst of hybridization \eqn{v_H}.
 	}
	}

\value{
	the species mean.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Study trait evolution on phylogenetic networks}. Submitted.
	}

\examples{
	#number of nonhybrid species
	ntax.nonhybrid<-5
	#number of hybrid species
	ntax.hybrid<-2
	#simulate network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
    #set the parameter values
    sigma.sq <- 0.01;mu <- 1;SE <- 0;bt<-12;vh<-0
    par<-c(sigma.sq,mu,SE,bt,vh)
    names(par)<-c("sigma.sq","mu","SE","bt","vh")
    actual.params<-12
    names(actual.params)<-"bt"
	#calculate the species means
	GetMeansModified(par, network$phy, network$flow,actual.params)
	}
