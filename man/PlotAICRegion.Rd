\name{PlotAICRegion}
\alias{PlotAICRegion}
\title{
	Plot the region from which points were estimated
	}

\description{
	Plots all points evaluated.
	}

\usage{
PlotAICRegion(sims, show.finite.only=TRUE, true.params=NULL,  ...)
 }

\arguments{
  \item{sims}{sims object from BMhybGrid.}
  \item{show.finite.only}{If true, only show where there were finite likelihoods; if false, show colors for these too}
  \item{true.params}{if known the true values used in a simulation}
  \item{...}{other arguments}
  }

\details{
	This plots the results of a grid search. By default, it shows the regions where there are finite values without distinguishing the best region (future work will add color and contour lines). A red diamond shows the best value found; a blue asterisk shows the true value (if known). Additional arguments to plot can be passed.
  }

\value{
	plot of the space.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Trait evolution on phylogenetic networks}. Submitted.
	}
