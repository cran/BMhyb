\name{PlotConvexHull}
\alias{PlotConvexHull}
\title{
	Convex hull
	}

\description{
	Plotting the convex hall for given 2D datasets.
	}

\usage{
	PlotConvexHull(xcoord, ycoord, lcolor)
	}
\arguments{
  		\item{xcoord}{coordinate for the horizontal axis.}
  		\item{ycoord}{coordinate for the vertical axis.}
  		\item{lcolor}{color of the convect hull.}
  }

\details{
	This function uses the base function \emph{chull()}, which takes x,y coordinates of data and outputs a vector(hpts variable) of integers that indicate the points in dataset that comprise the convex hull around data. It then uses the base function lines to draw the convex hull using the x,y coordinates stored in hpts.
	}

\value{
	The convex hull plot for a set of 2D points.
 	}

\author{
	Ken Takagi
	}

\references{
	https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
	}

\examples{
	# Create a set of random data to plot convex hull around
	x<-rnorm(100,0.8,0.3)
	y<- rnorm(100,0.8,0.3)
	#get max and min of all x and y data for nice plotting
	xrange<-range(x)
	yrange<-range(y)
	#plot it up!
	plot(x,y,type="p",pch=1,col='black',xlim=c(xrange),ylim=c(yrange))
	PlotConvexHull(xcoord=x,ycoord=y,lcolor='red')
	}
