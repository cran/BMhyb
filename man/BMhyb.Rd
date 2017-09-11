\name{BMhyb}
\alias{BMhyb}
\title{
 	Comparative method for studying hybridization using Brownian motion for trait evolution
	}

\description{
	This function fits the Brownian motion model of continuous character to investigate hybrid species through the hybrid vigor \eqn{\beta}, and the variation at the burst of hybridization \eqn{v_H}. Measurement error SE is also considered as well as the  parameters including the over all mean \eqn{\mu} and the overall variance \eqn{\sigma^2} in the typical Brownian motion model.
	}

\usage{
	BMhyb(data, phy, flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=TRUE,
	      get.se=TRUE, plot.se=TRUE, store.sims=FALSE, precision=2, auto.adjust=FALSE,
	      likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000,
              measurement.error = 0, do.kappa.check=FALSE, number.of.proportions = 101,
  	      number.of.proportions.adaptive = 101, allow.restart=TRUE,
  	      lower.bounds = c(0, -Inf, 0.000001, 0, 0),
 	      upper.bounds=c(10,Inf,100,100,100),
   	      badval.if.not.positive.definite=TRUE, attempt.deletion.fix=FALSE,
	      starting.values=NULL, n.random.start.points=5000,
	      do.Brissette.correction=FALSE,
	      do.Higham.correction=TRUE, do.DE.correction=FALSE)
	      }

\arguments{
  \item{data}{continuous trait data containing species information in vector format.}
  \item{phy}{a tree in phylo class.}
  \item{flow}{a structure of gene flow.}
  \item{opt.method}{the method for used for optimization. The default is Nelder-Mead.}
  \item{models}{the model used for analysis (see details).}
  \item{verbose}{a TRUE/FALSE argument to start optimization.}
  \item{get.se}{a TRUE/FALSE argument estimation for doing simulation to estimate parameter uncertainty (see details).}
  \item{plot.se}{a TRUE/FALSE argument for output the uncertainty plot for the model (see details).}
  \item{store.sims}{a TRUE/FALSE argument to record the the parameter estimates and relevant values.}
  \item{precision}{a numeric value to present the cutoff at which the user thinks the estimates become unreliable due to ill conditioned matrix.}
  \item{auto.adjust}{a TRUE/FALSE argument to adjust the the phylogeny.}
  \item{likelihood.precision}{a numerical value used for verifying the convergence of the estimation.}
  \item{allow.extrapolation}{a TRUE/FALSE argument. If TRUE, and the VCV matrix was ill-conditioned, would use splines to estimate its likelihood.}
  \item{n.points}{How many simulations to do to estimate confidence interval.}
  \item{measurement.error}{Estimate or fixed (see details).}
  \item{do.kappa.check}{Check for VCV matrix condition, approximate if poor. Leave FALSE unless you're cleverer than us.}
  \item{number.of.proportions}{Control how many points to use for extrapolation for ill-conditioned matrices.}
  \item{number.of.proportions.adaptive}{Control how many points to use for extrapolation for ill-conditioned matrices when calculating confidence (you may want to do fewer to speed things up).}
  \item{allow.restart}{if TRUE, stop and go back to calling function if this finds a better value.}
  \item{lower.bounds}{minimum values of all possible parameters (non-free ones will be deleted).}
	\item{upper.bounds}{maximum values of all possible parameters (non-free ones will be deleted).}
  \item{badval.if.not.positive.definite}{check to see if the network is numerically well-behaved).}
  \item{attempt.deletion.fix}{if the network is not well-behaved, try deleting taxa until it is.}
  \item{starting.values}{user given starting values (perhaps from BMhybGrid)}
  \item{n.random.start.points}{number of points to use for a starting grid search}
  \item{do.Brissette.correction}{If VCV is not positive definite, do Brissette et al. (2007) correction.}
  \item{do.Higham.correction}{If VCV is not positive definite, do Higham (2002).}
  \item{do.DE.correction}{If VCV is not positive definite, find a nearby matrix that is.}
  }

\details{Function BMhyb fits likelihood models for continuous characters. As an input, \pkg{BMhyb} requires a phylogenetic tree of the phylo class, a structure of gene
flow and a comparative data. Currently the 	method is developed for univariate analysis where the comparative data includes a single
trait for analyses. The full likelihood model includes several parameters: the ancestral state \eqn{\mu}, the overall rate of evolution \eqn{\sigma^2}, the hybrid vigor \eqn{\beta}, the hybrid burst variation at formation \eqn{v_H} and the measurement error SE.
The struture of gene flow is a five column table where the first and the second column contain the donor and recipient information. The third column is the information about the heritability factor \eqn{gamma} which is a fraction of the recipient trait that comes from the source.
For example, if one thought that 10% of a species' genome came from introgression from a particular species, \eqn{m} would be 0.1. The fourth column, recipient time, records time from the root of the recipient that counting forward from the root when the gene flow happened from the donor.
The fifth column, recipient time, records time from the root of the recipient that counting forward from the root when the gene flow happened from the donor. For more details, please see our paper. Note that gene flow out of and into a species occurs at a single time in this model, but it doesn't have to be at the same time for donor and recipient: for example, the genes could flow into an unsampled species, then later get passed to a recipient.

 The function allows some fixed values of parameters and treats others as free parameters: model 1 fixes \eqn{\beta}
	at 1 but 	allows \eqn{v_H} to vary; model 2 allows \eqn{\beta} to vary but fixes \eqn{v_H} at 0; model 3 fixes \eqn{\beta} at 1
	and \eqn{v_H} at 0; and model 4 allows both to vary.  BMhyb fits the model through maximum likelihood. When setting \eqn{get.se} to TRUE, the BMhyb will estimate
  parameter uncertainty through a sampling procedure (see help for \eqn{AdaptiveConfidenceIntervalSampling}. If \eqn{plot.se} is TRUE, a PDF will be generated
  showing these regions. Model averaged parameter estimates are calculated by the Akaike weight.

  \eqn{measurement.error} tells the program how to deal with measurement error. If set to NULL, it estimates a fixed one shared by all species.
  If you enter a single value (which could be zero if you measure your species really, really well), it uses this as fixed measurement error
  for all species. You can also enter a vector of measurement errors to allow species to have different ones; order should match that of taxa
  in vcv.phylo(phy). By default, it assumes zero measurement error.

  Parameter units matter; by convention, they are not reported in comparative methods, but that's unfortunate. If your analysis is of species
  body length, let's say you measure it in meters. You should probably log transform the measurements to fit the assumptions of Brownian
  motion (equal chance of a change by 1% per MY, regardless of scale, rather than equal chance of gaining a meter in length, regardless of
  whether the initial size is 10 m or 0.001 m; values can go arbitrarily bigger or smaller, rather than being bounded by zero).

  \eqn{data}: you should give this to the function in log(meters). I.e., if you measure sharks of length 1 m, 1.2 m, and 3 m, you
  should pass in 0.00, 0.18, and 1.10 as the traits.
  \eqn{measurement.error}: if you're not having the program calculate it, this should also be in units of log(meters). If you can
  measure shark length +/- 0.01 m, then this would be -4.61.

  The program returns variables that also have units. Assuming your tree has branch lengths in millions of years:
  \eqn{sigma-squared}, the rate of evolution, is in (log(meters))^2 / MY
  \eqn{mu}, the ancestral state, is in log(meters)
  \eqn{beta}, the scalar, is unitless. A value of 1.3 suggests that species formed from hybridization are 1.3 times larger than
  parent species. Note this is true for traits in log space (so in this case, for example, we'd expect a hybrid coming from
  species 2 m long to be 2.6 m long), but if you're using non-log transformed data the meaning of beta is different.
  \eqn{vH} is the variance that comes from a hybridization event. In this case, (log(meters))^2.
  \eqn{SE} is the inferred measurement error. It's in log(meters) units.

  If you wanted to compare things on the same scale, take into account the units. For example, you could compare the
  amount of variance for a species generated by one hybridization event:
  \eqn{sigma-squared} * tree height = variance coming from Brownian motion on the tree
  \eqn{v_H} = variance coming from the hybridization event
  \eqn{SE^2} = variance coming from measurement error.

  Though the model seems straightforward (Brownian motion on a basic network) there are numerical issues that can lead to difficulty calculating the likelihood (look up ill-conditioned matrices if you're curious).
  We have spent a lot of development time trying to address this issue and have various approximations in place to deal with this, including scaling the tree and trying transformations of the tree and then interpolation to get an estimate of the likelihood.
  There is a tradeoff. A cleaner approach would be to just return NA if the likelihood could not be calculated with high precision, and you can set the program to do this with some of the settings here (\eqn{allow.extrapolation}=FALSE, \eqn{auto.adjust}=FALSE, decreasing \eqn{precision} (a threshold for the condition number (kappa) of the matrix)).
  However, the real risk is that your hard-won dataset just generates an NA (or the equivalent, a very bad fixed negative log likelihood value, set to be something like 1e307 on many computers). Our defaults let you get an approximate answer; pay heed to the confidence intervals, as well. However, we recommend a "sniff test" for your results: the negative log likelihood values for your various models will vary, but should be something like 68.2 for one model, 74.3 for another, etc. If you get one model that is orders of magnitude different in log likelihood, something is likely deeply wrong.
 	}

\value{
	A summarized table including the type of model, the corresponding number of parameter, the parameter estimates, the likelihood values, the upper bound and lower bound of 	the parameters, and the Akaike weights for model averaging.
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2016. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted.

	Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. 	Springer, New York.
	}

\examples{
	#set up the number of non hybrid
	ntax.nonhybrid<-2
	#set up the number of non hybrid
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid,
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
    #simulate the tips data
	tips<-rnorm(ntax.nonhybrid+ntax.hybrid)
	names(tips)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")
	#run the analysis uses model 4 under optimization method (more than 5 seconds)
	\donttest{
	BMhyb(tips,network$phy,network$flow, opt.method="Nelder-Mead", models=4, verbose=TRUE,
	get.se=FALSE, plot.se=FALSE, store.sims=FALSE, precision=2, auto.adjust=FALSE,
	likelihood.precision=0.001, allow.extrapolation=TRUE)
	}
  #run the analysis uses model 4 under Grid method
  result<-BMhybGrid(tips,network$phy,network$flow, models=4, verbose=FALSE,
          get.se=FALSE, plot.se=FALSE, store.sims=FALSE, precision=2,
          auto.adjust=FALSE,likelihood.precision=0.001,
          allow.extrapolation=FALSE, n.points=10,
          measurement.error = 0, do.kappa.check=FALSE,
          number.of.proportions = 101, number.of.proportions.adaptive = 101,
          allow.restart=FALSE, lower.bounds = c(0, -Inf, 0.000001, 0, 0),
          upper.bounds=c(10,Inf,100,100,100),
          badval.if.not.positive.definite=TRUE,
          attempt.deletion.fix=FALSE, starting.values=NULL,
          do.Brissette.correction=FALSE, do.Higham.correction=TRUE,
          do.DE.correction=FALSE)
  #>result
  #  Model   sigma.sq         mu       bt       vh SE      AICc  NegLogL K
  #     4 0.02950803 0.02499673 8.064703 1.127043  0 -6.570164 2.714918 4
  #  sigma.sq.lower sigma.sq.upper mu.lower mu.upper bt.lower bt.upper vh.lower
  #             NA             NA       NA       NA       NA       NA       NA
  #  vh.upper SE.lower SE.upper penalty deltaAICc AkaikeWeight
  #       NA       NA       NA       0         0            1
	}
