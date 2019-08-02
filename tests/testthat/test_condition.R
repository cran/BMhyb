# test_that("Can convert from old style to new", {
#   utils::data("nicotiana")
#   p <- BMhyb:::ConvertPhyAndFlowToPhygraph(nicotiana$phy, nicotiana$flow)
#   plot(p)
# })

test_that("Basic run can go", {
  utils::data("cichlid")
  skip_on_travis()
  result <- BMhyb(phy.graph=cichlid$phy.graph, traits=cichlid$trait, confidence.points=100, max.steps=2, control=list(reltol=1e-1, trace=6))
  expect_equal(class(result), "BMhybResult")
  expect_gte(result$best$NegLogLik, 10)
})

test_that("Likelihood of BM works", {
  utils::data("nicotiana")
  nicotiana$phy.graph$edge.length[which(nicotiana$phy.graph$edge.length<0)] <- 0
  phy <- nicotiana$phy.graph
  class(phy) <- "phylo"
  phy$reticulation <- NULL
  p2 <- nicotiana$phy.graph
  p2$reticulation <- p2$reticulation[-sequence(nrow(p2$reticulation)),]
  result_geiger <- geiger::fitContinuous(phy, nicotiana$trait, model="BM", bounds=c(1,1), ncores=1)
  result_BMhyb <- BMhyb:::TryComputeLikelihood(c(sigma.sq=result_geiger$opt$sigsq, mu=result_geiger$opt$z0), phy.graph=p2, traits=nicotiana$trait)
  expect_equal(-result_geiger$opt$lnL, result_BMhyb, tolerance = .01)
})

test_that("Simulation works", {
  network <- SimulateNetwork(ntax=5, nhybridizations=2)
  tips <- SimulateTips(network, mu=1.1, bt=3, vh=1.1, SE=1)
  result <- BMhyb(phy.graph=network, traits=tips, confidence.points=100, max.steps=2,  control=list(reltol=1e-1, trace=6))
  expect_equal(class(result), "BMhybResult")
  expect_gte(result$best$NegLogLik, 2)
  expect_equal(nrow(result$bad.region) + nrow(result$good.region), 101) #100 sim points, 1 starting point
})


test_that("VCV well formed", {
  utils::data("cichlid")
  VCV <- ComputeVCV(cichlid$phy.graph)
  expect_gte(kappa(VCV), 1)
})

test_that("Compute means works", {
  utils::data("cichlid")
  means <- ComputeMeans(cichlid$phy.graph, bt=2)
  expect_equal(unname(means['Neolamprologus_brevis']),0)
  expect_equal(unname(means['Neolamprologus_fasciatus']),log(2))
})

test_that("Issue 13 resolved", {
  create_paper_network <- function(gamma, t1, t2, t3){
    phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
    network <- list(phy = phy,
                    flow = data.frame(donor = "X",
                                      recipient = "R",
                                      gamma = gamma,
                                      time.from.root.donor = t1,
                                      time.from.root.recipient = t1 + t2))
    network$flow$donor <- as.character(network$flow$donor)
    network$flow$recipient <- as.character(network$flow$recipient)
    return(network)
  }
  gamma <- 0.5
  t1 <- 0.3; t2 <- 0.4; t3 <- 0.3; # unit height
  network <- create_paper_network(gamma, t1, t2, t3)
  phy.graph<-ConvertPhyAndFlowToPhygraph(phy=network$phy,flow=network$flow)
  VCV <- ComputeVCV(phy.graph)
  expect_equal(VCV['X', 'R'], 0.15)
  expect_equal(VCV['Y', 'R'], 0.35)
  expect_equal(VCV['R', 'Y'], 0.35)
})

test_that("Issue 14 resolved", {
  t1 <- 0.3; t2 <- 0.4; t3 <- 0.3; gamma <- 0.5
  phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
  ## Network
  don_recp <- expand.grid(c("X"), c("Y", "R"))
  network <- list(phy = phy,
                  flow = data.frame(donor = don_recp[,1],
                                    recipient = don_recp[,2],
                                    gamma = rep(gamma, 2),
                                    time.from.root.donor = rep(t1, 2),
                                    time.from.root.recipient = rep(t1, 2)))
  network$flow$donor <- as.character(network$flow$donor)
  network$flow$recipient <- as.character(network$flow$recipient)
  phy.graph<-ConvertPhyAndFlowToPhygraph(phy=network$phy,flow=network$flow)
  VCV <- ComputeVCV(phy.graph)
  expect_equal(VCV['R', 'Y'], 0.55)
  expect_equal(VCV['R', 'R'], 0.85)
})

test_that("Issue 15 resolved", {
  gamma1 <- 0.5; gamma2 <- 0.5;
  t1 <- 0.2; t2 <- 0.2; t3 <- 0.2; t4 <- 0.2; t5 <- 0.2;
  phy <- ape::read.tree(text = paste0("(((R:",t4+t5,",Y:",t4+t5,"):",t3,",X:",t3+t4+t5,"):",t1+t2,",Z:",t1+t2+t3+t4+t5,");"))
  ## Network
  don_recp <- rbind(expand.grid(c("Z"), c("Y", "R", "X")),
                    expand.grid(c("X"), c("R")))
  network <- list(phy = phy,
                  flow = data.frame(donor = don_recp[,1],
                                    recipient = don_recp[,2],
                                    gamma = c(rep(gamma1, 3), gamma2),
                                    time.from.root.donor = c(rep(t1, 3), t1+t2+t3+t4),
                                    time.from.root.recipient = c(rep(t1, 3), t1+t2+t3+t4)))
  network$flow$donor <- as.character(network$flow$donor)
  network$flow$recipient <- as.character(network$flow$recipient)
  phy.graph <- ConvertPhyAndFlowToPhygraph(phy=network$phy,flow=network$flow)
  VCV <- ComputeVCV(phy.graph)
  expect_equal(VCV['R', 'R'], 0.7)
})


# test_that("ConditionBadCichlid", {
#   utils::data("cichlid")
#   parameters.to.try <- c(sigma.sq=6, mu=mean(cichlid$data)) #terrible values
#   expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=parameters.to.try , data=cichlid$data, phy=cichlid$phy, flow=cichlid$flow, allow.extrapolation=FALSE, do.kappa.check=TRUE, measurement.error=0))
#   #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
# #  matrix.condition <- kappa(V.modified, exact=TRUE) high
# })
#
# test_that("ConditionBadNicotiana", {
#   utils::data("nicotiana")
#   parameters.to.try <- c(sigma.sq=6, mu=mean(nicotiana$data)) #terrible values
#   expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=parameters.to.try, data=nicotiana$data, phy=nicotiana$phy, flow=nicotiana$flow, allow.extrapolation=FALSE, do.kappa.check=TRUE, measurement.error=0))
#   #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
# #  matrix.condition <- kappa(V.modified, exact=TRUE) high
# })
#
# test_that("BasicRun",{
#   utils::data("cichlid")
#   result <- BMhyb(cichlid$data, cichlid$phy, cichlid$flow, n.points=100, get.se=TRUE, plot.se=FALSE, measurement.error=0,n.random.start.points=100)
#   expect_equal(class(result), "data.frame")
#   expect_equal(nrow(result), 4)
#   expect_equal(ncol(result), 22)
# })
#
# test_that("Issue 13 is solved", {
#   create_paper_network <- function(gamma, t1, t2, t3){
#       phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
#       network <- list(phy = phy,
#                       flow = data.frame(donor = "X",
#                                         recipient = "R",
#                                         gamma = gamma,
#                                         time.from.root.donor = t1,
#                                         time.from.root.recipient = t1 + t2))
#       network$flow$donor <- as.character(network$flow$donor)
#       network$flow$recipient <- as.character(network$flow$recipient)
#       return(network)
#   }
#
#   gamma <- 0.5
#   t1 <- 0.3; t2 <- 0.4; t3 <- 0.3; # unit height
#   network <- create_paper_network(gamma, t1, t2, t3)
#
#   #PlotNetwork(network$phy, network$flow)
#   #axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))
#
#   sigma2 = 1
#   x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
#   vcv_BMhyb <- GetVModified(x, network$phy, network$flow, measurement.error=0)
#   expect_equal(vcv_BMhyb, structure(c(0.65, 0.7, 0.15, 0.7, 1, 0, 0.15, 0, 1), .Dim = c(3L,
# + 3L), .Dimnames = list(c("R", "Y", "X"), c("R", "Y", "X")))) #from Issue 13
# })
#
# test_that("Issue 14 is solved", {
#   gamma <- 0.5
#   t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
#   phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
#   ## Network
#   don_recp <- expand.grid(c("X"), c("Y", "R"))
#   network <- list(phy = phy,
#                   flow = data.frame(donor = don_recp[,1],
#                                     recipient = don_recp[,2],
#                                     gamma = rep(gamma, 2),
#                                     time.from.root.donor = rep(t1, 2),
#                                     time.from.root.recipient = rep(t1, 2)))
#   network$flow$donor <- as.character(network$flow$donor)
#   network$flow$recipient <- as.character(network$flow$recipient)
#   ## Plot
#
#   sigma2 = 1
#   x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
#   expect_equal(GetVModified(x, network$phy, network$flow, measurement.error=0), structure(c(0.85, 0.55, 0.15, 0.55, 0.85, 0.15, 0.15, 0.15, 1
# ), .Dim = c(3L, 3L), .Dimnames = list(c("R", "Y", "X"), c("R",
# "Y", "X"))))
# })
#
# test_that("VH matters", {
#   gamma <- 0.5
#   t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
#   phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
#   ## Network
#   don_recp <- expand.grid(c("X"), c("Y", "R"))
#   network <- list(phy = phy,
#                   flow = data.frame(donor = don_recp[,1],
#                                     recipient = don_recp[,2],
#                                     gamma = rep(gamma, 2),
#                                     time.from.root.donor = rep(t1, 2),
#                                     time.from.root.recipient = rep(t1, 2)))
#   network$flow$donor <- as.character(network$flow$donor)
#   network$flow$recipient <- as.character(network$flow$recipient)
#   ## Plot
#
#   sigma2 = 1
#   x <- c(sigma.sq = sigma2, mu = 0, vh=0, SE=0)
#
#
#   V0 <- GetVModified(x, network$phy, network$flow, measurement.error=0)
#   vh.add = 3
#   y <- c(sigma.sq = sigma2, mu = 0, vh=vh.add,SE = 0)
#   V1 <- GetVModified(y, network$phy, network$flow, measurement.error=0)
#   expect_equal(V1[1,1]-vh.add, V0[1,1])
#   expect_equal(V1[3,3], V0[3,3])
# })
#
# #
# # test_that("MergingTrees", {
# #   # idea is a network, decomposed into trees -- does it come back
# #   phy1 <- ape::read.tree(text="(((A:2,B:2):1,C:3):5,D:8);")
# #   phy2 <- ape::read.tree(text="(((B:2,C:2):1,A:3):5,D:8);")
# #   phy.graph <- MergeTreesIntoPhyGraph(c(phy1, phy2))
# #   free.parameters<-rep(TRUE, 5)
# #   names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
# #   free.parameters[which(names(free.parameters)=="bt")]<-FALSE
# #   free.parameters[which(names(free.parameters)=="vh")]<-FALSE
# #   x <- c(1, 0, 1)
# #   results <- GetVandMFromIgraph(x, phy.graph, free.parameters)
# #   expect_equal(nrow(results$V.modified), 4)
# # })
# #
# # test_that("NonNetworkWorks", {
# #   # idea is a network, decomposed into trees -- does it come back
# #   phy1 <- ape::read.tree(text="(((A:2,B:2):1,C:3):5,D:8);")
# #   phy2 <- phy1
# #   phy.graph <- MergeTreesIntoPhyGraph(c(phy1, phy2))
# #   free.parameters<-rep(TRUE, 5)
# #   names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
# #   free.parameters[which(names(free.parameters)=="bt")]<-FALSE
# #   free.parameters[which(names(free.parameters)=="vh")]<-FALSE
# #   x <- c(1, 0, 1)
# #   results <- GetVandMFromIgraph(x, phy.graph, free.parameters)
# #   expect_lte(max(abs(results$V.modified - ape::vcv(phy1))),1e-8)
# # })
