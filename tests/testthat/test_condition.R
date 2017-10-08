test_that("ConditionBadCichlid", {
  utils::data("cichlid")
  free.parameters<-rep(TRUE, 5)
  names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
  free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  x <- c(6, mean(cichlid$data)) #terrible values
  expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=x, data=cichlid$data, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], allow.extrapolation=FALSE, do.kappa.check=TRUE))
  #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
#  matrix.condition <- kappa(V.modified, exact=TRUE) high
})

test_that("ConditionBadNicotiana", {
  utils::data("nicotiana")
  free.parameters<-rep(TRUE, 5)
  names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
  free.parameters[which(names(free.parameters)=="bt")]<-FALSE
  free.parameters[which(names(free.parameters)=="vh")]<-FALSE
  x <- c(6, mean(nicotiana$data)) #terrible values
  expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=x, data=nicotiana$data, phy=nicotiana$phy, flow=nicotiana$flow, actual.params=free.parameters[which(free.parameters)], allow.extrapolation=FALSE, do.kappa.check=TRUE))
  #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
#  matrix.condition <- kappa(V.modified, exact=TRUE) high
})

#test_that("BasicRun",{
#  utils::data("cichlid")
#  result <- BMhyb(cichlid$data, cichlid$phy, cichlid$flow, n.points=100)
#  expect_equal(class(result), "data.frame")
#  expect_equal(nrow(result), 4)
#  expect_equal(ncol(result), 22)
#})
