context('Simulating data from Gaussian mixture model')

# Generate test data (all components from same distribution)
n.comps = 4
p = 5 # dimensionality
Mu = matrix(rep(0, p), p, n.comps)
Sigma = array(diag(p), c(p, p, n.comps))
mix.prob = rep(0.25, n.comps)

test_that("No errors when using one Mu.",
          expect_true({simMIX(100, n.comps, mix.prob, Mu[,1], Sigma); TRUE}))

test_that("No errors when using one Sigma.",
          expect_true({simMIX(100, n.comps, mix.prob, Mu, Sigma[,,1]); TRUE}))

test_that("Errors when using wrong number of Sigma or Mu.", {
           expect_error(simMIX(100, n.comps, mix.prob, Mu, Sigma[,,1:3]))
           expect_error(simMIX(100, n.comps, mix.prob, Mu[,1:3], Sigma))}) 

test_that("Errors with misspecified mix.prob", {
          expect_error(simMIX(100, n.comps, rep(0.25, 3), Mu, Sigma))
          expect_error(simMIX(100, n.comps, c(mix.prob[1:3], 1.25) , Mu, Sigma))
          expect_error(simMIX(100, n.comps, c(mix.prob[1:3], -0.25) , Mu, Sigma))
          expect_error(simMIX(100, n.comps, c(mix.prob[1:3], 0.26) , Mu, Sigma))})

test.data = simMIX(1e4, n.comps, mix.prob, Mu, Sigma)
emp.mean = sapply(1:4, function(x) colMeans(test.data$X[test.data$S==x,]))
emp.cov = sapply(1:4, function(x) var(test.data$X[test.data$S==x,]), simplify='array')

test_that("Correct Mean and Sigma", {
          expect_less_than(sum(abs(emp.mean-Mu)), 0.1*p*n.comps)
          expect_less_than(sum(abs(emp.cov-Sigma)), 0.1*p*n.comps*n.comps)})

emp.mix = sapply(1:4, function(x) sum(test.data$S==x)/length(test.data$S))

test_that("Correct mixture components",
          expect_less_than(sum(abs(emp.mix-mix.prob)), 0.1))

test_that("t-distribution throws no error",
          expect_true({simMIX(100, n.comps, mix.prob, Mu[,1], Sigma, 
                                      dist='t', df=2);TRUE}))

context('Simulating means, covariances, and simulating data from Gaussian mixture model')

test.covariance = generateInvCov(10, 0.8)

test_that("Inverting inverse covariance",
          expect_true({solve(test.covariance); TRUE}))

test_that("Sparsity of inverse covariance",
          expect_equal((sum(test.covariance != 0) - 10)/(10*9), 0.2))

test.data = sim_mix_networks(1e4, p, n.comps)
emp.mean = sapply(1:4, function(x) colMeans(test.data$data[test.data$comp==x,]))
emp.cov = sapply(1:4, function(x) var(test.data$data[test.data$comp==x,]), simplify='array')

test_that("Correct Mean and Sigma", {
  expect_less_than(sum(abs(emp.mean-test.data$Mu)), 0.1*p*n.comps)
  expect_less_than(sum(abs(emp.cov-test.data$Sig)), 0.1*p*n.comps*n.comps)})
