context('Running cv.glasso on heterogeneous dataset with fixed grouping')

# Generate test data (all components from same distribution)
n.comps = 4
p = 5
Mu = matrix(rep(0, p), p, n.comps)
Sigma = array(diag(p), c(p, p, n.comps))

test.data = simMIX(1000, n.comps, rep(0.25, n.comps), Mu, Sigma)

one.group.result = het.cv.glasso(test.data$X[test.data$S==1,])

test_that("No errors when running with one group.",
          expect_true(!is.null(one.group.result)))

full.result = het.cv.glasso(test.data$X, test.data$S)
 
test_that("Correct means", 
           expect_less_than(sum(abs(full.result$Mu-Mu)), 0.1*p*n.comps))
test_that("Correct variances", 
         expect_less_than(sum(
           abs(full.result$Sigma.diag-sapply(1:n.comps, 
                                             function(x) diag(Sigma[,,x])))), 0.1*p*n.comps))
 
test_that("Correct covariances", 
           expect_less_than(sum(abs(full.result$Sig-Sigma)), 0.1*p*p*n.comps))

context('Running mixglasso on heterogeneous dataset with unknown grouping')

# Generate heterogeneous dataset
Mu = sapply(1:n.comps, function(n.comp) rnorm(p,0,5))
test.data = sim.mix.networks(1000, p, n.comps, Mu=Mu)

# Try single n.comp
mixglasso.single = mixglasso(test.data$data, n.comps)

test_that("Rand index on inferred components with mixglasso is > 0.7", 
					expect_more_than(adjustedRandIndex(test.data$comp, mixglasso.single$comp), 0.7))

# Try multiple n.comp
mixglasso.mult = mixglasso(test.data$data, 2:6)

test_that("Correct BIC is selected", 
					expect_equal(mixglasso.mult$n.comp[mixglasso.mult$bic.opt], n.comps))