# Partial log-likelihood written directly from the raw data (no risk set
# summary), used as an independent reference
pl_direct <- function(b, time, status, trt) {
  sum(vapply(which(status == 1L), function(i) {
    at_risk <- time >= time[i]
    b * trt[i] - log(sum(exp(b * trt[at_risk])))
  }, numeric(1L)))
}

test_that("estimates maximize the partial likelihood", {
  set.seed(21)
  for (rep_i in 1:20) {
    n   <- sample(10:60, 1L)
    trt <- rbinom(n, 1, 0.5)
    tt  <- rexp(n, 0.1 * 0.6^trt)
    cc  <- runif(n, 0, 25)
    time   <- pmin(tt, cc)
    status <- as.integer(tt <= cc)
    b_hat <- .rcpgsd_cox1(.rcpgsd_riskset(time, status, trt))
    if (!is.finite(b_hat)) next
    b_ref <- optimize(pl_direct, c(-10, 10), time = time, status = status,
                      trt = trt, maximum = TRUE, tol = 1e-10)$maximum
    expect_equal(b_hat, b_ref, tolerance = 1e-6)
  }
})

test_that("estimates match survival::coxph", {
  skip_if_not_installed("survival")
  set.seed(22)
  for (rep_i in 1:10) {
    n   <- 80
    trt <- rep(0:1, each = n / 2)
    tt  <- rexp(n, 0.08 * 0.7^trt)
    cc  <- runif(n, 0, 30)
    time   <- pmin(tt, cc)
    status <- as.integer(tt <= cc)
    b_hat <- .rcpgsd_cox1(.rcpgsd_riskset(time, status, trt))
    fit   <- survival::coxph(survival::Surv(time, status) ~ trt)
    expect_equal(b_hat, unname(stats::coef(fit)), tolerance = 1e-6)
  }
})

test_that("monotone and flat partial likelihoods are handled", {
  # All events while both groups are at risk occur in the control group
  expect_equal(.rcpgsd_cox1(.rcpgsd_riskset(c(1, 2, 3, 4),
                                            c(1L, 1L, 0L, 0L),
                                            c(0L, 0L, 1L, 1L))), -Inf)
  # All such events occur in the experimental group
  expect_equal(.rcpgsd_cox1(.rcpgsd_riskset(c(1, 2, 3, 4),
                                            c(1L, 1L, 0L, 0L),
                                            c(1L, 1L, 0L, 0L))), Inf)
  # No events
  expect_true(is.na(.rcpgsd_cox1(.rcpgsd_riskset(c(1, 2), c(0L, 0L),
                                                 c(0L, 1L)))))
  # A single group only
  expect_true(is.na(.rcpgsd_cox1(.rcpgsd_riskset(c(1, 2), c(1L, 1L),
                                                 c(1L, 1L)))))
})
