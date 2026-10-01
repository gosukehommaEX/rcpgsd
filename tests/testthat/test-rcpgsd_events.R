# Scenario: 600 patients over 18 months, Japan (region 1) enrolling only
# after month 6, median survival 12 months in the control group. Expected
# values were computed independently in Python (scipy).
nacc_E <- rbind(c(0, 60), c(100, 170), c(100, 170))

test_that("expected times and proportions match independent values", {
  ev <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                      lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                      n_accrual = nacc_E)
  expect_s3_class(ev, "rcpgsd_events")
  expect_equal(ev$analysis$time, c(15.18413572, 23.75530021, 38.09367258),
               tolerance = 1e-7)
  expect_equal(ev$fs[1, ], c(0.06012201, 0.08629776, 0.09551110),
               tolerance = 1e-6)
  expect_equal(ev$events_region[, 3], c(42.97999592, 203.51000204,
                                        203.51000204), tolerance = 1e-7)
  expect_equal(colSums(ev$fs), rep(1, 3), tolerance = 1e-12)
  expect_equal(colSums(ev$events_region), c(150, 300, 450),
               tolerance = 1e-8)
  # Accrual ends at month 18: before that, 200 patients from months 0-6
  # plus 400 patients uniformly over months 6-18
  t1 <- ev$analysis$time[1]
  expect_equal(ev$analysis$enrolled, c(200 + 400 * (t1 - 6) / 12, 600, 600))
})

test_that("expected events agree with numerical integration", {
  # Two regions, one accrual interval (0, 10), evaluated at t = 20 and t = 5
  n_acc <- matrix(c(100, 100), ncol = 1)
  lam   <- c(0.1, 0.1)
  hr    <- c(0.5, 1)
  by_int <- function(t, n, lambda, a = 0, b = 10) {
    if (t <= a) return(0)
    n / (b - a) * integrate(function(u) 1 - exp(-lambda * (t - u)),
                            a, min(b, t))$value
  }
  for (t in c(20, 5)) {
    e <- .rcpgsd_expected_events(t, n_acc, c(0, 10), lam, hr, r = 1)
    ref1 <- by_int(t, 50, 0.05) + by_int(t, 50, 0.1)
    ref2 <- by_int(t, 100, 0.1)
    expect_equal(e$events, c(ref1, ref2), tolerance = 1e-8)
    expect_equal(e$enrolled, rep(100 * min(t, 10) / 10, 2))
  }
})

test_that("expected proportions agree with the patient-level simulation", {
  design <- list(upper = c(3.7103, 2.5114, 1.9930))
  ev  <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                       lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                       n_accrual = nacc_E)
  sim <- rcpgsd_sim(design, events = c(150, 300, 450), hr = 0.72,
                    lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                    n_accrual = nacc_E, nsim = 200, seed = 9)
  expect_lt(max(abs(ev$fs - sim$fs_mean)), 0.01)
  expect_lt(max(abs(ev$analysis$time - sim$analysis$time_mean)), 0.3)
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd_events(events = c(150, 600), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E), "smaller than")
  expect_error(rcpgsd_events(events = c(300, 150), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E), "increasing")
})
