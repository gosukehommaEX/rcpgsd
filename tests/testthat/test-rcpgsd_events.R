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

test_that("drop-out: expected times and proportions match independent values", {
  # Expected values computed in Python (scipy) by double numerical
  # integration over entry time and event time, with Brent's method
  ev <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                      lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                      n_accrual = nacc_E, lambda_d = -log(1 - 0.05) / 12)
  expect_equal(ev$analysis$time, c(15.3483597624, 24.3879408832,
                                   41.6548528189), tolerance = 1e-7)
  expect_equal(ev$fs[1, ], c(0.061374706, 0.0875828568, 0.0967244115),
               tolerance = 1e-6)
  expect_equal(ev$events_region[, 3], c(43.5259851591, 203.2370074204,
                                        203.2370074204), tolerance = 1e-7)
  expect_equal(colSums(ev$events_region), c(150, 300, 450),
               tolerance = 1e-8)
  # Group-specific drop-out (control, experimental)
  ev_g <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                        lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                        n_accrual = nacc_E, lambda_d = c(0.004, 0.006))
  expect_equal(ev_g$analysis$time, c(15.3722613466, 24.4891971162,
                                     42.3836851936), tolerance = 1e-7)
  expect_equal(ev_g$fs[1, ], c(0.061559049, 0.0877835066, 0.0969343705),
               tolerance = 1e-6)
  # Drop-out delays every analysis
  ev0 <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                       lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                       n_accrual = nacc_E)
  expect_true(all(ev$analysis$time > ev0$analysis$time))
})

test_that("drop-out: expected events agree with numerical integration", {
  # Two regions, one accrual interval (0, 10); drop-out hazards by region
  # (rows) for the control and experimental groups (columns)
  n_acc <- matrix(c(100, 100), ncol = 1)
  lam   <- c(0.1, 0.1)
  hr    <- c(0.5, 1)
  ld    <- rbind(c(0.02, 0.03), c(0, 0.05))
  by_int <- function(t, n, lambda, eta, a = 0, b = 10) {
    if (t <= a) return(0)
    p_event <- function(u) {
      vapply(u, function(v) {
        integrate(function(x) lambda * exp(-(lambda + eta) * x), 0, t - v,
                  rel.tol = 1e-10)$value
      }, numeric(1))
    }
    n / (b - a) * integrate(p_event, a, min(b, t), rel.tol = 1e-10)$value
  }
  for (t in c(20, 5)) {
    e <- .rcpgsd_expected_events(t, n_acc, c(0, 10), lam, hr, r = 1,
                                 lambda_d = ld)
    ref1 <- by_int(t, 50, 0.05, 0.03) + by_int(t, 50, 0.1, 0.02)
    ref2 <- by_int(t, 50, 0.1, 0.05) + by_int(t, 50, 0.1, 0)
    expect_equal(e$events, c(ref1, ref2), tolerance = 1e-8)
  }
})

test_that("drop-out hazards can be given in three equivalent forms", {
  args <- list(events = c(150, 300, 450), hr = 0.72, lambda_c = log(2) / 12,
               accrual_time = c(0, 6, 18), n_accrual = nacc_E)
  ev_1 <- do.call(rcpgsd_events, c(args, list(lambda_d = 0.01)))
  ev_2 <- do.call(rcpgsd_events, c(args, list(lambda_d = c(0.01, 0.01))))
  ev_m <- do.call(rcpgsd_events, c(args, list(lambda_d = matrix(0.01, 3, 2))))
  expect_identical(ev_1, ev_2)
  expect_identical(ev_1, ev_m)
  ev_g <- do.call(rcpgsd_events, c(args, list(lambda_d = c(0.004, 0.006))))
  expect_equal(ev_g$lambda_d, rbind(c(0.004, 0.006), c(0.004, 0.006),
                                    c(0.004, 0.006)))
  # Zero drop-out reproduces the default exactly
  ev_0 <- do.call(rcpgsd_events, c(args, list(lambda_d = 0)))
  ev_d <- do.call(rcpgsd_events, args)
  expect_identical(ev_0, ev_d)
  # Region-specific drop-out: more drop-out in region 1 lowers its share
  ld_r <- rbind(c(0.03, 0.03), c(0, 0), c(0, 0))
  ev_r <- do.call(rcpgsd_events, c(args, list(lambda_d = ld_r)))
  expect_true(all(ev_r$fs[1, ] < ev_d$fs[1, ]))
})

test_that("drop-out: expected proportions agree with the simulation", {
  # Expected values computed independently in Python
  ld  <- rbind(c(0.03, 0.02), c(0.005, 0.01))
  ev  <- rcpgsd_events(events = c(100, 200), hr = 0.7, lambda_c = 0.05,
                       accrual_time = c(0, 12), n_accrual = c(40, 260),
                       lambda_d = ld)
  expect_equal(ev$analysis$time, c(16.46124438, 39.56792814),
               tolerance = 1e-7)
  expect_equal(ev$fs[1, ], c(0.12314155, 0.11148181), tolerance = 1e-6)
  sim <- rcpgsd_sim(list(upper = c(2.963, 1.969)), events = c(100, 200),
                    hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
                    n_accrual = c(40, 260), lambda_d = ld, nsim = 150,
                    seed = 11)
  expect_lt(max(abs(ev$fs - sim$fs_mean)), 0.01)
  # Difference in mean analysis time within 4 Monte Carlo standard errors
  se_t <- sim$analysis$time_sd / sqrt(150)
  expect_true(all(abs(ev$analysis$time - sim$analysis$time_mean) <
                    4 * se_t))
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd_events(events = c(150, 600), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E), "smaller than")
  expect_error(rcpgsd_events(events = c(300, 150), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E), "increasing")
  expect_error(rcpgsd_events(events = c(150, 300), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E, lambda_d = -0.01),
               "lambda_d")
  expect_error(rcpgsd_events(events = c(150, 300), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E, lambda_d = c(0.01, 0.01,
                                                              0.01)),
               "lambda_d")
  expect_error(rcpgsd_events(events = c(150, 300), hr = 0.72,
                             lambda_c = 0.05, accrual_time = c(0, 6, 18),
                             n_accrual = nacc_E,
                             lambda_d = matrix(0.01, 2, 2)), "lambda_d")
  # Expected total number of events with drop-out is about 136.8
  expect_error(rcpgsd_events(events = c(100, 140), hr = 0.7,
                             lambda_c = 0.05, accrual_time = c(0, 12),
                             n_accrual = c(40, 260), lambda_d = 0.05),
               "expected total number of events")
})
