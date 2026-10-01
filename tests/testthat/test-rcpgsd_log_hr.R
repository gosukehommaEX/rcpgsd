test_that("each estimator returns the corresponding estimate", {
  set.seed(31)
  for (rep_i in 1:10) {
    n   <- 60
    trt <- rep(0:1, each = n / 2)
    tt  <- rexp(n, 0.1 * 0.7^trt)
    cc  <- runif(n, 0, 20)
    rs  <- .rcpgsd_riskset(pmin(tt, cc), as.integer(tt <= cc), trt)
    expect_identical(.rcpgsd_log_hr(rs, "cox"), .rcpgsd_cox1(rs))
    expect_identical(.rcpgsd_log_hr(rs, "peto"),
                     .rcpgsd_logrank(rs, "peto")$log_hr)
    expect_identical(.rcpgsd_log_hr(rs, "pike"),
                     .rcpgsd_logrank(rs, "pike")$log_hr)
  }
})

test_that("Peto and Pike estimates follow their definitions", {
  # Events at time 1 (experimental; 2 experimental and 2 control at risk),
  # time 2 (control; 1 and 2 at risk) and time 3 (experimental; 1 and 1 at
  # risk); the control patient followed up to time 4 is censored
  rs <- .rcpgsd_riskset(time = c(1, 2, 3, 4), status = c(1, 1, 1, 0),
                        trt = c(1, 0, 1, 0))
  o_t <- 2
  e_t <- 2 / 4 + 1 / 3 + 1 / 2
  v   <- 2 * 2 / 16 + 1 * 2 / 9 + 1 * 1 / 4
  expect_equal(.rcpgsd_log_hr(rs, "peto"), (o_t - e_t) / v)
  expect_equal(.rcpgsd_log_hr(rs, "pike"),
               log(o_t / e_t) - log((3 - o_t) / (3 - e_t)))
})

test_that("infinite and missing estimates are handled consistently", {
  # All events in the control group
  rs_c <- .rcpgsd_riskset(time = c(1, 2, 3, 4), status = c(1, 1, 0, 0),
                          trt = c(0, 0, 1, 1))
  expect_identical(.rcpgsd_log_hr(rs_c, "cox"), -Inf)
  expect_identical(.rcpgsd_log_hr(rs_c, "pike"), -Inf)
  expect_true(is.finite(.rcpgsd_log_hr(rs_c, "peto")))
  expect_lt(.rcpgsd_log_hr(rs_c, "peto"), 0)
  # No event while both groups are at risk
  rs_n <- .rcpgsd_riskset(time = c(1, 2, 3), status = c(0, 1, 1),
                          trt = c(0, 1, 1))
  for (e in c("cox", "peto", "pike")) {
    expect_identical(.rcpgsd_log_hr(rs_n, e), NA_real_)
  }
  # No patients
  rs_0 <- .rcpgsd_riskset(numeric(0), integer(0), integer(0))
  for (e in c("cox", "peto", "pike")) {
    expect_identical(.rcpgsd_log_hr(rs_0, e), NA_real_)
  }
})
