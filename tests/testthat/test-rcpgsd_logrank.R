test_that("log-rank quantities match a hand calculation", {
  rs <- .rcpgsd_riskset(time   = c(1, 2, 3, 4, 5),
                        status = c(1, 0, 1, 1, 0),
                        trt    = c(1, 0, 0, 1, 1))
  e_t <- 3 / 5 + 2 / 3 + 2 / 2
  v   <- 3 * 2 / 25 + 2 * 1 / 9
  lr_pike <- .rcpgsd_logrank(rs, "pike")
  lr_peto <- .rcpgsd_logrank(rs, "peto")
  expect_equal(lr_pike$d, 3L)
  expect_equal(lr_pike$o_t, 2)
  expect_equal(lr_pike$e_t, e_t)
  expect_equal(lr_pike$v, v)
  expect_equal(lr_pike$z, (e_t - 2) / sqrt(v))
  expect_equal(lr_pike$log_hr, log(2 / e_t) - log(1 / (3 - e_t)))
  expect_equal(lr_peto$log_hr, (2 - e_t) / v)
})

test_that("log-rank quantities match survival::survdiff", {
  skip_if_not_installed("survival")
  set.seed(11)
  n   <- 120
  trt <- rep(0:1, each = n / 2)
  tt  <- rexp(n, 0.1 * 0.7^trt)
  cc  <- runif(n, 0, 20)
  time   <- pmin(tt, cc)
  status <- as.integer(tt <= cc)
  lr <- .rcpgsd_logrank(.rcpgsd_riskset(time, status, trt))
  sd <- survival::survdiff(survival::Surv(time, status) ~ trt)
  expect_equal(lr$o_t, unname(sd$obs[2]))
  expect_equal(lr$e_t, unname(sd$exp[2]), tolerance = 1e-10)
  expect_equal(lr$v, sd$var[2, 2], tolerance = 1e-10)
  expect_equal(lr$z^2, sd$chisq, tolerance = 1e-10)
})

test_that("no events give missing statistics", {
  lr <- .rcpgsd_logrank(.rcpgsd_riskset(c(1, 2), c(0L, 0L), c(0L, 1L)))
  expect_equal(lr$d, 0L)
  expect_true(is.na(lr$z))
  expect_true(is.na(lr$log_hr))
})
