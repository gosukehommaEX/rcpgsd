# Overall plan: 600 patients at 50 patients per month over 12 months.
# Region 1 (Japan, 60 patients) may start late; regions 2 and 3 (270
# patients each) start at month 0. Expected values were computed
# independently in Python.
n_reg <- c(60, 270, 270)

test_that("a single delayed region is compensated by the other regions", {
  acc <- rcpgsd_accrual(n_reg, start = c(6, 0, 0), accrual_time = c(0, 12),
                        rate_overall = 50)
  expect_s3_class(acc, "rcpgsd_accrual")
  expect_equal(acc$accrual_time, c(0, 6, 12))
  expect_equal(acc$n_accrual, rbind(c(0, 60), c(150, 120), c(150, 120)))
  # Japan enrols uniformly from month 6 to month 12; the overall rate is kept
  expect_equal(acc$rate, rbind(c(0, 10), c(25, 20), c(25, 20)))
  expect_equal(colSums(acc$rate), c(50, 50))
  expect_equal(rowSums(acc$n_accrual), n_reg)
  # Without delay every region enrols uniformly over the whole period
  acc0 <- rcpgsd_accrual(n_reg, accrual_time = c(0, 12), rate_overall = 50)
  expect_equal(acc0$accrual_time, c(0, 12))
  expect_equal(acc0$n_accrual, matrix(n_reg, ncol = 1))
})

test_that("expected analysis times are those of the plan without delay", {
  acc  <- rcpgsd_accrual(n_reg, start = c(6, 0, 0), accrual_time = c(0, 12),
                         rate_overall = 50)
  acc0 <- rcpgsd_accrual(n_reg, accrual_time = c(0, 12), rate_overall = 50)
  args <- list(events = c(150, 300, 450), hr = 0.72, lambda_c = log(2) / 12)
  ev  <- do.call(rcpgsd_events, c(args, list(accrual_time = acc$accrual_time,
                                             n_accrual = acc$n_accrual)))
  ev0 <- do.call(rcpgsd_events, c(args, list(accrual_time = acc0$accrual_time,
                                             n_accrual = acc0$n_accrual)))
  expect_equal(ev$analysis$time, ev0$analysis$time, tolerance = 1e-8)
  expect_equal(ev$analysis$time, c(12.1162328769, 20.3833502791,
                                   34.7285491399), tolerance = 1e-7)
  # The delayed region has a smaller share of the events, mostly at the
  # interim analyses
  expect_equal(ev$fs[1, ], c(0.055963321, 0.0854786959, 0.0952486277),
               tolerance = 1e-6)
  expect_equal(ev0$fs[1, ], rep(0.1, 3), tolerance = 1e-8)
})

test_that("the feasibility condition is applied at its boundary", {
  # A single delay is feasible if start <= (1 - 60 / 600) * 12 = 10.8
  acc <- rcpgsd_accrual(n_reg, start = c(10.8, 0, 0),
                        accrual_time = c(0, 12), rate_overall = 50)
  expect_equal(acc$n_accrual, rbind(c(0, 60), c(270, 0), c(270, 0)))
  expect_equal(acc$rate[, 2], c(50, 0, 0))
  expect_error(rcpgsd_accrual(n_reg, start = c(11, 0, 0),
                              accrual_time = c(0, 12), rate_overall = 50),
               "cannot be compensated")
  expect_error(rcpgsd_accrual(n_reg, start = c(11, 0, 0),
                              accrual_time = c(0, 12), rate_overall = 50),
               "10.8")
})

test_that("several regions may start at different times", {
  # Region 2 (240 patients) starts at month 2, region 1 (60) at month 6
  acc <- rcpgsd_accrual(c(60, 240, 300), start = c(6, 2, 0),
                        accrual_time = c(0, 12), rate_overall = 50)
  expect_equal(acc$accrual_time, c(0, 2, 6, 12))
  expect_equal(acc$n_accrual, rbind(c(0, 0, 60), c(0, 96, 144),
                                    c(100, 104, 96)))
  expect_equal(acc$rate, rbind(c(0, 0, 10), c(0, 24, 24), c(50, 26, 16)))
  # 400 patients from month 4 (50 per month) plus 60 from month 6 exceed
  # the overall rate of 50 per month after month 6
  expect_error(rcpgsd_accrual(c(60, 400, 140), start = c(6, 4, 0),
                              accrual_time = c(0, 12), rate_overall = 50),
               "between times 6 and 12, region\\(s\\) 1, 2")
  expect_error(rcpgsd_accrual(n_reg, start = c(6, 2, 1),
                              accrual_time = c(0, 12), rate_overall = 50),
               "At least one region")
})

test_that("a piecewise overall accrual plan is kept", {
  # Ramp-up: 20 patients per month in months 0-3, 60 per month thereafter
  acc  <- rcpgsd_accrual(n_reg, start = c(6, 0, 0),
                         accrual_time = c(0, 3, 12), rate_overall = c(20, 60))
  acc0 <- rcpgsd_accrual(n_reg, accrual_time = c(0, 3, 12),
                         rate_overall = c(20, 60))
  expect_equal(acc$accrual_time, c(0, 3, 6, 12))
  expect_equal(acc$n_accrual, rbind(c(0, 0, 60), c(30, 90, 150),
                                    c(30, 90, 150)))
  expect_equal(acc$rate_overall, c(20, 60, 60))
  args <- list(events = c(150, 300, 450), hr = 0.72, lambda_c = log(2) / 12)
  ev  <- do.call(rcpgsd_events, c(args, list(accrual_time = acc$accrual_time,
                                             n_accrual = acc$n_accrual)))
  ev0 <- do.call(rcpgsd_events, c(args, list(accrual_time = acc0$accrual_time,
                                             n_accrual = acc0$n_accrual)))
  expect_equal(ev$analysis$time, ev0$analysis$time, tolerance = 1e-8)
  expect_equal(ev$analysis$time, c(12.9481500008, 21.2160034538,
                                   35.5624879505), tolerance = 1e-7)
  # Numbers and rates of the overall plan give the same result
  acc_n <- rcpgsd_accrual(n_reg, start = c(6, 0, 0),
                          accrual_time = c(0, 3, 12), n_overall = c(60, 540))
  expect_equal(acc_n$n_accrual, acc$n_accrual)
  expect_equal(acc_n$rate, acc$rate)
})

test_that("rounding keeps the regional totals", {
  acc <- rcpgsd_accrual(c(61, 269, 270), start = c(5, 0, 0),
                        accrual_time = c(0, 12), rate_overall = 50)
  expect_equal(acc$n_accrual_exact[2, ], c(124.7680890538, 144.2319109462),
               tolerance = 1e-9)
  expect_equal(acc$n_accrual, rbind(c(0, 61), c(125, 144), c(125, 145)))
  expect_equal(rowSums(acc$n_accrual), c(61, 269, 270))
  expect_true(all(acc$n_accrual == round(acc$n_accrual)))
  cum_diff <- t(apply(acc$n_accrual, 1, cumsum)) -
    t(apply(acc$n_accrual_exact, 1, cumsum))
  expect_lte(max(abs(cum_diff)), 0.5)
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd_accrual(c(60.5, 539.5), accrual_time = c(0, 12),
                              rate_overall = 50), "n_region")
  expect_error(rcpgsd_accrual(c(60, 500), accrual_time = c(0, 12),
                              rate_overall = 50), "must equal")
  expect_error(rcpgsd_accrual(n_reg, accrual_time = c(0, 12)),
               "Exactly one")
  expect_error(rcpgsd_accrual(n_reg, accrual_time = c(0, 12),
                              n_overall = 600, rate_overall = 50),
               "Exactly one")
  expect_error(rcpgsd_accrual(n_reg, accrual_time = c(0, 6, 12),
                              n_overall = 600), "n_overall")
  expect_error(rcpgsd_accrual(n_reg, accrual_time = c(0, 6, 12),
                              rate_overall = -1), "rate_overall")
  expect_error(rcpgsd_accrual(n_reg, accrual_time = c(1, 12),
                              rate_overall = 50), "accrual_time")
  expect_error(rcpgsd_accrual(n_reg, start = c(12, 0, 0),
                              accrual_time = c(0, 12), rate_overall = 50),
               "start")
  expect_error(rcpgsd_accrual(n_reg, start = c(1, 0),
                              accrual_time = c(0, 12), rate_overall = 50),
               "start")
})
