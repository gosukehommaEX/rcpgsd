# Expected values for group sequential designs were computed independently
# in Python (scipy) with Brent's root-finding method.
design_F <- list(upper = c(2.963, 1.969), alpha = 0.025)

test_that("Table 1 of Hayashi and Itoh (2018) is reproduced", {
  # Hazard reduction scale, unconditional probability 0.8, pi = 0.5,
  # numerical integration; required numbers of Japanese events
  tab <- data.frame(hr     = c(0.8, 0.7, 0.6, 0.5),
                    events = c(844, 330, 161, 87),
                    e_jp   = c(173.5, 64.6, 29.8, 15.2))
  for (i in seq_len(nrow(tab))) {
    req <- rcpgsd_get_fs(list(upper = qnorm(0.975)), events = tab$events[i],
                         fs = c(0.2, 0.8), hr = tab$hr[i], method = "M1",
                         type = "unconditional", scale = "hr")
    expect_equal(round(req$result$events_r1, 1), tab$e_jp[i])
  }
})

test_that("the log hazard ratio scale gives a constant fraction", {
  # Hayashi and Itoh (2018), Section 6, and formula (8) of Homma (2024):
  # with 90% power at two-sided 5%, the required fraction is about 22.4%
  # irrespective of the hazard ratio
  zz <- (qnorm(0.975) + qnorm(0.9)) * 0.5 / qnorm(0.8)
  w  <- 1 / (zz^2 + 0.75)
  for (hr in c(0.6, 0.8)) {
    d <- 4 * (qnorm(0.975) + qnorm(0.9))^2 / log(hr)^2
    req <- rcpgsd_get_fs(list(upper = qnorm(0.975)), events = d,
                         fs = c(0.2, 0.8), hr = hr, method = "M1",
                         type = "unconditional")
    expect_equal(req$result$fs_r1, w, tolerance = 1e-5)
  }
  expect_equal(round(w, 3), 0.224)
})

test_that("group sequential results match independent values", {
  req_o <- rcpgsd_get_fs(design_F, events = c(200, 400), fs = c(0.15, 0.85),
                         hr = 0.7)
  req_s <- rcpgsd_get_fs(design_F, events = c(200, 400), fs = c(0.15, 0.85),
                         hr = 0.7, type = "stagewise")
  expect_equal(req_o$result$fs_r1, c(0.198365, 0.058259), tolerance = 1e-3)
  expect_equal(req_s$result$fs_r1, c(0.209081, 0.061991), tolerance = 1e-3)
  # The criterion equals the target at the solution
  expect_equal(req_o$fit$M1$overall$rcp, 0.8, tolerance = 1e-4)
  expect_equal(min(req_s$fit$M1$stage$rcp), 0.8, tolerance = 1e-4)
  expect_true(all(req_o$result$value_ceiling >= 0.8))
  expect_equal(req_o$result$events_r1, req_o$result$fs_r1 * 400)
})

test_that("the shape of the template is kept", {
  tm  <- cbind(c(0.05, 0.95), c(0.15, 0.85))
  req_o <- rcpgsd_get_fs(design_F, events = c(200, 400), fs = tm, hr = 0.7,
                         method = "M1")
  req_s <- rcpgsd_get_fs(design_F, events = c(200, 400), fs = tm, hr = 0.7,
                         method = "M1", type = "stagewise")
  expect_equal(req_o$result$fs_r1, 0.234651, tolerance = 1e-3)
  expect_equal(req_s$result$fs_r1, 0.540028, tolerance = 1e-3)
  f <- req_o$fit$M1$fs
  expect_equal(f[1, 1] / f[1, 2], 0.05 / 0.15)
  expect_equal(colSums(f), c(1, 1))
  # The search stays where the other region's events do not decrease:
  # 0.95 * 200 / 0.95 = 200 and 0.85 * 400 / 0.85 = 400 give
  # x <= (400 - 200) / (400 - 200 / 3) = 0.6
  ev_s <- sweep(req_s$fit$M1$fs, 2, c(200, 400), "*")
  expect_true(ev_s[2, 2] >= ev_s[2, 1])
  expect_lt(req_s$result$fs_r1, 0.6)
})

test_that("limits of the search interval are reported", {
  req <- rcpgsd_get_fs(design_F, events = c(200, 400), fs = c(0.15, 0.85),
                       hr = 0.7, method = "M1", target = 0.3)
  expect_true(req$result$at_lower)
  expect_equal(req$result$fs_r1, 0.01)
  expect_warning(
    req2 <- rcpgsd_get_fs(design_F, events = c(200, 400),
                          fs = c(0.15, 0.85), hr = 0.7, method = "M2",
                          target = 0.999),
    "not reached"
  )
  expect_true(is.na(req2$result$fs_r1))
  expect_null(req2$fit$M2)
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd_get_fs(design_F, c(200, 400), c(0.15, 0.85), hr = 0.7,
                             target = 1), "target")
  expect_error(rcpgsd_get_fs(design_F, c(200, 400), c(0.15, 0.85), hr = 0.7,
                             interval = c(0.5, 0.2)), "interval")
  expect_error(rcpgsd_get_fs(design_F, c(200, 400), c(0.15, 0.85)), "hr")
})
