# Scenario B: three analyses with O'Brien-Fleming-type boundaries,
# S = 3 regions with fs = (0.1, 0.45, 0.45), 785 events at the final
# analysis (Method 1). Scenario 2: two analyses, S = 2 regions with
# fs = (0.15, 0.85), 400 events at the final analysis (Method 2, which is
# more expensive to evaluate). Expected values were computed independently
# in Python (scipy) with Brent's root-finding method.
design_B <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
events_B <- c(262, 523, 785)
fs_B     <- c(0.1, 0.45, 0.45)
adj_B    <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1")
design_2 <- list(upper = c(2.963, 1.969), alpha = 0.025)
adj_2    <- rcpgsd_get_gamma(design_2, c(200, 400), c(0.15, 0.85),
                             method = "M2")

test_that("reference null RCP equals the fixed-design value", {
  # Homma (2026) reports 0.651 for Method 1 with K = 1; the Method 2 value
  # 0.742 of Homma (2026) is checked in test-rcpgsd.R
  expect_equal(unname(adj_B$ref[["M1"]]), 0.6513492, tolerance = 1e-5)
  expect_equal(unname(adj_2$ref[["M2"]]), 0.8345168, tolerance = 2e-4)
})

test_that("unadjusted null RCP is inflated at early analyses", {
  expect_equal(adj_B$result$p_eff0, c(0.0001035, 0.0059452, 0.0189569),
               tolerance = 1e-4)
  expect_equal(adj_B$result$rcp0_M1_unadj,
               c(0.7446936, 0.6811454, 0.6497045), tolerance = 1e-4)
  expect_equal(adj_2$result$p_eff0, c(0.0015233, 0.0234529),
               tolerance = 1e-4)
  expect_equal(adj_2$result$rcp0_M2_unadj, c(0.9125367, 0.8340365),
               tolerance = 5e-4)
  expect_true(adj_B$result$rcp0_M1_unadj[1] > adj_B$ref[["M1"]])
  expect_true(adj_2$result$rcp0_M2_unadj[1] > adj_2$ref[["M2"]])
})

test_that("stagewise adjusted thresholds restore the reference value", {
  expect_equal(adj_B$result$gamma_M1_adj, c(0.704522, 0.587161, 0.5),
               tolerance = 2e-3)
  expect_equal(adj_2$result$gamma_M2_adj[1], 0.130102, tolerance = 5e-3)
  # The last analysis is not inflated, so the unadjusted thresholds are kept
  expect_equal(adj_B$result$gamma_M1_adj[3], 0.5)
  expect_equal(adj_2$result$gamma_M2_adj[2], 0)
  expect_equal(adj_B$result$rcp0_M1_adj[1:2], rep(adj_B$ref[["M1"]], 2),
               tolerance = 1e-4)
  expect_equal(adj_2$result$rcp0_M2_adj[1], adj_2$ref[["M2"]],
               tolerance = 5e-4)
})

test_that("relaxation is possible only when requested", {
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1",
                          allow_relax = TRUE)
  expect_equal(adj$result$gamma_M1_adj, c(0.704522, 0.587161, 0.494210),
               tolerance = 2e-3)
  expect_true(adj$result$gamma_M1_adj[3] < 0.5)
  # A reference value that cannot be reached with a non-negative threshold
  # One warning per analysis is expected; all are collected here
  msgs <- character(0)
  adj2 <- withCallingHandlers(
    rcpgsd_get_gamma(design_2, c(200, 400), c(0.15, 0.85), method = "M2",
                     allow_relax = TRUE, rcp_ref = c(M2 = 0.99)),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(msgs, 2L)
  expect_true(all(grepl("No adjusted threshold", msgs)))
  expect_true(all(is.na(adj2$result$gamma_M2_adj)))
  # Without relaxation the unadjusted threshold is kept instead
  adj3 <- rcpgsd_get_gamma(design_2, c(200, 400), c(0.15, 0.85),
                           method = "M2", rcp_ref = c(M2 = 0.99))
  expect_equal(adj3$result$gamma_M2_adj, c(0, 0))
})

test_that("overall adjustment gives a single threshold", {
  # Two analyses (lighter than scenario B); expected values computed
  # independently in Python (scipy)
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85),
                          type = "overall")
  expect_equal(unname(adj$ref), c(0.6878309, 0.8345168), tolerance = 2e-4)
  expect_equal(adj$result$gamma_M1_adj, rep(0.509934, 2), tolerance = 2e-3)
  expect_equal(adj$result$gamma_M2_adj, rep(0.004255, 2), tolerance = 0.1)
  expect_equal(adj$overall$rcp0_unadj, c(0.6913113, 0.8388244),
               tolerance = 5e-4)
  expect_equal(adj$overall$rcp0_adj, unname(adj$ref), tolerance = 5e-4)
})

test_that("adjusted thresholds can be passed to rcpgsd()", {
  adj <- adj_B
  res0 <- rcpgsd(design_B, events_B, fs_B, hr = 1, method = "M1",
                 gamma_M1 = adj$result$gamma_M1_adj)
  expect_equal(res0$stage$rcp, adj$result$rcp0_M1_adj, tolerance = 1e-8)
  res1 <- rcpgsd(design_B, events_B, fs_B, hr = exp(-0.2), method = "M1",
                 gamma_M1 = adj$result$gamma_M1_adj)
  res_u <- rcpgsd(design_B, events_B, fs_B, hr = exp(-0.2), method = "M1")
  expect_true(all(res1$stage$rcp <= res_u$stage$rcp + 1e-8))
})

test_that("a user-specified reference value is honoured", {
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1",
                          rcp_ref = c(M1 = 0.7))
  expect_equal(unname(adj$ref[["M1"]]), 0.7)
  expect_equal(adj$result$rcp0_M1_adj[1], 0.7, tolerance = 1e-4)
  expect_equal(adj$result$gamma_M1_adj[2:3], c(0.5, 0.5))
})

test_that("hazard reduction scale gives a valid adjustment", {
  # Two analyses; expected values computed independently in Python (scipy)
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85),
                          method = "M1", scale = "hr")
  expect_equal(unname(adj$ref[["M1"]]), 0.6979385, tolerance = 1e-4)
  expect_equal(adj$result$gamma_M1_adj, c(0.672895, 0.5), tolerance = 2e-3)
  expect_equal(adj$result$rcp0_M1_adj[1], adj$ref[["M1"]], tolerance = 1e-4)
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd_get_gamma(list(upper = design_B$upper), events_B,
                                fs_B), "significance level")
  expect_error(rcpgsd_get_gamma(design_B, events_B, fs_B,
                                interval_M1 = c(0.6, 0.9)), "interval_M1")
  expect_error(rcpgsd_get_gamma(design_B, events_B, fs_B,
                                gamma_M1 = c(0.5, 0.5)), "gamma_M1")
  expect_error(rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1",
                                rcp_ref = c(M2 = 0.7)), "rcp_ref")
})
