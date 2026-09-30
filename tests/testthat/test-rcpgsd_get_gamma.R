# Scenario B: three analyses with O'Brien-Fleming-type boundaries,
# S = 3 regions with fs = (0.1, 0.45, 0.45), 785 events at the final
# analysis. Expected values were computed independently in Python (scipy)
# with Brent's root-finding method.
design_B <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
events_B <- c(262, 523, 785)
fs_B     <- c(0.1, 0.45, 0.45)
adj_B    <- rcpgsd_get_gamma(design_B, events_B, fs_B)

test_that("reference null RCP equals the fixed-design value", {
  adj <- adj_B
  # Homma (2026) reports 0.651 (Method 1) and 0.742 (Method 2) for K = 1
  expect_equal(unname(adj$ref[["M1"]]), 0.6513492, tolerance = 1e-5)
  expect_equal(unname(adj$ref[["M2"]]), 0.74198, tolerance = 2e-4)
})

test_that("unadjusted null RCP is inflated at early analyses", {
  adj <- adj_B
  expect_equal(adj$result$p_eff0, c(0.0001035, 0.0059452, 0.0189569),
               tolerance = 1e-4)
  expect_equal(adj$result$rcp0_M1_unadj, c(0.7446936, 0.6811454, 0.6497045),
               tolerance = 1e-4)
  expect_equal(adj$result$rcp0_M2_unadj, c(0.9051884, 0.8138122, 0.7393343),
               tolerance = 5e-4)
  expect_true(adj$result$rcp0_M1_unadj[1] > adj$ref[["M1"]])
  expect_true(adj$result$rcp0_M2_unadj[1] > adj$ref[["M2"]])
})

test_that("stagewise adjusted thresholds restore the reference value", {
  adj <- adj_B
  expect_equal(adj$result$gamma_M1_adj, c(0.704522, 0.587161, 0.5),
               tolerance = 2e-3)
  expect_equal(adj$result$gamma_M2_adj[1:2], c(0.20212, 0.04392),
               tolerance = 5e-3)
  # Stage 3 is not inflated, so the unadjusted thresholds are kept
  expect_equal(adj$result$gamma_M1_adj[3], 0.5)
  expect_equal(adj$result$gamma_M2_adj[3], 0)
  expect_equal(adj$result$rcp0_M1_adj[1:2], rep(adj$ref[["M1"]], 2),
               tolerance = 1e-4)
  expect_equal(adj$result$rcp0_M2_adj[1:2], rep(adj$ref[["M2"]], 2),
               tolerance = 5e-4)
})

test_that("relaxation is possible only when requested", {
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1",
                          allow_relax = TRUE)
  expect_equal(adj$result$gamma_M1_adj[3], 0.494210, tolerance = 2e-3)
  expect_true(adj$result$gamma_M1_adj[3] < 0.5)
  # For Method 2 the solution at stage 3 would be negative
  expect_warning(
    adj2 <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M2",
                             allow_relax = TRUE),
    "No adjusted threshold"
  )
  expect_true(is.na(adj2$result$gamma_M2_adj[3]))
})

test_that("overall adjustment gives a single threshold", {
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, type = "overall")
  expect_equal(adj$result$gamma_M1_adj, rep(0.520876, 3), tolerance = 2e-3)
  expect_equal(adj$result$gamma_M2_adj, rep(0.00685, 3), tolerance = 0.1)
  expect_equal(adj$overall$rcp0_unadj, c(0.6575729, 0.7577615),
               tolerance = 5e-4)
  expect_equal(adj$overall$rcp0_adj, unname(adj$ref), tolerance = 5e-4)
})

test_that("adjusted thresholds can be passed to rcpgsd()", {
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1")
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
  adj <- rcpgsd_get_gamma(design_B, events_B, fs_B, method = "M1",
                          scale = "hr")
  expect_true(adj$result$gamma_M1_adj[1] > 0.5)
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
