# Scenario: three analyses with O'Brien-Fleming-type boundaries, 600
# patients, Japan (region 1) enrolling only after month 6
design_S <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
events_S <- c(150, 300, 450)
nacc_S   <- rbind(c(0, 60), c(100, 170), c(100, 170))
lam_S    <- log(2) / 12
sim_S    <- rcpgsd_sim(design_S, events_S, hr = 0.72, lambda_c = lam_S,
                       accrual_time = c(0, 6, 18), n_accrual = nacc_S,
                       nsim = 1000, seed = 2)

test_that("output structure is consistent", {
  sim <- sim_S
  expect_s3_class(sim, "rcpgsd_sim")
  expect_s3_class(sim, "rcpgsd")
  expect_equal(nrow(sim$stage), 6L)
  expect_equal(sim$overall$method, c("M1", "M2"))
  expect_equal(dim(sim$fs_mean), c(3L, 3L))
  # Every simulated trial has exactly events[k] events at analysis k
  expect_equal(colSums(sim$fs_mean), rep(1, 3), tolerance = 1e-12)
  expect_equal(sim$analysis$events_r1_mean, sim$fs_mean[1, ] * events_S)
  expect_equal(sim$stage$events_r1[1:3], sim$analysis$events_r1_mean)
  expect_true(all(diff(sim$analysis$time_mean) > 0))
  expect_true(all(sim$stage$p_joint <= sim$stage$p_eff))
  expect_equal(sim$stage$rcp, sim$stage$p_joint / sim$stage$p_eff)
  for (m in c("M1", "M2")) {
    sub <- sim$stage[sim$stage$method == m, ]
    expect_equal(sim$overall$power[sim$overall$method == m], sum(sub$p_eff))
  }
})

test_that("results are reproducible with a seed", {
  args <- list(design = design_S, events = events_S, hr = 0.72,
               lambda_c = lam_S, accrual_time = c(0, 6, 18),
               n_accrual = nacc_S, nsim = 50)
  sim1 <- do.call(rcpgsd_sim, c(args, seed = 5))
  sim2 <- do.call(rcpgsd_sim, c(args, seed = 5))
  sim3 <- do.call(rcpgsd_sim, c(args, seed = 6))
  expect_identical(sim1$stage, sim2$stage)
  expect_false(identical(sim1$analysis$time_mean, sim3$analysis$time_mean))
  # Zero drop-out reproduces the results without drop-out exactly
  sim4 <- do.call(rcpgsd_sim, c(args, seed = 5,
                                list(lambda_d = matrix(0, 3, 2))))
  expect_identical(sim1$stage, sim4$stage)
  expect_identical(sim1$analysis, sim4$analysis)
  expect_equal(sim1$analysis$p_not_reached, rep(0, 3))
})

test_that("simulation agrees with the closed form", {
  cf <- rcpgsd(design_S, events_S, fs = sim_S$fs_mean, hr = 0.72)
  # Monte Carlo standard error is at most 0.5 / sqrt(1000) = 0.016;
  # absolute differences are compared against about 3.8 standard errors
  for (col in c("p_eff", "p_joint", "rcp_uncond")) {
    expect_lt(max(abs(sim_S$stage[[col]] - cf$stage[[col]])), 0.06)
  }
})

test_that("the consistency scale changes only the Method 1 criterion", {
  args <- list(design = design_S, events = events_S, hr = 0.72,
               lambda_c = lam_S, accrual_time = c(0, 6, 18),
               n_accrual = nacc_S, nsim = 100, seed = 8)
  sim_l <- do.call(rcpgsd_sim, c(args, scale = "loghr"))
  sim_h <- do.call(rcpgsd_sim, c(args, scale = "hr"))
  expect_equal(sim_l$stage$p_eff, sim_h$stage$p_eff)
  m2 <- sim_l$stage$method == "M2"
  expect_equal(sim_l$stage$p_joint[m2], sim_h$stage$p_joint[m2])
  # With HR < 1 the hazard reduction criterion is met at least as often
  m1 <- sim_l$stage$method == "M1"
  expect_true(all(sim_h$stage$rcp_uncond[m1] >= sim_l$stage$rcp_uncond[m1] -
                    0.05))
})

test_that("Pike and Peto estimators give similar results", {
  args <- list(design = design_S, events = events_S, hr = 0.72,
               lambda_c = lam_S, accrual_time = c(0, 6, 18),
               n_accrual = nacc_S, nsim = 150, seed = 3)
  sim_pike <- do.call(rcpgsd_sim, c(args, overall_est = "pike"))
  sim_peto <- do.call(rcpgsd_sim, c(args, overall_est = "peto"))
  expect_equal(sim_pike$stage$p_eff, sim_peto$stage$p_eff)
  expect_lt(max(abs(sim_pike$stage$rcp_uncond - sim_peto$stage$rcp_uncond)),
            0.03)
})

test_that("small regions produce missing or infinite estimates", {
  sim <- rcpgsd_sim(list(upper = c(2.963, 1.969)), events = c(60, 120),
                    hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
                    n_accrual = c(10, 390), nsim = 200, seed = 4)
  expect_true(sim$analysis$p_no_est_r1[1] + sim$analysis$p_inf_r1[1] > 0)
  expect_true(all(sim$p_no_est >= 0 & sim$p_no_est <= 1))
  expect_true(all(sim$p_inf >= 0 & sim$p_inf <= 1))
})

test_that("non-binding futility can be obeyed or ignored", {
  design_f <- list(upper = c(2.963, 1.969), lower = 0.5, binding = FALSE)
  args <- list(design = design_f, events = c(100, 200), hr = 0.8,
               lambda_c = 0.05, accrual_time = c(0, 12),
               n_accrual = c(40, 260), nsim = 150, seed = 7)
  sim_o <- do.call(rcpgsd_sim, c(args, futility = "obey"))
  sim_i <- do.call(rcpgsd_sim, c(args, futility = "ignore"))
  expect_true(sim_o$stage$p_fut[1] > 0)
  expect_equal(sim_i$stage$p_fut, rep(0, 4))
  expect_true(sim_i$stage$p_eff[2] >= sim_o$stage$p_eff[2])
})

test_that("analyses not reached because of drop-out are reported", {
  # The expected total number of events is about 136.8, so the second
  # analysis (135 events) is not reached in a part of the trials
  args <- list(design = list(upper = c(2.963, 1.969)), events = c(100, 135),
               hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
               n_accrual = c(40, 260), lambda_d = 0.05, nsim = 30, seed = 10)
  expect_warning(sim <- do.call(rcpgsd_sim, args), "not reached")
  p_nr <- sim$analysis$p_not_reached
  expect_equal(p_nr[1], 0)
  expect_gt(p_nr[2], 0)
  # Such trials neither stop nor meet the consistency criteria
  st2 <- sim$stage[sim$stage$stage == 2, ]
  expect_true(all(st2$p_eff + st2$p_fut <= 1 - p_nr[2]))
  expect_true(all(st2$rcp_uncond <= 1 - p_nr[2]))
  # Means are taken over the trials that reached the analysis
  expect_equal(colSums(sim$fs_mean), c(1, 1), tolerance = 1e-12)
  expect_true(all(is.finite(sim$analysis$time_mean)))
  expect_equal(sim$lambda_d, matrix(0.05, 2, 2))
})

test_that("invalid inputs are rejected", {
  base <- list(design = design_S, events = events_S, hr = 0.72,
               lambda_c = lam_S, accrual_time = c(0, 6, 18),
               n_accrual = nacc_S, nsim = 10)
  mod <- function(...) {
    a <- base
    a[names(list(...))] <- list(...)
    a
  }
  expect_error(do.call(rcpgsd_sim, mod(events = c(150.5, 300, 450))),
               "events")
  expect_error(do.call(rcpgsd_sim, mod(accrual_time = c(1, 6, 18))),
               "accrual_time")
  expect_error(do.call(rcpgsd_sim, mod(n_accrual = nacc_S[, 1])),
               "n_accrual")
  expect_error(do.call(rcpgsd_sim, mod(events = c(150, 300, 700))),
               "total number of patients")
  expect_error(do.call(rcpgsd_sim, mod(hr = c(0.7, 0.8))), "hr")
  expect_error(do.call(rcpgsd_sim, mod(lambda_c = -1)), "lambda_c")
  expect_error(do.call(rcpgsd_sim, mod(n_accrual = rbind(c(0, 0),
                                                         c(300, 300)))),
               "at least one patient")
  expect_error(do.call(rcpgsd_sim, mod(lambda_d = -0.01)), "lambda_d")
  expect_error(do.call(rcpgsd_sim, mod(lambda_d = matrix(0.01, 2, 2))),
               "lambda_d")
  expect_error(do.call(rcpgsd_sim, mod(lambda_d = 0.05)),
               "expected total number of events")
})
