# Helper: extract a single value for a method from a data frame
get_m <- function(df, m, col) {
  v <- df[[col]][df[["method"]] == m]
  expect_true(length(v) >= 1L)
  v
}

# Scenario A used below (three analyses, binding futility, delayed
# enrolment in region 1, region-specific hazard ratios). Expected values
# were computed independently in Python (scipy) with adaptive quadrature.
design_A <- list(upper = c(2.8, 2.4, 2.1), lower = c(0.0, 0.5),
                 binding = TRUE, alpha = 0.025)
events_A <- c(150, 300, 450)
fs_A <- cbind(c(0.05, 0.45, 0.50),
              c(0.08, 0.44, 0.48),
              c(0.10, 0.43, 0.47))
hr_A <- c(0.75, 0.70, 0.72)

test_that("fixed design reproduces the published RCPs of Homma (2026)", {
  # Homma (2026), Pharmaceutical Statistics 25:e70111, Section 3 with K = 1:
  # N = 785, delta = 0.2, sigma = 1, r = 1, S = 3. The normal model with
  # information N / 4 corresponds to 785 events and HR = exp(-0.2).
  z <- qnorm(0.975)
  expected <- data.frame(f1 = c(0.1, 0.1, 0.2, 0.2),
                         hr = c(exp(-0.2), 1, exp(-0.2), 1),
                         M1 = c(0.699, 0.651, 0.780, 0.720),
                         M2 = c(0.832, 0.742, 0.906, 0.812))
  for (i in seq_len(nrow(expected))) {
    f1  <- expected$f1[i]
    res <- rcpgsd(list(upper = z), events = 785,
                  fs = c(f1, (1 - f1) / 2, (1 - f1) / 2),
                  hr = expected$hr[i])
    expect_equal(round(get_m(res$overall, "M1", "rcp"), 3), expected$M1[i])
    expect_equal(round(get_m(res$overall, "M2", "rcp"), 3), expected$M2[i])
  }
})

test_that("fixed design reproduces Hayashi and Itoh (2018)", {
  # Hayashi and Itoh (2018), Appendix B: HR = 0.6, E = 161, E_JP = 30,
  # theta = 0.5, pi = 0.5; numerical integration gives 0.8009655 and the
  # alternative (delta-method) approximation gives 0.799.
  z <- qnorm(0.975)
  res_exact <- rcpgsd(list(upper = z), events = 161,
                      fs = c(30, 131) / 161, hr = 0.6, method = "M1",
                      scale = "hr", hr_method = "exact")
  res_delta <- rcpgsd(list(upper = z), events = 161,
                      fs = c(30, 131) / 161, hr = 0.6, method = "M1",
                      scale = "hr", hr_method = "delta")
  expect_equal(res_exact$stage$rcp_uncond, 0.8009655, tolerance = 1e-6)
  expect_equal(round(res_delta$stage$rcp_uncond, 3), 0.799)
})

test_that("fixed design agrees with the closed forms of Homma (2024)", {
  # Formulas (5), (6) and the joint probability behind (10) of Homma
  # (2024), evaluated here with pnorm() and integrate() only
  d  <- 400
  I  <- d / 4
  f  <- c(0.15, 0.35, 0.50)
  hr <- 0.75
  g  <- 0.5
  z  <- qnorm(0.975)
  dl <- -log(hr)
  lam <- dl * sqrt(I)
  unc1 <- pnorm(sqrt(I) * dl * (1 - g) / sqrt(1 / f[1] + g * (g - 2)))
  unc2 <- prod(pnorm(sqrt(f * I) * dl))
  v_d  <- (1 / f[1] + g * (g - 2)) / I
  c_dz <- (1 - g) / sqrt(I)
  joint1 <- integrate(function(u) {
    dnorm(u - lam) *
      pnorm(((1 - g) * dl + c_dz * (u - lam)) / sqrt(v_d - c_dz^2))
  }, z, Inf, rel.tol = 1e-10)$value

  res <- rcpgsd(list(upper = z), events = d, fs = f, hr = hr)
  expect_equal(get_m(res$stage, "M1", "rcp_uncond"), unc1, tolerance = 1e-8)
  expect_equal(get_m(res$stage, "M2", "rcp_uncond"), unc2, tolerance = 1e-6)
  expect_equal(get_m(res$stage, "M1", "p_joint"), joint1, tolerance = 1e-6)
  expect_equal(get_m(res$stage, "M1", "p_eff"), pnorm(lam - z),
               tolerance = 1e-10)
})

test_that("group sequential results match independent reference values", {
  res <- rcpgsd(design_A, events_A, fs_A, hr_A)
  expect_equal(get_m(res$stage, "M1", "p_eff"),
               c(0.2347782, 0.4710712, 0.2171700), tolerance = 1e-4)
  expect_equal(get_m(res$stage, "M1", "p_fut"),
               c(0.0189101, 0.0043910, 0), tolerance = 1e-4)
  expect_equal(get_m(res$stage, "M1", "p_joint"),
               c(0.1464548, 0.3066338, 0.1412061), tolerance = 1e-4)
  expect_equal(get_m(res$stage, "M1", "rcp_uncond"),
               c(0.5654699, 0.6180289, 0.6614594), tolerance = 1e-6)
  expect_equal(get_m(res$stage, "M2", "p_joint"),
               c(0.1774369, 0.3746013, 0.1718668), tolerance = 5e-4)
  expect_equal(get_m(res$stage, "M2", "rcp_uncond"),
               c(0.5595496, 0.7259987, 0.8202545), tolerance = 1e-4)
  expect_equal(get_m(res$overall, "M1", "power"), 0.9230194,
               tolerance = 1e-4)
  expect_equal(get_m(res$overall, "M1", "rcp"), 0.6438594,
               tolerance = 1e-4)
  expect_equal(get_m(res$overall, "M2", "rcp"), 0.7842792,
               tolerance = 5e-4)
})

test_that("hazard reduction scale matches independent reference values", {
  res_ex <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M1",
                   scale = "hr", hr_method = "exact")
  res_dm <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M1",
                   scale = "hr", hr_method = "delta")
  res_m2 <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M2",
                   scale = "hr", gamma_M2 = 0.05)
  expect_equal(res_ex$stage$p_joint,
               c(0.1511802, 0.3141100, 0.1438625), tolerance = 1e-4)
  expect_equal(res_ex$stage$rcp_uncond,
               c(0.5748784, 0.6326266, 0.6797401), tolerance = 1e-6)
  expect_equal(res_dm$stage$p_joint,
               c(0.1504370, 0.3138725, 0.1436447), tolerance = 1e-4)
  expect_equal(res_dm$stage$rcp_uncond,
               c(0.5731453, 0.6311311, 0.6784374), tolerance = 1e-6)
  expect_equal(res_m2$stage$p_joint,
               c(0.1713199, 0.3513126, 0.1544083), tolerance = 5e-4)
  expect_equal(res_m2$stage$rcp_uncond,
               c(0.4966546, 0.6569710, 0.7559045), tolerance = 1e-4)
})

test_that("stage-specific thresholds are applied stage by stage", {
  res <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M1",
                gamma_M1 = c(0.6, 0.55, 0.5))
  expect_equal(res$stage$gamma, c(0.6, 0.55, 0.5))
  expect_equal(res$stage$p_joint,
               c(0.1394801, 0.2984504, 0.1412061), tolerance = 1e-4)
})

test_that("log hazard ratio and hazard ratio scales coincide for M2 with 0", {
  res_l <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M2")
  res_h <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M2",
                  scale = "hr")
  expect_equal(res_l$stage$p_joint, res_h$stage$p_joint, tolerance = 1e-10)
})

test_that("closed form agrees with Monte Carlo simulation", {
  for (sc in c("loghr", "hr")) {
    res_f <- rcpgsd(design_A, events_A, fs_A, hr_A, scale = sc)
    res_s <- rcpgsd(design_A, events_A, fs_A, hr_A, scale = sc,
                    approach = "simulation", nsim = 200000, seed = 123)
    # Monte Carlo standard error is at most 0.5 / sqrt(200000) = 0.0011;
    # absolute differences are compared against about 4.5 standard errors
    for (col in c("p_eff", "p_fut", "p_joint", "rcp_uncond")) {
      expect_lt(max(abs(res_s$stage[[col]] - res_f$stage[[col]])), 0.005)
    }
    expect_true(all(res_s$stage$p_joint > 0))
  }
})

test_that("output structure is consistent", {
  res <- rcpgsd(design_A, events_A, fs_A, hr_A)
  expect_s3_class(res, "rcpgsd")
  expect_equal(nrow(res$stage), 6L)
  expect_equal(res$overall$method, c("M1", "M2"))
  expect_equal(res$stage$rcp, res$stage$p_joint / res$stage$p_eff)
  expect_true(all(res$stage$p_joint <= res$stage$p_eff))
  expect_equal(res$stage$events_r1[1:3], c(7.5, 24, 45))
  for (m in c("M1", "M2")) {
    sub <- res$stage[res$stage$method == m, ]
    expect_equal(get_m(res$overall, m, "power"), sum(sub$p_eff))
    expect_equal(get_m(res$overall, m, "rcp"),
                 sum(sub$p_joint) / sum(sub$p_eff))
  }
})

test_that("non-binding futility can be obeyed or ignored", {
  design_nb <- design_A
  design_nb$binding <- FALSE
  res_obey <- rcpgsd(design_nb, events_A, fs_A, hr_A, method = "M1")
  res_ign  <- rcpgsd(design_nb, events_A, fs_A, hr_A, method = "M1",
                     futility = "ignore")
  expect_equal(res_obey$design$futility, "non-binding (obeyed)")
  expect_equal(res_ign$design$futility, "non-binding (ignored)")
  expect_equal(res_ign$stage$lower, rep(-Inf, 3))
  expect_equal(res_ign$stage$p_fut, rep(0, 3))
  expect_equal(res_ign$stage$p_eff,
               c(0.2347782, 0.4720514, 0.2267290), tolerance = 1e-4)
  expect_equal(res_ign$stage$p_joint,
               c(0.1464548, 0.3073478, 0.1479996), tolerance = 1e-4)
  expect_equal(res_obey$stage$p_joint,
               c(0.1464548, 0.3066338, 0.1412061), tolerance = 1e-4)
})

test_that("binding futility cannot be ignored", {
  expect_warning(
    res <- rcpgsd(design_A, events_A, fs_A, hr_A, method = "M1",
                  futility = "ignore"),
    "binding"
  )
  expect_equal(res$design$futility, "binding")
  expect_equal(res$stage$lower, c(0, 0.5, -Inf))
})

test_that("a numeric design is equivalent to a list design", {
  res_n <- rcpgsd(c(2.8, 2.4, 2.1), events_A, fs_A, hr_A, method = "M1")
  res_l <- rcpgsd(list(upper = c(2.8, 2.4, 2.1)), events_A, fs_A, hr_A,
                  method = "M1")
  expect_equal(res_n$stage$p_joint, res_l$stage$p_joint)
})

test_that("a timing mismatch triggers a warning", {
  design_t <- list(upper = c(2.8, 2.4, 2.1), timing = c(0.5, 0.75, 1))
  expect_warning(rcpgsd(design_t, events_A, fs_A, hr_A, method = "M1"),
                 "information fractions")
})

test_that("invalid inputs are rejected", {
  expect_error(rcpgsd(design_A, c(150, 300), fs_A, hr_A), "events")
  expect_error(rcpgsd(design_A, c(300, 150, 450), fs_A, hr_A), "increasing")
  expect_error(rcpgsd(design_A, events_A, c(0.1, 0.5, 0.5), hr_A), "sum")
  expect_error(rcpgsd(design_A, events_A, c(0.1, 0.9), c(0.7, 0.8, 0.9)),
               "hr")
  expect_error(rcpgsd(design_A, events_A,
                      cbind(c(0.2, 0.8), c(0.05, 0.95), c(0.1, 0.9)), 0.7),
               "non-decreasing")
  expect_error(rcpgsd(design_A, events_A, fs_A), "hr")
  expect_error(rcpgsd(design_A, events_A, fs_A, hr_A, gamma_M1 = 1), "gamma_M1")
  expect_error(rcpgsd(design_A, events_A, fs_A, hr_A, gamma_M1 = c(0.5, 0.5)),
               "gamma_M1")
  expect_error(rcpgsd(design_A, events_A, fs_A, hr_A, scale = "hr",
                      gamma_M2 = 1), "gamma_M2")
  expect_error(rcpgsd(design_A, events_A, 1, hr_A), "two regions")
  expect_error(rcpgsd("design", events_A, fs_A, hr_A), "design")
})
