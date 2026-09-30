test_that("list designs are standardized", {
  b <- .rcpgsd_bounds(list(upper = c(3, 2.5, 2), lower = c(0, -10),
                           binding = TRUE, alpha = 0.025))
  expect_equal(b$K, 3L)
  expect_equal(b$upper, c(3, 2.5, 2))
  expect_equal(b$lower, c(0, -Inf, -Inf))
  expect_true(b$binding)
  expect_true(b$has_futility)
  expect_equal(b$alpha, 0.025)
  expect_equal(b$source, "list")

  b2 <- .rcpgsd_bounds(list(upper = c(3, 2)))
  expect_equal(b2$lower, c(-Inf, -Inf))
  expect_false(b2$binding)
  expect_false(b2$has_futility)
  expect_true(is.na(b2$alpha))

  b3 <- .rcpgsd_bounds(c(3, 2))
  expect_equal(b3$upper, c(3, 2))
  expect_equal(b3$source, "numeric")
})

test_that("invalid designs are rejected", {
  expect_error(.rcpgsd_bounds(list(lower = 0)), "upper")
  expect_error(.rcpgsd_bounds(list(upper = c(3, 2), lower = c(3.5, 0))),
               "futility")
  expect_error(.rcpgsd_bounds(list(upper = c(3, 2, 1.9), lower = 1:4)),
               "length")
  expect_error(.rcpgsd_bounds("a"), "design")
})

test_that("gsDesign objects are converted", {
  message("gsDesign available for testing: ",
          requireNamespace("gsDesign", quietly = TRUE))
  skip_if_not_installed("gsDesign")
  gs4 <- gsDesign::gsDesign(k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
                            sfu = gsDesign::sfLDOF, sfl = gsDesign::sfHSD,
                            sflpar = -2)
  b4 <- .rcpgsd_bounds(gs4)
  expect_equal(b4$upper, gs4$upper$bound)
  expect_equal(b4$lower[1:2], gs4$lower$bound[1:2])
  expect_equal(b4$lower[3], -Inf)
  expect_false(b4$binding)
  expect_equal(b4$alpha, 0.025)
  expect_equal(b4$timing, gs4$timing)
  expect_equal(b4$source, "gsDesign")

  gs3 <- gsDesign::gsDesign(k = 3, test.type = 3, alpha = 0.025, beta = 0.1,
                            sfu = gsDesign::sfLDOF, sfl = gsDesign::sfHSD,
                            sflpar = -2)
  expect_true(.rcpgsd_bounds(gs3)$binding)

  gs1 <- gsDesign::gsDesign(k = 3, test.type = 1, alpha = 0.025, beta = 0.1,
                            sfu = gsDesign::sfLDOF)
  b1 <- .rcpgsd_bounds(gs1)
  expect_equal(b1$lower, rep(-Inf, 3))
  expect_false(b1$has_futility)
})

test_that("rpact objects are converted", {
  message("rpact available for testing: ",
          requireNamespace("rpact", quietly = TRUE))
  skip_if_not_installed("rpact")
  rp <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025, beta = 0.1,
                                        typeOfDesign = "asOF",
                                        futilityBounds = c(0, 0.5),
                                        bindingFutility = TRUE)
  b <- .rcpgsd_bounds(rp)
  expect_equal(b$upper, as.numeric(rp$criticalValues))
  expect_equal(b$lower, c(0, 0.5, -Inf))
  expect_true(b$binding)
  expect_equal(b$alpha, 0.025)
  expect_equal(b$timing, as.numeric(rp$informationRates))
  expect_equal(b$source, "rpact")

  rp0 <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025,
                                         typeOfDesign = "asOF")
  b0 <- .rcpgsd_bounds(rp0)
  expect_equal(b0$lower, rep(-Inf, 3))
  expect_false(b0$has_futility)
})

test_that("gsDesign and rpact boundaries give the same RCPs", {
  avail <- c(gsDesign = requireNamespace("gsDesign", quietly = TRUE),
             rpact    = requireNamespace("rpact", quietly = TRUE))
  message("Packages available for cross-check: ",
          paste(names(avail), avail, sep = " = ", collapse = ", "))
  skip_if_not_installed("gsDesign")
  skip_if_not_installed("rpact")
  gs <- gsDesign::gsDesign(k = 3, test.type = 1, alpha = 0.025, beta = 0.1,
                           sfu = gsDesign::sfLDOF)
  rp <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025,
                                        typeOfDesign = "asOF")
  expect_equal(.rcpgsd_bounds(gs)$upper, .rcpgsd_bounds(rp)$upper,
               tolerance = 1e-3)
  res_gs <- rcpgsd(gs, events = c(150, 300, 450), fs = c(0.1, 0.45, 0.45),
                   hr = 0.72)
  res_rp <- rcpgsd(rp, events = c(150, 300, 450), fs = c(0.1, 0.45, 0.45),
                   hr = 0.72)
  expect_equal(res_gs$overall$rcp, res_rp$overall$rcp, tolerance = 1e-3)
})
