test_that("print.rcpgsd prints all sections and returns x invisibly", {
  res <- rcpgsd(list(upper = c(2.963, 1.969), lower = 0, binding = TRUE),
                events = c(200, 400), fs = c(0.15, 0.85), hr = 0.7)
  expect_output(out <- print(res), "Group Sequential Design")
  expect_identical(out, res)
  expect_output(print(res), "Method 1")
  expect_output(print(res), "Method 2")
  expect_output(print(res), "binding")
  expect_output(print(res), "Overall")
})

test_that("print.rcpgsd reports the approach and the scale", {
  res_s <- rcpgsd(list(upper = c(2.963, 1.969)), events = c(200, 400),
                  fs = c(0.15, 0.85), hr = 0.7, method = "M1",
                  scale = "hr", approach = "simulation", nsim = 1000)
  expect_output(print(res_s), "Monte Carlo simulation")
  expect_output(print(res_s), "hazard reduction")
  res_d <- rcpgsd(list(upper = c(2.963, 1.969)), events = c(200, 400),
                  fs = c(0.15, 0.85), hr = 0.7, method = "M1",
                  scale = "hr", hr_method = "delta")
  expect_output(print(res_d), "delta method")
})
