test_that("print.rcpgsd_sim prints all sections and returns x invisibly", {
  sim <- rcpgsd_sim(list(upper = c(2.963, 1.969)), events = c(100, 200),
                    hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
                    n_accrual = c(40, 260), nsim = 100)
  expect_output(out <- print(sim), "Patient-Level Simulation")
  expect_identical(out, sim)
  expect_output(print(sim), "Pike's method")
  expect_output(print(sim), "Pr\\(inf est R1\\)")
  expect_output(print(sim), "Method 1")
  expect_output(print(sim), "Method 2")
  expect_output(print(sim), "Overall")
})
