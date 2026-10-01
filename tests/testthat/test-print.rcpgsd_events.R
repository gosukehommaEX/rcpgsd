test_that("print.rcpgsd_events prints all sections and returns x invisibly", {
  ev <- rcpgsd_events(events = c(100, 200), hr = 0.7, lambda_c = 0.05,
                      accrual_time = c(0, 12), n_accrual = c(40, 260))
  expect_output(out <- print(ev), "Expected Analysis Times")
  expect_identical(out, ev)
  expect_output(print(ev), "Region 1")
  expect_output(print(ev), "Stage 2")
  expect_output(print(ev), "True HR           : \\(0.7, 0.7\\)")
})
