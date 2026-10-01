test_that("print.rcpgsd_accrual prints all sections and returns x invisibly", {
  acc <- rcpgsd_accrual(c(60, 270, 270), start = c(6, 0, 0),
                        accrual_time = c(0, 12), rate_overall = 50)
  expect_output(out <- print(acc), "Regional Accrual")
  expect_identical(out, acc)
  expect_output(print(acc), "Start by region   : \\(6, 0, 0\\)")
  expect_output(print(acc), "\\[6, 12\\)")
  expect_output(print(acc), "Overall")
  expect_output(print(acc), "Accrual rate by interval")
})
