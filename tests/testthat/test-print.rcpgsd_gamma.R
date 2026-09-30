test_that("print.rcpgsd_gamma prints all sections and returns x invisibly", {
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85))
  expect_output(out <- print(adj), "Adjusted Consistency Thresholds")
  expect_identical(out, adj)
  expect_output(print(adj), "gamma_M1_adj")
  expect_output(print(adj), "gamma_M2_adj")
  expect_output(print(adj), "Overall null RCP")
})

test_that("print.rcpgsd_gamma handles a single method", {
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85),
                          method = "M2")
  expect_output(print(adj), "gamma_M2_adj")
  expect_false(any(grepl("gamma_M1_adj", capture.output(print(adj)))))
})
