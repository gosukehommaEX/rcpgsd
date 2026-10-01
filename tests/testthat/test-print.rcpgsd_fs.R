test_that("print.rcpgsd_fs prints all sections and returns x invisibly", {
  req <- rcpgsd_get_fs(list(upper = c(2.963, 1.969)), events = c(200, 400),
                       fs = c(0.15, 0.85), hr = 0.7, method = "M1")
  expect_output(out <- print(req), "Required Proportion of Events")
  expect_identical(out, req)
  expect_output(print(req), "Method 1 at the required proportion")
  expect_output(print(req), "overall RCP given success")
})

test_that("print.rcpgsd_fs reports the lower limit and missing solutions", {
  req <- rcpgsd_get_fs(list(upper = c(2.963, 1.969)), events = c(200, 400),
                       fs = c(0.15, 0.85), hr = 0.7, method = "M1",
                       target = 0.3)
  expect_output(print(req), "lower limit")
  # The Method 2 RCP does not exceed about 0.996 for any proportion, so the
  # target 0.999 is not reached (checked in Python)
  req2 <- suppressWarnings(
    rcpgsd_get_fs(list(upper = c(2.963, 1.969)), events = c(200, 400),
                  fs = c(0.15, 0.85), hr = 0.7, method = "M2",
                  target = 0.999)
  )
  expect_true(is.na(req2$result$fs_r1))
  expect_output(print(req2), "NA")
})
