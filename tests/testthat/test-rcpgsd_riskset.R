test_that("risk sets are counted correctly in a small example", {
  # Sorted times 1-5; events at times 1 (experimental), 3 (control) and
  # 4 (experimental)
  rs <- .rcpgsd_riskset(time   = c(4, 2, 5, 1, 3),
                        status = c(1, 0, 0, 1, 1),
                        trt    = c(1, 0, 1, 1, 0))
  expect_equal(rs$n_t, c(3, 2, 2))
  expect_equal(rs$n_c, c(2, 1, 0))
  expect_equal(rs$d_t, c(1, 0, 1))
})

test_that("empty data give empty risk sets", {
  rs <- .rcpgsd_riskset(numeric(0), integer(0), integer(0))
  expect_length(rs$n_t, 0L)
  expect_length(rs$d_t, 0L)
})
