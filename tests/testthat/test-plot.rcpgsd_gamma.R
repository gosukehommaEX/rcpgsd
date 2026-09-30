test_that("plot.rcpgsd_gamma returns a ggplot that can be built", {
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85))
  p <- plot(adj)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  # 2 methods x (1 threshold series + 2 null RCP series) x 2 stages
  expect_equal(nrow(p$data), 12L)
})

test_that("plot.rcpgsd_gamma works with a single method", {
  adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
                          events = c(200, 400), fs = c(0.15, 0.85),
                          method = "M1")
  p <- plot(adj)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(nrow(p$data), 6L)
})
