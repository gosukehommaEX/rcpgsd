res_plot <- rcpgsd(list(upper = c(3.7103, 2.5114, 1.9930)),
                   events = c(150, 300, 450), fs = c(0.1, 0.45, 0.45),
                   hr = 0.72)

test_that("plot.rcpgsd returns a ggplot that can be built", {
  for (xa in c("stage", "info_frac", "events")) {
    p <- plot(res_plot, x_axis = xa)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p2 <- plot(res_plot, what = "prob")
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("plot.rcpgsd options change the layers", {
  p_full <- plot(res_plot)
  p_min  <- plot(res_plot, show_uncond = FALSE, show_overall = FALSE)
  expect_equal(length(p_full$layers), 3L)
  expect_equal(length(p_min$layers), 2L)
  expect_equal(nrow(p_min$data), 6L)
  expect_equal(nrow(p_full$data), 12L)
})

test_that("plot.rcpgsd works with a single method", {
  res1 <- rcpgsd(list(upper = c(2.963, 1.969)), events = c(200, 400),
                 fs = c(0.15, 0.85), hr = 0.7, method = "M2")
  p <- plot(res1)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot.rcpgsd works for rcpgsd_sim objects", {
  sim <- rcpgsd_sim(list(upper = c(2.963, 1.969)), events = c(100, 200),
                    hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
                    n_accrual = c(40, 260), nsim = 100)
  p <- plot(sim)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot(sim, what = "prob")
  expect_no_error(ggplot2::ggplot_build(p2))
})
