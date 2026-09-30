test_that("Gauss-Legendre rule integrates polynomials exactly", {
  gl <- .rcpgsd_gauss_legendre(8, -1, 2)
  expect_length(gl$x, 8L)
  expect_true(all(diff(gl$x) > 0))
  expect_equal(sum(gl$w), 3, tolerance = 1e-12)
  # A rule with n nodes is exact for polynomials of degree 2n - 1 = 15
  expect_equal(sum(gl$w * gl$x^15), (2^16 - 1) / 16, tolerance = 1e-12)
})

test_that("Gauss-Legendre rule integrates the normal density", {
  gl <- .rcpgsd_gauss_legendre(64, -1, 1.5)
  expect_equal(sum(gl$w * dnorm(gl$x)), pnorm(1.5) - pnorm(-1),
               tolerance = 1e-12)
})
