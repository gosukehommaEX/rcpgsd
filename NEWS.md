# rcpgsd 0.1.0

* Initial version.
* `rcpgsd()` computes stagewise and overall regional consistency
  probabilities (Method 1 and Method 2) for group sequential designs with a
  time-to-event endpoint, on the log hazard ratio scale or on the hazard
  reduction scale (exact numerical integration or delta method), by closed
  form or by Monte Carlo simulation. Boundaries can be supplied as `gsDesign`
  or `rpact` objects, as a list, or as a numeric vector; binding and
  non-binding futility boundaries are supported.
* `rcpgsd_get_gamma()` computes stagewise or overall consistency thresholds
  that restore the null RCP given efficacy stopping to its fixed-design
  reference value.
* `rcpgsd_sim()` simulates event-driven trials at the patient level
  (piecewise uniform accrual by region, exponential event times) and
  estimates the same probabilities as `rcpgsd()`, using the log-rank test,
  the Pike or Peto estimate of the overall hazard ratio and regional Cox
  models. It reports the mean regional event proportions, which can be
  passed to `rcpgsd()` for comparison with the closed form, and the
  frequencies of missing or infinite regional estimates.
* `print()` methods for all result classes and `plot()` methods for
  `rcpgsd` (also used for `rcpgsd_sim`) and `rcpgsd_gamma` objects.
