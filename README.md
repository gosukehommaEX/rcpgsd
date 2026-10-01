# rcpgsd

**Regional Consistency Probability under Group Sequential Designs**

## Overview

`rcpgsd` is an R package for computing regional consistency probabilities (RCPs) in multi-regional clinical trials (MRCTs) with a time-to-event primary endpoint analysed by a group sequential design (GSD).

When a trial stops for efficacy at an interim analysis, regional consistency is evaluated with the data available at that analysis. The relevant quantity is therefore the RCP conditional on stopping for efficacy at each analysis, which depends on the stopping boundaries, on the number of events in each region at each analysis, and on the history of the trial before that analysis. This package implements:

- **Closed-form stagewise RCPs** conditional on efficacy stopping at each analysis, and the overall RCP given success, based on the joint multivariate normal distribution of the overall and regional log hazard ratio estimates across analyses
- **Two consistency scales** for Method 1: the log hazard ratio (Ikeda and Bretz, 2010; Homma, 2024) and the hazard reduction `1 - HR` (Hayashi and Itoh, 2018), the latter by exact numerical integration or by the delta method
- **Group sequential boundaries** taken from `gsDesign` or `rpact` design objects (or supplied directly), with **binding or non-binding** futility boundaries
- **Delayed enrolment** in the region of interest through regional event proportions that change across analyses, and **region-specific hazard ratios**
- **Adjusted consistency thresholds** that restore the null RCP given efficacy stopping to its fixed-design reference value, extending the adjustment of Homma (2026)
- **Planning tools**: accrual by region with delayed regional starts under the overall accrual plan of the protocol, expected regional event proportions derived from the accrual pattern by region (with optional exponential drop-out), and the smallest proportion of events in the region of interest that achieves a target RCP (overall, stagewise, or unconditional)
- **Monte Carlo simulation** of the same normal model for verification of the numerical evaluation
- **Patient-level simulation** of event-driven trials (piecewise uniform accrual by region, exponential event times, optional exponential drop-out; log-rank test; Cox, Peto or Pike estimate of the hazard ratio, the same method in the overall population and in each region) to assess the accuracy of the normal approximation
- `print()` and `plot()` methods for all results

Regional consistency is evaluated following the Japanese Ministry of Health, Labour and Welfare (MHLW, 2007) guidance, with region 1 as the region of interest:

- **Method 1**: region 1 retains at least a fraction `gamma_M1` of the overall treatment effect
- **Method 2**: all regional treatment effects exceed the threshold `gamma_M2`

## Installation

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/rcpgsd", build_vignettes = TRUE)
```

## Vignettes

| Vignette | Contents |
|---|---|
| `vignette("rcpgsd", package = "rcpgsd")` | Model, stagewise RCPs, delayed enrolment, hazard reduction scale, planning of the proportion of Japanese events, boundaries from `gsDesign` and `rpact` |
| `vignette("thresholds-and-simulation", package = "rcpgsd")` | Null RCP inflation and adjusted thresholds, simulation of the normal model, patient-level simulation |

## Main Functions

| Function | Description |
|---|---|
| `rcpgsd()` | Stagewise and overall RCPs (Method 1 and Method 2) under a GSD, by closed form or simulation |
| `rcpgsd_get_gamma()` | Stagewise or overall consistency thresholds adjusted for the inflation of the null RCP caused by early efficacy stopping |
| `rcpgsd_sim()` | Patient-level simulation of event-driven trials estimating the same probabilities as `rcpgsd()` |
| `rcpgsd_accrual()` | Accrual by region when some regions start late, keeping the overall accrual plan (delays absorbed by the regions starting at time 0) |
| `rcpgsd_events()` | Expected analysis times and regional event proportions from the accrual pattern by region |
| `rcpgsd_get_fs()` | Smallest proportion (and number) of events in the region of interest achieving a target RCP |

## Model

Let `delta_s = -log(HR_s)` be the treatment effect in region `s`, `d_k` the cumulative number of events at analysis `k`, and `f_{s,k}` the proportion of these events in region `s`. With allocation ratio `r`, the regional information is `I_{s,k} = r f_{s,k} d_k / (1 + r)^2`. The regional log hazard ratio estimates are normally distributed with variance `1 / I_{s,k}`, have independent increments within each region, and are independent across regions. The overall estimate is the information-weighted average of the regional estimates, and the Wald statistic at analysis `k` is compared with the efficacy and futility boundaries of the design.

## Examples

### 1. Stagewise RCPs under a GSD

```r
library(rcpgsd)

# O'Brien-Fleming-type efficacy boundaries and a binding futility boundary
design <- list(upper   = c(3.7103, 2.5114, 1.9930),
               lower   = c(0, 0.5),
               binding = TRUE,
               alpha   = 0.025)

# Japan (region 1) enrols late: 6%, 9%, and 10% of the events
fs <- cbind(c(0.06, 0.47, 0.47),
            c(0.09, 0.455, 0.455),
            c(0.10, 0.45, 0.45))

res <- rcpgsd(design, events = c(150, 300, 450), fs = fs, hr = 0.72)
print(res)
plot(res)
plot(res, what = "prob")
```

### 2. Hazard reduction scale (Hayashi and Itoh, 2018)

```r
res_hr <- rcpgsd(design, events = c(150, 300, 450), fs = fs, hr = 0.72,
                 method = "M1", scale = "hr", hr_method = "exact")
print(res_hr)
```

### 3. Boundaries from gsDesign or rpact

```r
gsd <- gsDesign::gsDesign(k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
                          sfu = gsDesign::sfLDOF,
                          sfl = gsDesign::sfHSD, sflpar = -2)

# Non-binding futility: obeyed (default) or ignored
res_obey   <- rcpgsd(gsd, events = c(150, 300, 450),
                     fs = c(0.10, 0.45, 0.45), hr = 0.72)
res_ignore <- rcpgsd(gsd, events = c(150, 300, 450),
                     fs = c(0.10, 0.45, 0.45), hr = 0.72,
                     futility = "ignore")

rpd <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025,
                                       typeOfDesign = "asOF",
                                       futilityBounds = c(0, 0.5),
                                       bindingFutility = TRUE)
res_rpact <- rcpgsd(rpd, events = c(150, 300, 450),
                    fs = c(0.10, 0.45, 0.45), hr = 0.72)
```

### 4. Adjusted consistency thresholds

```r
design_B <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)

# Stagewise thresholds that restore the null RCP to the fixed-design value
adj <- rcpgsd_get_gamma(design_B, events = c(262, 523, 785),
                        fs = c(0.1, 0.45, 0.45))
print(adj)
plot(adj)

# RCPs under the alternative with the adjusted thresholds
res_adj <- rcpgsd(design_B, events = c(262, 523, 785),
                  fs = c(0.1, 0.45, 0.45), hr = exp(-0.2),
                  gamma_M1 = adj$result$gamma_M1_adj,
                  gamma_M2 = adj$result$gamma_M2_adj)
print(res_adj)
```

### 5. Fixed design and simulation check

```r
# A single-stage design gives the fixed-design results
res_fix <- rcpgsd(list(upper = qnorm(0.975)), events = 450,
                  fs = c(0.10, 0.45, 0.45), hr = 0.72)

# Monte Carlo simulation of the same model
res_sim <- rcpgsd(design, events = c(150, 300, 450), fs = fs, hr = 0.72,
                  approach = "simulation", nsim = 1e6, seed = 1)
```

### 6. Patient-level simulation of an event-driven trial

```r
design_B <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)

# 600 patients enrolled over 18 months; Japan (region 1) starts after
# month 6; median survival 12 months in the control group
n_acc <- rbind(c(0, 60),      # region 1: months 0-6, 6-18
               c(100, 170),   # region 2
               c(100, 170))   # region 3
sim <- rcpgsd_sim(design_B, events = c(150, 300, 450), hr = 0.72,
                  lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                  n_accrual = n_acc, nsim = 10000, seed = 1)
print(sim)
plot(sim)

# Closed form evaluated at the mean regional event proportions
cf <- rcpgsd(design_B, events = c(150, 300, 450), fs = sim$fs_mean,
             hr = 0.72)
cbind(simulation = sim$stage$rcp, closed_form = cf$stage$rcp)
```

### 7. Planning: expected proportions and required proportion of Japanese events

```r
# Overall plan of the protocol: 600 patients at 50 patients per month over
# 12 months. Japan (60 patients) starts 6 months late and the other
# regions absorb the delay, so the analysis times are unchanged
acc <- rcpgsd_accrual(n_region = c(60, 270, 270), start = c(6, 0, 0),
                      accrual_time = c(0, 12), rate_overall = 50)
print(acc)
ev_acc <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                        lambda_c = log(2) / 12,
                        accrual_time = acc$accrual_time,
                        n_accrual = acc$n_accrual)

# Expected regional event proportions with delayed enrolment in Japan
n_acc <- rbind(c(0, 60), c(100, 170), c(100, 170))
ev <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                    lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                    n_accrual = n_acc)
print(ev)

# Drop-out probability of 5% per year in both groups
ev_d <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
                      lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
                      n_accrual = n_acc, lambda_d = -log(1 - 0.05) / 12)

# Smallest proportion of Japanese events (keeping the delayed-enrolment
# shape) such that the overall RCP given success is at least 0.8
design_B <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
req <- rcpgsd_get_fs(design_B, events = c(150, 300, 450), fs = ev$fs,
                     hr = 0.72, target = 0.8, type = "overall")
print(req)

# Fixed design, hazard reduction scale, unconditional probability
# (reproduces Table 1 of Hayashi and Itoh, 2018)
rcpgsd_get_fs(list(upper = qnorm(0.975)), events = 161, fs = c(0.2, 0.8),
              hr = 0.6, method = "M1", type = "unconditional",
              scale = "hr")$result
```

## Verification

The test suite checks the closed-form results against:

- the RCPs reported by Homma (2026) for a single primary endpoint (fixed design);
- the consistency probability on the hazard reduction scale reported by Hayashi and Itoh (2018) (numerical integration and delta method);
- the closed-form expressions for the fixed design given in Homma (2024);
- the required numbers of Japanese events in Table 1 of Hayashi and Itoh (2018) and the constant required fraction on the log hazard ratio scale;
- numerical integration of the expected numbers of events (with and without drop-out), and the patient-level simulation for the expected regional proportions;
- independent reference values for group sequential designs computed in Python (SciPy), and Monte Carlo simulation;
- boundaries obtained from `gsDesign` and from `rpact` for the same design;
- for the patient-level simulation, `survival::survdiff()` (log-rank statistic), `survival::coxph()` and direct maximization of the partial likelihood (regional Cox estimates), and the closed form evaluated at the simulated regional event proportions.

## References

- Ministry of Health, Labour and Welfare (2007). Basic principles on global clinical trials.
- Ikeda K, Bretz F (2010). Sample size and proportion of Japanese patients in multi-regional trials. *Pharmaceutical Statistics*, 9, 207-216.
- Hayashi N, Itoh Y (2018). A re-examination of Japanese sample size calculation for multi-regional clinical trial evaluating survival endpoint. *Japanese Journal of Biometrics*, 38, 79-92.
- Homma G (2024). Cautionary note on regional consistency evaluation in multiregional clinical trials with binary outcomes. *Pharmaceutical Statistics*, 23, 385-398.
- Homma G (2026). A cautionary note on regional consistency evaluation with multiple primary endpoints in multiregional clinical trials. *Pharmaceutical Statistics*, 25, e70111.

## License

GPL (>= 2)
