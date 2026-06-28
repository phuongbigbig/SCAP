# Survival Curves Analysis Platform (SCAP)

**Version:** 1.0  
**Author:** Nguyen Thanh Phuong  
**License:** Copyright 2026

A single-file, client-side web application for survival analysis. SCAP runs entirely in the browser — no server, no dependencies, no installation. Open the HTML file, paste your data, and get publication-ready Kaplan-Meier curves, log-rank tests, Cox proportional hazards models, and forest plots.

---

## Table of Contents

1. [Features](#features)  
2. [Getting Started](#getting-started)  
3. [Data Formats](#data-formats)  
4. [Statistical Methods and Mathematics](#statistical-methods-and-mathematics)  
   - [Kaplan-Meier Estimator](#1-kaplan-meier-estimator)  
   - [Greenwood's Formula for Confidence Intervals](#2-greenwoods-formula-for-confidence-intervals)  
   - [Median Survival Time](#3-median-survival-time-and-confidence-interval)  
   - [Log-Rank Test](#4-log-rank-test-mantel-cox)  
   - [Cox Proportional Hazards Model](#5-cox-proportional-hazards-model)  
   - [Mathematical Utility Functions](#6-mathematical-utility-functions)  
5. [Curve Types](#curve-types)  
6. [Visualizations](#visualizations)  
7. [Export Options](#export-options)  
8. [Customization](#customization)  
9. [Technical Notes](#technical-notes)  

---

## Features

- Upload or paste survival data (tab- or comma-delimited)
- Two input formats: aggregate (tank-level death counts) and individual (one row per subject)
- Automatic group detection and reference group selection
- Kaplan-Meier survival curves, cumulative mortality, and cumulative hazard plots
- 95% confidence intervals via Greenwood's formula
- Log-rank (Mantel-Cox) test for K-group comparison with proper variance-covariance matrix inversion
- Cox proportional hazards regression with hazard ratios, Wald confidence intervals, and p-values
- Forest plot of hazard ratios on a log2 scale
- Number-at-risk table beneath the KM plot
- SVG and PNG export at 3x resolution
- CSV and XLSX export of statistical results
- Seven page color themes and seven journal color palettes (NPG, AAAS, NEJM, Lancet, JCO, Classic, D3)
- Fully client-side: zero network requests, zero dependencies

---

## Getting Started

1. Open `survival_analysis_platform.html` in any modern browser.
2. Click **Load Example (Aggregate)** or **Load Example (Individual)** to see sample data.
3. Or paste your own data into the text area and click **Run Analysis**.
4. Select a reference (control) group from the dropdown in the Groups card.
5. Adjust appearance settings (line width, font scale, plot width, aspect ratio, curve type, color palette).
6. Export plots and statistics using the Export card.

---

## Data Formats

### Aggregate Format

Each row represents one time point for one tank/replicate, with a count of deaths at that time. Required columns (case-insensitive, flexible naming):

| Column | Aliases | Description |
|--------|---------|-------------|
| Group | `exp`, `group`, `condition`, `treatment` | Treatment group name |
| Tank | `tank`, `replicate`, `tank_id`, `rep` | Replicate identifier (optional) |
| Day | `day`, `time`, `days`, `timepoint` | Time point |
| Start N | `start_number`, `n`, `n_start`, `start_n`, `total` | Initial number of subjects in this tank |
| Deaths | `daily_dead`, `deaths`, `dead`, `n_dead`, `events` | Number of deaths at this time point |

The platform automatically expands aggregate data to individual-level records: each death becomes one row with `status=1` at the recorded day, and surviving subjects become `status=0` (censored) at the maximum observed day for that tank.

### Individual Format

One row per subject. Required columns:

| Column | Aliases | Description |
|--------|---------|-------------|
| Group | `group`, `exp`, `condition`, `treatment` | Treatment group name |
| Time | `time`, `day`, `days`, `survival_time` | Survival or censoring time |
| Status | `status`, `event`, `dead`, `censored` | 1 = event (death), 0 = censored |
| Tank | `tank`, `replicate`, `tank_id` | Replicate identifier (optional) |

---

## Statistical Methods and Mathematics

### 1. Kaplan-Meier Estimator

The Kaplan-Meier (KM) estimator is a non-parametric method for estimating the survival function from right-censored data (Kaplan & Meier, 1958).

**Definition.** Let t₁ < t₂ < ... < tₘ be the ordered, distinct event (death) times. At each event time tⱼ, let:

- nⱼ = number of subjects at risk just before time tⱼ (alive and not yet censored)
- dⱼ = number of events (deaths) at time tⱼ

The KM estimate of the survival function is:

```
S(t) = ∏ (1 - dⱼ / nⱼ)    for all j where tⱼ ≤ t
```

That is, S(t) is the product of the conditional survival probabilities at each event time up to and including time t.

**Handling of censored observations.** Subjects censored between two consecutive event times tⱼ and tⱼ₊₁ are included in the risk set at tⱼ but removed before tⱼ₊₁. This is computed by counting censored subjects with observation times strictly between the previous event time and the current event time.

**Initial conditions.** S(0) = 1 (all subjects alive at time zero). The risk set at the first event time includes all n subjects minus any censored before that time.

### 2. Greenwood's Formula for Confidence Intervals

The variance of the KM estimator is calculated using Greenwood's formula (Greenwood, 1926), which gives a pointwise estimate based on the cumulative hazard contributions at each event time.

**Variance.** The estimated variance of S(t) is:

```
Var[S(t)] = S(t)² × Σ dⱼ / [nⱼ × (nⱼ - dⱼ)]    summed over all tⱼ ≤ t
```

The standard error is:

```
SE[S(t)] = S(t) × √( Σ dⱼ / [nⱼ × (nⱼ - dⱼ)] )
```

**95% Confidence Interval.** The pointwise confidence interval uses the normal approximation:

```
S(t) ± 1.96 × SE[S(t)]
```

clamped to the interval [0, 1].

**Note:** This is the "plain" or "linear" CI. Alternative transformations (log, log-log, arcsine-square-root) exist but are not used here. The linear CI can occasionally produce bounds outside [0, 1], which is handled by clamping.

### 3. Median Survival Time and Confidence Interval

**Median survival** is the smallest time t at which S(t) ≤ 0.5:

```
median = min{ t : S(t) ≤ 0.5 }
```

If S(t) never drops to 0.5, the median is reported as "NR" (Not Reached).

**Confidence interval for the median** uses the Brookmeyer-Crowley method (Brookmeyer & Crowley, 1982):

- **Lower bound:** the smallest time where the *upper* CI of S(t) drops to or below 0.5
- **Upper bound:** the smallest time where the *lower* CI of S(t) drops to or below 0.5

Intuitively, if even the optimistic (upper) survival bound has dropped below 50%, we are confident the true median is at least that early; if even the pessimistic (lower) bound has dropped below 50%, the true median is at most that late.

### 4. Log-Rank Test (Mantel-Cox)

The log-rank test is a non-parametric hypothesis test for comparing survival distributions across K ≥ 2 groups (Mantel, 1966; Cox, 1972). The null hypothesis is that all groups have identical survival functions.

**Setup.** At each distinct event time tⱼ across all groups, define:

- N = total number at risk across all groups
- D = total number of events across all groups
- nⱼ⁽ᵍ⁾ = number at risk in group g
- dⱼ⁽ᵍ⁾ = number of events in group g

**Expected events.** Under the null hypothesis (equal survival), the expected number of events in group g at time tⱼ is:

```
eⱼ⁽ᵍ⁾ = nⱼ⁽ᵍ⁾ × D / N
```

This is the hypergeometric expectation: events are distributed proportionally to group size in the risk set.

**Observed minus Expected.** For each group g, accumulate over all event times:

```
(O - E)ᵍ = Σⱼ [ dⱼ⁽ᵍ⁾ - eⱼ⁽ᵍ⁾ ]
```

**Variance-Covariance Matrix.** The variance-covariance of the (O-E) vector comes from the hypergeometric distribution:

```
Vᵍʰ = Σⱼ  D × (N - D) × nⱼ⁽ᵍ⁾ × (δᵍʰ × N - nⱼ⁽ʰ⁾) / [N² × (N - 1)]
```

where δᵍʰ is the Kronecker delta (1 if g = h, 0 otherwise).

For the **diagonal** entries (g = h), this simplifies to:

```
Vᵍᵍ = Σⱼ  D × (N - D) × nⱼ⁽ᵍ⁾ × (N - nⱼ⁽ᵍ⁾) / [N² × (N - 1)]
```

For the **off-diagonal** entries (g ≠ h):

```
Vᵍʰ = Σⱼ  -D × (N - D) × nⱼ⁽ᵍ⁾ × nⱼ⁽ʰ⁾ / [N² × (N - 1)]
```

**Test Statistic.** Since the K group (O-E) values sum to zero (Σ (O-E)ᵍ = 0), only K-1 are independent. The chi-squared test statistic is:

```
χ² = (O-E)ᵀ V⁻¹ (O-E)
```

where (O-E) is the vector of the first K-1 group values and V is the corresponding (K-1) × (K-1) submatrix.

For **K = 2**, this reduces to:

```
χ² = [(O-E)₁]² / V₁₁
```

For **K > 2**, the (K-1) × (K-1) matrix V is inverted using Gauss-Jordan elimination with partial pivoting, and the quadratic form is computed directly.

**Degrees of freedom:** K - 1.

**P-value:** from the chi-squared distribution with K-1 degrees of freedom (see [Section 6](#6-mathematical-utility-functions)).

### 5. Cox Proportional Hazards Model

The Cox PH model (Cox, 1972) is a semi-parametric regression model that estimates the effect of covariates on the hazard (instantaneous risk of event). SCAP uses it with a single binary covariate to compare each non-reference group against the reference group.

**Model.** The hazard for subject i at time t is:

```
h(t | xᵢ) = h₀(t) × exp(β × xᵢ)
```

where h₀(t) is an unspecified baseline hazard, xᵢ ∈ {0, 1} is the group indicator (0 = reference, 1 = treatment), and β is the log hazard ratio.

**Partial Likelihood.** Cox's key insight is that β can be estimated without specifying h₀(t). The partial likelihood is:

```
L(β) = ∏ᵢ  [ exp(β × xᵢ) / Σⱼ∈R(tᵢ) exp(β × xⱼ) ]
```

where the product is over all event times, and R(tᵢ) is the risk set at time tᵢ (all subjects still alive and uncensored just before tᵢ).

**Log Partial Likelihood.** Taking the log:

```
ℓ(β) = Σᵢ [ β × xᵢ - log( Σⱼ∈R(tᵢ) exp(β × xⱼ) ) ]
```

**Newton-Raphson Optimization.** The maximum partial likelihood estimator β̂ is found iteratively:

```
β_{k+1} = β_k + U(β_k) / I(β_k)
```

where U is the score function (first derivative of ℓ) and I is the observed information (negative second derivative of ℓ).

**Score Function.** At each event time tᵢ, define the weighted risk set sums:

```
S₀(β, tᵢ) = Σⱼ∈R(tᵢ) exp(β × xⱼ)
S₁(β, tᵢ) = Σⱼ∈R(tᵢ) xⱼ × exp(β × xⱼ)
S₂(β, tᵢ) = Σⱼ∈R(tᵢ) xⱼ² × exp(β × xⱼ)
```

Then:

```
U(β) = Σᵢ [ xᵢ - S₁(β, tᵢ) / S₀(β, tᵢ) ]
```

**Observed Information:**

```
I(β) = Σᵢ [ S₂(β, tᵢ) / S₀(β, tᵢ) - (S₁(β, tᵢ) / S₀(β, tᵢ))² ]
```

**Tie Handling (Breslow Method).** When multiple events occur at the same time, the Breslow approximation (Breslow, 1974) treats each event as if it has access to the full risk set at that time. This is the default in most software and is implemented by computing cumulative sums from the tail of the time-sorted data.

**Implementation Detail.** The risk set sums S₀, S₁, S₂ are computed efficiently using reverse cumulative sums: sort subjects by time, then accumulate exp(β × xⱼ) from the last subject backward. This avoids an O(n²) loop over risk sets.

**Convergence.** The algorithm iterates until |Δβ| < 10⁻⁹ or 50 iterations.

**Outputs:**

| Quantity | Formula |
|----------|---------|
| Hazard Ratio (HR) | HR = exp(β̂) |
| Standard Error | SE = 1 / √I(β̂) |
| 95% Confidence Interval | exp(β̂ ± 1.96 × SE) |
| Wald z-statistic | z = \|β̂\| / SE |
| P-value | p = 2 × [1 - Φ(z)] |

**Interpretation.** HR > 1 means the treatment group has a higher hazard (dies faster) compared to the reference. HR < 1 means the treatment group has a lower hazard (lives longer). HR = 1 means no difference.

### 6. Mathematical Utility Functions

The statistical tests require computing p-values from the chi-squared and standard normal distributions. Since the platform is entirely client-side JavaScript with no libraries, these are implemented from scratch using classical numerical methods.

#### Log-Gamma Function (Lanczos Approximation)

The gamma function Γ(x) is needed for the incomplete gamma function used in chi-squared p-values. The implementation uses the Lanczos approximation (Lanczos, 1964) with 6 terms:

```
log Γ(x) ≈ -tmp + log(√(2π) × ser / x)
```

where tmp = (x + 5.5) - (x + 0.5) × log(x + 5.5), and ser is a polynomial sum over the Lanczos coefficients. This provides ~15 digits of accuracy for x > 0.

#### Error Function (Abramowitz & Stegun Approximation)

The error function erf(x) is approximated using the rational approximation from Abramowitz & Stegun (1964), equation 7.1.26:

```
erf(x) ≈ 1 - (a₁t + a₂t² + a₃t³ + a₄t⁴ + a₅t⁵) × exp(-x²)
```

where t = 1 / (1 + 0.3275911 × |x|). Maximum error: 1.5 × 10⁻⁷.

#### Standard Normal CDF (Φ)

The cumulative distribution function of the standard normal distribution:

```
Φ(x) = ½ × [1 + erf(x / √2)]
```

Used for computing p-values from the Wald z-statistic in the Cox PH model.

#### Regularized Incomplete Gamma Functions

The chi-squared p-value requires the regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x), where:

```
p-value = P(χ² > observed | df) = Q(df/2, observed/2)
```

Two complementary algorithms are used depending on the regime:

**P(a, x) — Series expansion** (for x < a + 1):

```
P(a, x) = [exp(-x) × x^a / Γ(a)] × Σₙ x^n / [(a)(a+1)...(a+n)]
```

Converges rapidly when x is small relative to a.

**Q(a, x) — Continued fraction** (for x ≥ a + 1), using the Lentz-Thompson-Barnett algorithm:

```
Q(a, x) = [exp(-x) × x^a / Γ(a)] × CF(a, x)
```

where CF is the modified Lentz continued fraction. Converges rapidly when x is large relative to a.

Both use convergence tolerance of 10⁻¹² to 10⁻¹⁴.

#### Chi-Squared P-value

```
p = Q(df/2, χ²/2)
```

This maps the chi-squared test statistic with a given degrees of freedom to a tail probability.

---

## Curve Types

SCAP supports three curve types, selectable via the segmented control:

| Curve | Y-axis | Formula |
|-------|--------|---------|
| **Survival** | Survival probability | S(t) |
| **Mortality** | Cumulative mortality | 1 - S(t) |
| **Cumulative Hazard** | Nelson-Aalen cumulative hazard | H(t) = -log S(t) |

The Nelson-Aalen estimator H(t) = -log S(t) is derived from the relationship between the survival function and the cumulative hazard. Under the proportional hazards assumption, log H(t) plots for different groups should be parallel.

---

## Visualizations

### Kaplan-Meier Plot

- Step function for each group, colored by the selected journal palette
- Optional 95% CI shading (Greenwood)
- Optional censoring tick marks at censored observation times
- Optional median survival dashed lines (horizontal at S=0.5, vertical at median time)
- Optional log-rank p-value annotation (bottom-left corner)
- Number-at-risk table below the x-axis

### Forest Plot

- One row per non-reference group
- Point estimate (filled circle) at the hazard ratio
- Optional 95% CI bar (Wald interval)
- Log2-scale x-axis with dashed reference line at HR = 1
- Annotation text showing HR [lower-upper] and p-value
- Dynamic x-axis range that always includes HR = 1

---

## Export Options

| Button | Output |
|--------|--------|
| KM · SVG | Vector graphics of the Kaplan-Meier plot |
| KM · PNG | Raster image at 3× resolution |
| Forest · SVG | Vector graphics of the forest plot |
| Forest · PNG | Raster image at 3× resolution |
| Stats · CSV | Comma-separated statistical results |
| Stats · XLSX | Excel workbook (minimal zipStore format, no external dependencies) |

The XLSX export constructs a valid Open XML Spreadsheet package entirely in JavaScript using a minimal ZIP (stored, no compression) builder with CRC-32 checksums.

---

## Customization

### Appearance Controls

- **Line width:** Stroke thickness multiplier for plot lines (default 1)
- **Font scale:** Multiplier for all text in SVG plots (default 1)
- **Plot width:** Width in pixels (default 640)
- **Aspect ratio:** 7:5, 4:3, or 16:10 (default 16:10)
- **Show CI:** Toggle confidence interval shading on the KM plot
- **Show Forest CI:** Toggle confidence interval bars on the forest plot (independent of KM CI)
- **Show censor marks:** Toggle tick marks at censored times
- **Show median lines:** Toggle dashed median survival lines
- **Show p-value:** Toggle log-rank p-value annotation

### Color Palettes

Seven journal-inspired color palettes are available:

| Palette | Source |
|---------|--------|
| NPG | Nature Publishing Group |
| AAAS | American Association for the Advancement of Science |
| NEJM | New England Journal of Medicine |
| Lancet | The Lancet |
| JCO | Journal of Clinical Oncology |
| Classic | Seaborn classic |
| D3 | D3.js category10 |

### Page Themes

Seven page themes control the header, background, and accent colors: Blue, Teal, Violet, Crimson, Green, Amber, and Graphite. The default is Teal.

---

## Technical Notes

### Architecture

The entire application is a single self-contained HTML file (~830 lines) with inline CSS and JavaScript. There are no external dependencies, no build step, and no network requests. All statistical computations run in the browser's JavaScript engine.

### Numerical Precision

- The Lanczos approximation for log Γ(x) provides ~15 digits of accuracy.
- The error function approximation has a maximum error of 1.5 × 10⁻⁷.
- The incomplete gamma function uses series and continued fraction representations with convergence tolerances of 10⁻¹² to 10⁻¹⁴.
- Newton-Raphson for Cox PH converges to |Δβ| < 10⁻⁹.
- Chi-squared p-values have been verified against standard tables (χ²=3.84, df=1 → p=0.0500; χ²=5.99, df=2 → p=0.0500; χ²=7.81, df=3 → p=0.0501).

### Limitations

- **No pseudoreplication correction.** The platform treats each expanded individual as independent. For clustered designs (multiple fish per tank), use R with `survival::coxph(..., cluster=tank)` or frailty models. The companion R Markdown files (`KM_Survival_DHC_vs_DC.Rmd`) implement these corrections.
- **Breslow tie handling only.** The Efron approximation (slightly more accurate for heavy ties) is not implemented.
- **Single covariate.** Cox PH is fit with one binary group indicator. Multivariable models require R or other statistical software.
- **Linear confidence intervals.** Greenwood CIs use the normal approximation on the survival scale, which can produce bounds outside [0, 1] (clamped). Log-log transformed CIs are not implemented.
- **Complete separation.** If all events in one group occur before any events in another, the Cox model will not converge (beta → ±∞). This is a known property of maximum likelihood estimation with complete separation.

### Browser Compatibility

Tested in modern browsers (Chrome, Firefox, Safari, Edge). Requires ES5+ JavaScript support. The application uses `Blob`, `URL.createObjectURL`, `TextEncoder`, and `XMLSerializer` for export functionality.

---

## References

- Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*, 53(282), 457-481.
- Greenwood, M. (1926). *Reports on Public Health and Medical Subjects*, No. 33. London: HMSO.
- Mantel, N. (1966). Evaluation of survival data and two new rank order statistics arising in its consideration. *Cancer Chemotherapy Reports*, 50(3), 163-170.
- Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*, 34(2), 187-220.
- Breslow, N. E. (1974). Covariance analysis of censored survival data. *Biometrics*, 30(1), 89-99.
- Brookmeyer, R., & Crowley, J. (1982). A confidence interval for the median survival time. *Biometrics*, 38(1), 29-41.
- Lanczos, C. (1964). A precision approximation of the gamma function. *SIAM Journal on Numerical Analysis*, 1(1), 86-96.
- Abramowitz, M., & Stegun, I. A. (1964). *Handbook of Mathematical Functions*. National Bureau of Standards.
