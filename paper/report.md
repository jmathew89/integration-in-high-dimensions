# The Cost of Dimension
## A Computational Investigation of Deterministic and Monte Carlo Integration

**Author:** Jonah Mathew  
**Repository:** https://github.com/jmathew89/integration-in-high-dimensions

## Abstract

Numerical integration becomes increasingly difficult as the dimension of the integration domain grows. Tensor-product deterministic rules can be highly accurate in low dimensions, but their number of function evaluations grows rapidly with dimension. Monte Carlo integration has a slower convergence rate in low-dimensional settings, but its characteristic \(N^{-1/2}\) rate does not deteriorate directly with dimension. This investigation compares the midpoint rule and Monte Carlo integration under fixed maximum evaluation budgets.

The experiments use two functions over the unit hypercube, a sum of squares and a product of squares. Both have known exact integrals, allowing absolute error, relative error, and root mean square error (RMSE) to be evaluated. The study includes a one-dimensional baseline, convergence and variability experiments, comparisons across dimensions one through ten, and a separate runtime experiment.

The results show that the midpoint rule is substantially more accurate for the one-dimensional baseline. In three dimensions, observed convergence rates closely match theoretical predictions. In the multidimensional experiments at a budget of 10,000 evaluations, the observed crossover occurs between dimensions four and five for the sum-of-squares function, and between dimensions five and six for the product-of-squares function. These different crossover points demonstrate that dimension alone does not determine which method performs better. The integrand and available evaluation budget also matter.

---

## 1. Introduction

### 1.1 Background and Motivation

Definite integrals are used throughout mathematics, physics, statistics, and engineering. Many integrals can be evaluated exactly using symbolic methods, but practical problems often require numerical approximations. Numerical integration methods estimate an integral using a finite number of evaluations of the integrand.

The difficulty of numerical integration can increase substantially when the domain has many dimensions. A deterministic method such as the midpoint rule divides the integration domain into smaller regions and evaluates the function at selected points. In multiple dimensions, a straightforward extension uses a grid of points across every coordinate. If the grid has \(m\) subdivisions in each dimension, it requires \(m^d\) evaluations in dimension \(d\). Consequently, the number of evaluations grows exponentially with dimension when the resolution per coordinate is held fixed.

Monte Carlo integration takes a different approach. Rather than constructing a regular grid, it estimates the integral using randomly sampled points. Under appropriate assumptions, its root mean square error decreases proportionally to \(N^{-1/2}\), where \(N\) is the number of samples. This convergence rate does not contain the dimension in its exponent. However, the error also depends on the variance of the function being integrated, which can itself change significantly with dimension.

### 1.2 Research Question

**How does the dimension of an integration problem affect the relative accuracy and computational efficiency of midpoint-rule integration and Monte Carlo integration under a fixed evaluation budget?**

### 1.3 Scope and Approach

This investigation compares the composite midpoint rule with standard Monte Carlo integration over the unit hypercube \([0,1]^d\).

Two test functions are used:

1. The sum of squares, \(f_d(\mathbf{x})=\sum_{i=1}^{d}x_i^2\).
2. The product of squares, \(g_d(\mathbf{x})=\prod_{i=1}^{d}x_i^2\).

Both functions have known exact integrals. The experiments examine dimensions one through ten and nominal evaluation budgets of 100, 1,000, and 10,000. Additional experiments investigate one-dimensional accuracy, convergence rates, Monte Carlo variability, and runtime.

The midpoint method uses the largest complete tensor-product grid that does not exceed the specified evaluation budget. Monte Carlo uses the full budget for each trial. The methods are therefore compared under the same maximum budget, but their actual numbers of evaluations are not always identical.

### 1.4 Contribution and Limitations

The investigation combines theoretical derivations with computational experiments. Its purpose is not to introduce a new integration method, but to examine how the theoretical properties of two established methods appear in controlled examples.

The conclusions are limited to the two selected integrands, the midpoint rule implementation, the tested dimensions, and the tested budgets. The results should not be interpreted as a universal ranking of deterministic and stochastic integration methods.

---

## 2. Mathematical Foundation

### 2.1 The Multidimensional Integral

Let the integration domain be the \(d\)-dimensional unit hypercube:

\[
[0,1]^d.
\]

For a function \(f:[0,1]^d\rightarrow\mathbb{R}\), define

\[
I_d=\int_{[0,1]^d}f(\mathbf{x})\,d\mathbf{x},
\qquad
\mathbf{x}=(x_1,x_2,\ldots,x_d).
\]

The volume of the unit hypercube is one. Therefore, if \(\mathbf{X}\) is uniformly distributed over this domain, the integral can also be expressed as an expected value:

\[
I_d=\mathbb{E}[f(\mathbf{X})].
\]

This identity provides the mathematical basis for Monte Carlo integration.

### 2.2 The Midpoint Rule

In one dimension, divide the interval \([0,1]\) into \(m\) equal subintervals. The midpoint of subinterval \(j\), where \(j=0,1,\ldots,m-1\), is

\[
x_j=\frac{j+\frac12}{m}.
\]

The composite midpoint approximation is

\[
M_m=\frac{1}{m}\sum_{j=0}^{m-1}
f\left(\frac{j+\frac12}{m}\right).
\]

In dimension \(d\), the tensor-product midpoint rule uses a grid with \(m\) subdivisions along every coordinate:

\[
M_{d,m}=\frac{1}{m^d}
\sum_{j_1=0}^{m-1}\cdots\sum_{j_d=0}^{m-1}
f\left(
\frac{j_1+\frac12}{m},\ldots,
\frac{j_d+\frac12}{m}
\right).
\]

The total number of function evaluations is

\[
N=m^d.
\]

This exponential dependence on dimension is the central computational limitation of the tensor-product midpoint rule.

For a maximum evaluation budget \(B\), the largest permitted grid uses

\[
m=\left\lfloor B^{1/d}\right\rfloor,
\qquad N=m^d\leq B.
\]

Because \(m\) must be an integer, the actual midpoint evaluation count can be substantially smaller than the nominal budget, particularly in higher dimensions.

### 2.3 Monte Carlo Integration

Let \(\mathbf{X}_1,\ldots,\mathbf{X}_N\) be independent random vectors, each uniformly distributed over \([0,1]^d\). The Monte Carlo estimator is

\[
\widehat I_N=
\frac{1}{N}\sum_{i=1}^{N}f(\mathbf{X}_i).
\]

Since each sample has expected value \(I_d\), the estimator is unbiased:

\[
\mathbb{E}[\widehat I_N]=I_d.
\]

If the variance of \(f(\mathbf{X})\) is \(\sigma_f^2\), independence gives

\[
\operatorname{Var}(\widehat I_N)=\frac{\sigma_f^2}{N}.
\]

The standard deviation of the estimator is therefore

\[
\operatorname{SD}(\widehat I_N)=\frac{\sigma_f}{\sqrt{N}}.
\]

This is the source of the characteristic \(N^{-1/2}\) Monte Carlo convergence rate. Increasing the sample size by a factor of 100 reduces the theoretical standard deviation by a factor of 10, assuming the integrand and sampling distribution remain unchanged.

For an unbiased estimator, the theoretical standard deviation is also its root mean square error. In finite experiments, RMSE can be estimated from repeated trials by

\[
\operatorname{RMSE}=
\sqrt{\frac{1}{T}\sum_{t=1}^{T}
(\widehat I_N^{(t)}-I_d)^2},
\]

where \(T\) is the number of independent trials.

### 2.4 Test Function A: Sum of Squares

The first test function is

\[
f_d(\mathbf{x})=\sum_{i=1}^{d}x_i^2.
\]

Its exact integral is

\[
I_d=\sum_{i=1}^{d}\int_0^1x_i^2\,dx_i
=\frac{d}{3}.
\]

For the one-dimensional midpoint rule applied to \(x^2\),

\[
M_m=\frac{1}{m}\sum_{j=0}^{m-1}
\left(\frac{j+\frac12}{m}\right)^2
=\frac13-\frac{1}{12m^2}.
\]

By additivity, the multidimensional midpoint approximation is

\[
M_{d,m}=d\left(\frac13-\frac{1}{12m^2}\right).
\]

Thus, its absolute error is

\[
|M_{d,m}-I_d|=\frac{d}{12m^2},
\]

and its relative error is

\[
\frac{|M_{d,m}-I_d|}{I_d}
=\frac{1}{4m^2}.
\]

For a uniform random variable \(X\) on \([0,1]\),

\[
\mathbb{E}[X^2]=\frac13,
\qquad
\mathbb{E}[X^4]=\frac15,
\qquad
\operatorname{Var}(X^2)=\frac4{45}.
\]

The sum of \(d\) independent squared coordinates therefore has variance \(4d/45\). The Monte Carlo estimator has variance

\[
\operatorname{Var}(\widehat I_N)=\frac{4d}{45N},
\]

and theoretical RMSE

\[
\operatorname{RMSE}=\sqrt{\frac{4d}{45N}}.
\]

Dividing by the exact integral gives the relative RMSE:

\[
\frac{\operatorname{RMSE}}{I_d}
=\sqrt{\frac{4}{5dN}}.
\]

For this particular function, the relative Monte Carlo RMSE decreases with dimension, even though the absolute RMSE increases.

### 2.5 Test Function B: Product of Squares

The second test function is

\[
g_d(\mathbf{x})=\prod_{i=1}^{d}x_i^2.
\]

Its exact integral is

\[
J_d=\prod_{i=1}^{d}\int_0^1x_i^2\,dx_i
=\left(\frac13\right)^d.
\]

Because the function is a product of independent coordinate terms, its midpoint approximation is

\[
G_{d,m}=
\left(\frac13-\frac{1}{12m^2}\right)^d
=\left(\frac13\right)^d
\left(1-\frac{1}{4m^2}\right)^d.
\]

Its relative absolute error is therefore

\[
\frac{|G_{d,m}-J_d|}{J_d}
=1-\left(1-\frac{1}{4m^2}\right)^d.
\]

For a uniformly distributed coordinate \(X\),

\[
\mathbb{E}[X^2]=\frac13,
\qquad
\mathbb{E}[X^4]=\frac15.
\]

By independence,

\[
\mathbb{E}[g_d(\mathbf{X})^2]=\left(\frac15\right)^d.
\]

Consequently,

\[
\operatorname{Var}(g_d(\mathbf{X}))
=\left(\frac15\right)^d-\left(\frac19\right)^d,
\]

and the Monte Carlo estimator has RMSE

\[
\operatorname{RMSE}
=\sqrt{\frac{(1/5)^d-(1/9)^d}{N}}.
\]

Its relative RMSE is

\[
\frac{\operatorname{RMSE}}{J_d}
=\sqrt{\frac{(9/5)^d-1}{N}}.
\]

Unlike the sum-of-squares function, the relative Monte Carlo RMSE for the product-of-squares function increases with dimension. The exact integral decreases exponentially, which makes relative error especially important when comparing the methods.

### 2.6 Theoretical Convergence Comparison

For the sum-of-squares function, the midpoint relative error can be expressed in terms of its actual evaluation count \(N=m^d\):

\[
E_{\mathrm{mid}}=\frac{1}{4N^{2/d}}.
\]

Its error exponent is therefore \(-2/d\). Monte Carlo retains its \(-1/2\) exponent. As dimension increases, the midpoint convergence exponent becomes less favorable.

For the product-of-squares function, the midpoint relative error is

\[
E_{\mathrm{mid}}=
1-\left(1-\frac{1}{4N^{2/d}}\right)^d,
\]

while the Monte Carlo relative RMSE is

\[
E_{\mathrm{MC}}=
\sqrt{\frac{(9/5)^d-1}{N}}.
\]

These expressions show why the relative performance of the two methods can differ between integrands. A convergence rate alone does not determine the error at a particular finite evaluation budget.

---

## 3. Experimental Methodology

### 3.1 Computational Environment

The experiments were implemented in Python using NumPy for numerical calculations and Matplotlib for plotting. Random-number generators were initialized with fixed seeds to make the stochastic experiments reproducible.

The experiment scripts are separated by purpose so that the baseline, convergence, variability, multidimensional, and runtime analyses can be examined independently.

### 3.2 Evaluation Budgets

The multidimensional comparison uses nominal budgets of 100, 1,000, and 10,000 evaluations.

For each dimension \(d\), the midpoint method selects the largest integer \(m\) satisfying

\[
m^d\leq B.
\]

It then performs \(m^d\) evaluations. Monte Carlo performs \(B\) evaluations per trial. This design ensures that neither method exceeds the nominal budget, but it does not guarantee equal evaluation counts.

This distinction is important when interpreting the results. In higher dimensions, the midpoint grid may use only a fraction of the available budget because a complete tensor-product grid is required.

### 3.3 Error Measurements

The experiments use several measures:

- **Absolute error:** the absolute difference between an approximation and the exact integral.
- **Relative error:** absolute error divided by the absolute value of the exact integral.
- **RMSE:** the square root of the mean squared error across Monte Carlo trials.
- **Standard deviation:** the spread of Monte Carlo estimates across independent trials.
- **Runtime:** elapsed computation time measured separately from the accuracy experiments.

For the multidimensional comparisons, midpoint absolute error is compared with Monte Carlo RMSE. Relative errors are also reported because the exact integral changes substantially between the two test functions.

### 3.4 Monte Carlo Trials

The multidimensional comparison uses 300 independent Monte Carlo trials per configuration. The variability experiment uses 500 trials for each sample size. The one-dimensional baseline uses 200 trials for each sample size.

The percentile range reported in the multidimensional experiment describes the central 95% of the observed trial estimates. It is an empirical distribution range, not a confidence interval for the exact integral.

### 3.5 Runtime Measurements

Runtime is measured in a separate experiment using repeated function evaluations. The reported values are median runtimes from nine timing repetitions, with a batch of evaluations used to reduce timing noise. A warm-up evaluation is performed before measurements.

Runtime results are implementation- and hardware-dependent. They are treated as secondary evidence rather than as universal comparisons of the mathematical methods.

---

## 4. One-Dimensional Baseline Experiment

### 4.1 Objective

The baseline experiment compares the two methods on

\[
f(x)=x^2,
\qquad
I=\int_0^1x^2\,dx=\frac13.
\]

This simple problem provides a reference point before considering the effect of dimension.

### 4.2 Theoretical Predictions

For the midpoint rule with \(N\) subdivisions,

\[
M_N=\frac13-\frac{1}{12N^2},
\]

so the absolute error is

\[
E_{\mathrm{mid}}=\frac{1}{12N^2}.
\]

The midpoint error therefore decreases proportionally to \(N^{-2}\).

For Monte Carlo integration, the variance of \(X^2\) is \(4/45\). Its theoretical RMSE is

\[
E_{\mathrm{MC}}=\sqrt{\frac{4}{45N}},
\]

which decreases proportionally to \(N^{-1/2}\).

### 4.3 Results

| Evaluations \(N\) | Midpoint absolute error | Monte Carlo RMSE |
|---:|---:|---:|
| 10 | 0.00083333333 | 0.095659271 |
| 30 | 0.000092592593 | 0.054478039 |
| 100 | 0.0000083333333 | 0.029704355 |
| 300 | 0.00000092592593 | 0.017166070 |
| 1,000 | 0.000000083333333 | 0.0090261356 |
| 3,000 | 0.0000000092592593 | 0.005244207 |
| 10,000 | 0.00000000083333324 | 0.0030899155 |

The observed log-log slopes were \(-2.0000\) for the midpoint rule and \(-0.5016\) for Monte Carlo. These closely match the theoretical slopes of \(-2\) and \(-0.5\), respectively.

### 4.4 Interpretation

The midpoint rule is substantially more accurate than Monte Carlo integration for this one-dimensional function at the tested sample sizes. This is consistent with its faster theoretical convergence rate for a smooth function.

The baseline illustrates why Monte Carlo integration should not be assumed to be more accurate simply because it has a dimension-independent convergence exponent. In low dimensions, a deterministic method can take advantage of the structure and smoothness of the integrand.

### 4.5 Limitations

This experiment uses a single smooth function in one dimension. It does not test the growth in grid size that motivates the high-dimensional comparison. Its purpose is to establish a baseline, not to determine which method is preferable for all one-dimensional integrals.

---

## 5. Convergence and Variability Experiments

### 5.1 Three-Dimensional Convergence

The convergence experiment uses the sum-of-squares function in three dimensions:

\[
f_3(x,y,z)=x^2+y^2+z^2,
\qquad I_3=1.
\]

The theoretical midpoint relative error is

\[
E_{\mathrm{mid}}=\frac{1}{4N^{2/3}},
\]

and the theoretical Monte Carlo RMSE is

\[
E_{\mathrm{MC}}=\sqrt{\frac{4}{15N}}.
\]

The observed results are:

| Evaluations \(N\) | Midpoint error | Monte Carlo RMSE |
|---:|---:|---:|
| 27 | 0.0277778 | 0.0950469 |
| 64 | 0.0156250 | 0.0655715 |
| 125 | 0.0100000 | 0.0453954 |
| 216 | 0.00694444 | 0.0350659 |
| 343 | 0.00510204 | 0.0260757 |
| 512 | 0.00390625 | 0.0236300 |
| 729 | 0.00308642 | 0.0193721 |
| 1,000 | 0.00250000 | 0.0165300 |
| 1,728 | 0.00173611 | 0.0118083 |
| 4,096 | 0.000976562 | 0.00792917 |
| 8,000 | 0.000625000 | 0.00593936 |

The observed midpoint log-log slope was \(-0.6667\), matching the theoretical value \(-2/3\). The observed Monte Carlo slope was \(-0.4953\), close to the theoretical value \(-0.5\).

The results support the predicted convergence behavior in three dimensions. The midpoint rule remains more accurate over the displayed sample sizes, but its convergence exponent is less favorable than in the one-dimensional baseline.

### 5.2 Monte Carlo Variability

A separate experiment investigates how the spread of Monte Carlo estimates changes with sample size. It uses the three-dimensional sum-of-squares function and 500 independent trials per sample size.

The theoretical standard deviation is

\[
\sigma_N=\sqrt{\frac{4d}{45N}}
=\sqrt{\frac{4}{15N}},
\qquad d=3.
\]

| Sample size \(N\) | Observed standard deviation | Theoretical standard deviation |
|---:|---:|---:|
| 100 | 0.049166 | 0.051640 |
| 300 | 0.029579 | 0.029814 |
| 1,000 | 0.017194 | 0.016330 |
| 3,000 | 0.009655 | 0.009428 |
| 10,000 | 0.005359 | 0.005164 |

The observed log-log slope was \(-0.4822\), compared with the theoretical value of \(-0.5\). Increasing the sample size from 100 to 10,000 reduced the observed standard deviation by a factor of approximately 9.2, close to the theoretical factor of 10.

### 5.3 Interpretation

The convergence experiment confirms that the observed error trends closely follow the theoretical rates. The variability experiment provides a separate check on the Monte Carlo standard-error formula.

The observed values do not match the theoretical predictions exactly because the Monte Carlo results depend on finite random samples. The overall trend, however, is consistent with the expected inverse-square-root relationship.

---

## 6. Multidimensional Comparison

### 6.1 Experimental Design

The main comparison examines both test functions across dimensions one through ten, with nominal budgets of 100, 1,000, and 10,000 evaluations. Each configuration uses 300 Monte Carlo trials.

For each configuration, the midpoint rule uses the largest complete grid that fits within the budget. Monte Carlo uses the full budget for every trial. The midpoint result is a deterministic approximation, while the Monte Carlo RMSE summarizes the distribution of stochastic estimates.

### 6.2 Sum-of-Squares Results

The exact integral is \(d/3\). At the largest budget, the results are:

| Dimension | Midpoint evaluations | Midpoint relative error | Monte Carlo relative RMSE |
|---:|---:|---:|---:|
| 1 | 10,000 | 0.00000025% | 0.9537% |
| 2 | 10,000 | 0.0025% | 0.5842% |
| 3 | 9,261 | 0.0567% | 0.5683% |
| 4 | 10,000 | 0.2500% | 0.4177% |
| 5 | 7,776 | 0.6944% | 0.4088% |
| 6 | 4,096 | 1.5625% | 0.3711% |
| 7 | 2,187 | 2.7778% | 0.3229% |
| 8 | 6,561 | 2.7778% | 0.2917% |
| 9 | 512 | 6.2500% | 0.2772% |
| 10 | 1,024 | 6.2500% | 0.2772% |

At a budget of 10,000, the midpoint rule has lower relative error through dimension four. Monte Carlo has lower relative RMSE beginning at dimension five. The observed crossover is therefore between dimensions four and five.

The midpoint errors do not increase perfectly smoothly with dimension. This is a consequence of the integer grid restriction. A small change in dimension can cause the largest permissible value of \(m\) to change abruptly, producing a substantial change in the number of midpoint evaluations.

For the sum-of-squares function, the exact integral increases linearly with dimension. The absolute Monte Carlo RMSE increases with dimension, while its relative RMSE decreases because the exact integral grows as well.

### 6.3 Product-of-Squares Results

The exact integral is \((1/3)^d\). At the largest budget, the results are:

| Dimension | Midpoint evaluations | Midpoint relative error | Monte Carlo relative RMSE |
|---:|---:|---:|---:|
| 1 | 10,000 | 0.00000025% | 0.8502% |
| 2 | 10,000 | 0.0050% | 1.5176% |
| 3 | 9,261 | 0.1700% | 2.2994% |
| 4 | 10,000 | 0.9963% | 3.0122% |
| 5 | 7,776 | 3.4243% | 4.1035% |
| 6 | 4,096 | 9.0163% | 5.9597% |
| 7 | 2,187 | 17.8970% | 7.8341% |
| 8 | 6,561 | 20.1777% | 10.4710% |
| 9 | 512 | 44.0575% | 13.8765% |
| 10 | 1,024 | 47.5540% | 18.0183% |

At a budget of 10,000, the midpoint rule has lower relative error through dimension five. Monte Carlo has lower relative RMSE beginning at dimension six. The observed crossover is therefore between dimensions five and six.

Unlike the sum-of-squares function, the exact integral of the product-of-squares function decreases exponentially with dimension. Its relative Monte Carlo RMSE consequently increases with dimension, even though the absolute variance of the integrand decreases.

### 6.4 Comparison of the Two Integrands

The two functions demonstrate that the effect of dimension depends on the mathematical structure of the integrand.

For the sum-of-squares function, the exact integral grows linearly with dimension. The relative Monte Carlo RMSE decreases with dimension, and the observed crossover occurs at a lower dimension.

For the product-of-squares function, the exact integral decreases exponentially. The relative Monte Carlo RMSE increases with dimension, and the midpoint rule remains competitive through a higher dimension.

The comparison shows why absolute error alone is insufficient when the scale of the exact integral changes. Relative error provides an additional perspective on the accuracy of each method.

### 6.5 Crossover Analysis

The crossover is defined as the point at which Monte Carlo RMSE becomes smaller than midpoint absolute error for a given configuration. Dividing both quantities by the same nonzero exact integral does not change which method has the smaller error, so the winner is identical when relative errors are compared.

At the 10,000-evaluation budget, the observed crossover occurs between dimensions four and five for the sum-of-squares function, and between dimensions five and six for the product-of-squares function.

At smaller budgets, the crossover can be less regular. The midpoint grid size changes in integer steps, so increasing the nominal budget does not always produce a proportional increase in actual midpoint evaluations. This can create local irregularities in the comparison.

### 6.6 Monte Carlo Bias and Variability

The multidimensional experiments also record the mean and standard deviation of the Monte Carlo estimates. The observed biases are generally small relative to the RMSE, consistent with the theoretical unbiasedness of the estimator.

The empirical standard deviations provide a check on the spread of the estimates. The finite number of trials means that these observed quantities fluctuate around their theoretical values.

### 6.7 Limitations

The multidimensional experiments use only two test functions. Although they have useful contrasting properties and exact integrals, they do not represent every type of high-dimensional integrand.

The deterministic method is specifically the tensor-product midpoint rule. Other deterministic methods, including sparse-grid methods and higher-order quadrature rules, may behave differently.

The nominal evaluation budgets are maximum budgets rather than equal actual evaluation counts. The midpoint method can use fewer evaluations than Monte Carlo because it requires a complete grid. This is an important feature of the comparison, but it also means that the methods are not always evaluated at identical sample counts.

Finally, the experiments cover dimensions one through ten and budgets no larger than 10,000. Conclusions beyond this range require additional experiments.

---

## 7. Runtime Comparison

### 7.1 Method

Runtime is measured separately from the accuracy experiments. For each dimension and budget, the implementation performs a warm-up evaluation and then records repeated batches of function evaluations. The median time per evaluation is reported.

The runtime experiment uses nine timing repetitions. Measurements are reported in milliseconds.

### 7.2 Results

| Dimension | Budget | Midpoint evaluations | Midpoint runtime (ms) | Monte Carlo runtime (ms) |
|---:|---:|---:|---:|---:|
| 1 | 100 | 100 | 0.010 | 0.005 |
| 1 | 1,000 | 1,000 | 0.013 | 0.007 |
| 1 | 10,000 | 10,000 | 0.031 | 0.026 |
| 2 | 100 | 100 | 0.016 | 0.005 |
| 2 | 1,000 | 961 | 0.024 | 0.014 |
| 2 | 10,000 | 10,000 | 0.097 | 0.108 |
| 3 | 100 | 64 | 0.021 | 0.006 |
| 3 | 1,000 | 1,000 | 0.030 | 0.017 |
| 3 | 10,000 | 9,261 | 0.117 | 0.127 |
| 4 | 100 | 81 | 0.024 | 0.006 |
| 4 | 1,000 | 625 | 0.032 | 0.017 |
| 4 | 10,000 | 10,000 | 0.137 | 0.128 |
| 5 | 100 | 32 | 0.027 | 0.006 |
| 5 | 1,000 | 243 | 0.031 | 0.019 |
| 5 | 10,000 | 7,776 | 0.135 | 0.149 |
| 6 | 100 | 64 | 0.033 | 0.007 |
| 6 | 1,000 | 729 | 0.047 | 0.021 |
| 6 | 10,000 | 4,096 | 0.105 | 0.169 |
| 7 | 100 | 1 | 0.021 | 0.007 |
| 7 | 1,000 | 128 | 0.041 | 0.024 |
| 7 | 10,000 | 2,187 | 0.081 | 0.431 |
| 8 | 100 | 1 | 0.020 | 0.007 |
| 8 | 1,000 | 256 | 0.046 | 0.025 |
| 8 | 10,000 | 6,561 | 0.190 | 0.210 |
| 9 | 100 | 1 | 0.023 | 0.007 |
| 9 | 1,000 | 512 | 0.061 | 0.026 |
| 9 | 10,000 | 512 | 0.064 | 0.519 |
| 10 | 100 | 1 | 0.025 | 0.007 |
| 10 | 1,000 | 1 | 0.026 | 0.029 |
| 10 | 10,000 | 1,024 | 0.083 | 0.556 |

### 7.3 Interpretation

Monte Carlo is generally faster in the smaller-budget measurements, while the results at the largest budget are mixed. In several higher-dimensional configurations, the midpoint method is faster because its complete grid contains far fewer points than the Monte Carlo sample count.

For example, at dimension nine and a nominal budget of 10,000, the midpoint grid contains only 512 points. Its runtime is therefore not a comparison against an equal number of function evaluations.

Runtime is also affected by the implementation, NumPy's vectorization, memory allocation, grid construction, hardware, and measurement precision. The values are specific to this computational environment and should not be interpreted as universal performance rankings.

The runtime experiment is supplementary. The accuracy and theoretical convergence analyses provide the primary basis for comparing the methods.

---

## 8. Discussion and Conclusion

### 8.1 Overall Findings

The experiments demonstrate that the relative performance of midpoint-rule integration and Monte Carlo integration changes with dimension, integrand, and evaluation budget.

In the one-dimensional baseline, the midpoint rule is substantially more accurate. Its observed convergence slope of \(-2.0000\) closely matches the theoretical rate of \(N^{-2}\), while the Monte Carlo slope of \(-0.5016\) matches its theoretical \(N^{-1/2}\) behavior.

In three dimensions, the observed slopes also agree closely with theory. The midpoint slope is \(-0.6667\), and the Monte Carlo slope is \(-0.4953\). The variability experiment independently supports the predicted Monte Carlo standard-error relationship.

In the multidimensional comparison at a budget of 10,000, the observed crossover occurs between dimensions four and five for the sum-of-squares function, and between dimensions five and six for the product-of-squares function.

### 8.2 The Role of the Integrand

The two test functions produce different relative-error trends. For the sum-of-squares function, the exact integral increases linearly with dimension. Its relative Monte Carlo RMSE decreases with dimension.

For the product-of-squares function, the exact integral decreases exponentially. Its relative Monte Carlo RMSE increases with dimension.

This contrast is central to interpreting high-dimensional integration. The dimension of the domain matters, but the function being integrated also determines the variance and scale of the estimator. A method that performs well for one integrand may not perform as well for another.

### 8.3 The Curse of Dimensionality

The tensor-product midpoint rule requires \(m^d\) evaluations for \(m\) subdivisions per coordinate. This exponential growth makes fine grids increasingly expensive as dimension increases.

Monte Carlo avoids constructing a full grid. Its theoretical convergence exponent remains \(-1/2\), although the constant in its error depends on the variance of the integrand. This makes Monte Carlo an important alternative when tensor-product grids become impractical.

However, a dimension-independent convergence exponent does not guarantee that Monte Carlo will always have lower error at a given budget. The product-of-squares experiment demonstrates that the relative error can still grow rapidly with dimension.

### 8.4 Accuracy and Computational Efficiency

The results do not identify a single method that is best in every setting. The midpoint rule performs very well on the smooth low-dimensional examples, while Monte Carlo becomes more competitive as the tensor-product grid becomes constrained by dimension.

The crossover points are not universal constants. They depend on the integrand, the evaluation budget, and the discrete grid sizes available to the midpoint rule.

The runtime experiment reinforces the need to distinguish theoretical complexity from measured execution time. Implementation details and actual evaluation counts influence observed runtime, so timing results should be interpreted alongside the accuracy measurements.

### 8.5 Limitations and Future Work

Several extensions could make the investigation more comprehensive:

- Test additional integrands with different smoothness, variance, and coordinate dependence.
- Compare the tensor-product midpoint rule with other deterministic methods, including sparse grids or higher-order quadrature.
- Examine larger dimensions and evaluation budgets where computational resources permit.
- Use equal actual evaluation counts in an additional comparison to complement the maximum-budget design.
- Investigate variance-reduction methods for Monte Carlo integration.
- Repeat runtime measurements on additional hardware and report the computational environment in greater detail.

These extensions would help determine how broadly the observed patterns apply beyond the two functions studied here.

### 8.6 Conclusion

This investigation compared the midpoint rule and Monte Carlo integration through mathematical analysis and computational experiments. The results support the theoretical convergence rates of both methods and show that their relative performance depends on more than dimension alone.

The midpoint rule is highly accurate in low-dimensional settings, but its tensor-product evaluation count grows exponentially with dimension. Monte Carlo maintains an \(N^{-1/2}\) convergence rate, making it increasingly competitive when deterministic grids become expensive. The two test functions produce different crossover points, demonstrating the importance of the integrand's mathematical structure.

The central conclusion is that choosing a numerical integration method requires considering the dimension, the integrand, the available evaluation budget, and the computational implementation together.

---

## 9. Reproducibility

### 9.1 Software

The experiments were written in Python and use NumPy for numerical computation and Matplotlib for figures.

Random-number generators are initialized with fixed seeds in the scripts. This makes the stochastic experiments reproducible when run with compatible software versions and numerical environments.

### 9.2 Experiment Configuration

| Experiment | Configuration |
|---|---|
| One-dimensional baseline | \(x^2\), sample sizes 10 through 10,000, 200 Monte Carlo trials |
| Three-dimensional convergence | Sum of squares, sample sizes 27 through 8,000 |
| Monte Carlo variability | Dimension 3, sample sizes 100 through 10,000, 500 trials |
| Multidimensional comparison | Dimensions 1 through 10, budgets 100, 1,000, and 10,000, 300 trials per configuration |
| Runtime comparison | Dimensions 1 through 10, budgets 100, 1,000, and 10,000, nine timing repetitions |

### 9.3 Source Files

The repository contains the following experiment scripts:

| File | Purpose |
|---|---|
| monte_carlo.py | One-dimensional baseline |
| multidimensional.py | Sum-of-squares multidimensional experiment |
| multidimensional_product.py | Product-of-squares multidimensional experiment |
| convergence_rates.py | Convergence-rate comparison |
| compare_experiments.py | Main multidimensional comparison |
| variability_analysis.py | Monte Carlo variability |
| runtime_comparison.py | Runtime measurements |

### 9.4 Running the Experiments

Clone the repository and install the dependencies:

    git clone https://github.com/jmathew89/integration-in-high-dimensions.git
    cd integration-in-high-dimensions
    pip install -r requirements.txt

Run an experiment from the repository root. For example:

    python src/compare_experiments.py

The other scripts can be run in the same way by replacing the filename.

### 9.5 Reproducibility Considerations

The fixed random seeds make it possible to reproduce the same pseudorandom sequences under compatible environments. Small differences in numerical output may occur across software versions or hardware.

The runtime measurements are particularly environment-dependent. They should be treated as observations from the machine on which the experiment was run, not as universal timing guarantees.

---

## References

This investigation uses standard results from numerical integration, probability, and Monte Carlo estimation. The mathematical derivations are provided in Section 2 so that the theoretical predictions used in the experiments can be checked directly.

No claim of a new numerical integration algorithm is made. The contribution of the project is the comparison of established methods through a controlled set of theoretical and computational experiments.
