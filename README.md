# statistics-in-java

## Explanation

The following gives a rough overview of how the statistics implementations in java and R differ. Particularly, the R
implementation has more capabilities out of the box while the Java implementation has to be used more explicitly.

## In the R source code regarding the t.test function and the pt function. What does alternative = greater mean. How is this verified in the source code?

To understand what `alternative = "greater"` means in the context of the `t.test` function in R, and how this is
implemented in the source code, we need to delve into the basics of hypothesis testing and the specifics of the t-test.

### Understanding `alternative = "greater"`

The `t.test` function in R performs a Student's t-test, which is used to determine whether the means of two groups
differ significantly. The `alternative` argument specifies the nature of the test being performed:

- `alternative = "two.sided"` (default): Tests the null hypothesis that the population means of the two groups are
  equal.
- `alternative = "less"`: Tests the null hypothesis that the population mean of one group is less than that of another.
- `alternative = "greater"`: Tests the null hypothesis that the population mean of one group is greater than that of
  another.

In simpler terms, `alternative = "greater"` is used when you hypothesize that the mean of one group is higher than the
mean of another group.

### Verifying in Source Code

To verify how `alternative = "greater"` is handled in the source code, we would typically look at the implementation of
the `t.test` function and related functions like `pt`, which computes the cumulative distribution function for the
t-distribution. However, without direct access to the current version of R's source code, I'll describe a general
approach based on common practices in statistical computing.

1. **Locate the Implementation**: First, find the source code for the `t.test` function. This might be located in a file
   named something like `stats.c` within the R source code directory.

2. **Identify the Test Type Handling**: Within the `t.test` function, there should be logic that handles
   the `alternative` argument. Look for conditional statements that check the value of `alternative`.

3. **Check for `alternative = "greater"` Handling**: Specifically, you're looking for a case where the code checks
   if `alternative == "greater"`. When this condition is met, the function will perform a one-sided t-test where the
   null hypothesis is that the mean of one group is not greater than the mean of the other group.

4. **Use of `pt` Function**: The `pt` function is likely used to compute the p-value for the test. The `pt` function
   takes a vector of quantiles and a degrees of freedom parameter and returns the cumulative distribution function
   values. For `alternative = "greater"`, the function would calculate the probability of observing a t-statistic as
   extreme or more so under the null hypothesis, given the alternative that the mean of one group is greater.

5. **Comparison with Null Hypothesis**: Finally, the result of the `pt` function (or similar) is compared against a
   significance level (commonly 0.05) to decide whether to reject the null hypothesis.

Without direct access to the R source code, this explanation provides a conceptual overview of
how `alternative = "greater"` might be implemented. For precise details, reviewing the actual source code of R would be
necessary.

Citations:

## this is the implementation of the pt function:

```c
/*

* R : A Computer Language for Statistical Data Analysis
* Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka
* Copyright (C) 2000-2007 The R Core Team
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, a copy is available at
* https://www.R-project.org/Licenses/
  */

#include "nmath.h"
#include "dpq.h"

double pt(double x, double n, int lower_tail, int log_p)
{
/* return P[ T <= x ]    where

* T ~ t_{n}  (t distrib. with n degrees of freedom).

* --> ./pnt.c for NON-central
  */
  double val, nx;
  #ifdef IEEE_754
  if (ISNAN(x) || ISNAN(n))
  return x + n;
  #endif
  if (n <= 0.0) ML_ERR_return_NAN;

if(!R_FINITE(x))
return (x < 0) ? R_DT_0 : R_DT_1;
if(!R_FINITE(n))
return pnorm(x, 0.0, 1.0, lower_tail, log_p);

#ifdef R_version_le_260
if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
/* Approx. from Abramowitz & Stegun 26.7.8 (p.949) */
val = 1./(4.*n);
return pnorm(x*(1. - val)/sqrt(1. + x*x*2.*val), 0.0, 1.0,
lower_tail, log_p);
}
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing rather than gaining precision,
     * now that pbeta(*, log_p = TRUE) is much better.
     * Note however that a version of this test *is* needed for x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
	/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
	   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
	   with z = 1/nx,  a = n/2,  b= 1/2 :
	*/
	double lval;
	lval = -0.5*n*(2*log(fabs(x)) - log(n))
		- lbeta(0.5*n, 0.5) - log(0.5*n);
	val = log_p ? lval : exp(lval);
    } else {
	val = (n > x * x)
	    ? pbeta (x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
	    : pbeta (1. / nx,             n / 2., 0.5, /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
    if(x <= 0.)
	lower_tail = !lower_tail;

    if(log_p) {
	if(lower_tail) return log1p(-0.5*exp(val));
	else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
	val /= 2.;
	return R_D_Cval(val);
    }

}
```

This is the implementation of the t.test function in r:

```r
# File src/library/stats/R/t.test.R on github

t.test <- function(x, ...) UseMethod("t.test")

t.test.default <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
...)
{
alternative <- match.arg(alternative)

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
        stop("'mu' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
    if( !is.null(y) ) {
	dname <- paste(deparse(substitute(x)),"and",
		       deparse(substitute(y)))
	if(paired)
	    xok <- yok <- complete.cases(x,y)
	else {
	    yok <- !is.na(y)
	    xok <- !is.na(x)
	}
	y <- y[yok]
    }
    else {
	dname <- deparse(substitute(x))
	if (paired) stop("'y' is missing for paired test")
	xok <- !is.na(x)
	yok <- NULL
    }
    x <- x[xok]
    if (paired) {
	x <- x-y
	y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if(is.null(y)) {
        if(nx < 2) stop("not enough 'x' observations")
	df <- nx-1
	stderr <- sqrt(vx/nx)
        if(stderr < 10 *.Machine$double.eps * abs(mx))
            stop("data are essentially constant")
	tstat <- (mx-mu)/stderr
	method <- if(paired) "Paired t-test" else "One Sample t-test"
	estimate <-
	    setNames(mx, if(paired)"mean of the differences" else "mean of x")
    } else {
	ny <- length(y)
        if(nx < 1 || (!var.equal && nx < 2))
            stop("not enough 'x' observations")
	if(ny < 1 || (!var.equal && ny < 2))
            stop("not enough 'y' observations")
        if(var.equal && nx+ny < 3) stop("not enough observations")
	my <- mean(y)
	vy <- var(y)
	method <- paste(if(!var.equal)"Welch", "Two Sample t-test")
	estimate <- c(mx,my)
	names(estimate) <- c("mean of x","mean of y")
	if(var.equal) {
	    df <- nx+ny-2
            v <- 0
            if(nx > 1) v <- v + (nx-1)*vx
            if(ny > 1) v <- v + (ny-1)*vy
	    v <- v/df
	    stderr <- sqrt(v*(1/nx+1/ny))
	} else {
	    stderrx <- sqrt(vx/nx)
	    stderry <- sqrt(vy/ny)
	    stderr <- sqrt(stderrx^2 + stderry^2)
	    df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
	}
        if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
	pval <- pt(tstat, df)
	cint <- c(-Inf, tstat + qt(conf.level, df) )
    }
    else if (alternative == "greater") {
	pval <- pt(tstat, df, lower.tail = FALSE)
	cint <- c(tstat - qt(conf.level, df), Inf)
    }
    else {
	pval <- 2 * pt(-abs(tstat), df)
	alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
	cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
	       conf.int = cint, estimate = estimate, null.value = mu,
	       alternative = alternative,
	       method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)

}

t.test.formula <-
function(formula, data, subset, na.action, ...)
{
if(missing(formula)
|| (length(formula) != 3L)
|| (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
stop("'formula' missing or incorrect")
m <- match.call(expand.dots = FALSE)
if(is.matrix(eval(m$data, parent.frame())))
m$data <- as.data.frame(data)

## need stats:: for non-standard evaluation

m[[1L]] <- quote(stats::model.frame)
m$... <- NULL
mf <- eval(m, parent.frame())
DNAME <- paste(names(mf), collapse = " by ")
names(mf) <- NULL
response <- attr(attr(mf, "terms"), "response")
g <- factor(mf[[-response]])
if(nlevels(g) != 2L)
stop("grouping factor must have exactly 2 levels")
DATA <- setNames(split(mf[[response]], g), c("x", "y"))
y <- do.call("t.test", c(DATA, list(...)))
y$data.name <- DNAME
if(length(y$estimate) == 2L)
names(y$estimate) <- paste("mean in group", levels(g))
y
}
```

How does the code make sure that one computed mean is bigger than the theoretical population mean?
The provided code snippets give us insight into how the `pt` function calculates the cumulative distribution function
for the t-distribution and how the `t.test` function uses this to perform a t-test, including handling the scenario
where the alternative hypothesis is that one mean is greater than the other (`alternative = "greater"`). Let's break
down how these pieces work together to ensure that one computed mean is considered larger than the theoretical
population mean.

### Understanding the `pt` Function

The `pt` function calculates the cumulative distribution function for the t-distribution. It takes four
parameters: `x` (the value at which to evaluate the CDF), `n` (degrees of freedom), `lower_tail` (a boolean indicating
whether to calculate the area to the left or right of `x`), and `log_p` (a boolean indicating whether to return the
result as a logarithm of the probability).

- When `lower_tail` is `TRUE` (the default), `pt` calculates the probability that a random variable following a
  t-distribution with `n` degrees of freedom is less than or equal to `x`.
- When `lower_tail` is `FALSE`, it calculates the probability that a random variable following a t-distribution with `n`
  degrees of freedom is greater than `x`.

### How `t.test` Uses `pt` for `alternative = "greater"`

In the `t.test` function, the calculation of the p-value depends on the `alternative` argument:

- For `alternative = "less"`, the p-value is calculated using the standard formula, assuming the null hypothesis that
  the means are equal.
- For `alternative = "greater"`, the function calculates the p-value differently. It uses the `pt` function
  with `lower_tail = FALSE` to calculate the probability that the absolute value of the t-statistic is greater than the
  observed value. This effectively tests the hypothesis that the mean of one group is greater than the mean of the
  other.

This is done in the following lines of the `t.test` function:

```r
if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df) )
}
else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
}
```

- For `alternative = "greater"`, `pt` is called with `lower.tail = FALSE`, meaning it calculates the probability that
  the t-statistic is greater than `x`. This is because the alternative hypothesis is that the mean of one group is
  greater than the mean of the other, so we want to know the probability of observing such a large t-statistic under the
  null hypothesis (that the means are equal).
- The confidence interval (`cint`) is then calculated accordingly, reflecting the fact that we are interested in values
  greater than the observed t-statistic.

### Ensuring One Computed Mean Is Greater Than the Theoretical Population Mean

The `t.test` function calculates the t-statistic as `(mx - my - mu) / stderr`, where `mx` and `my` are the sample means
of the two groups, `mu` is the hypothesized population mean under the null hypothesis, and `stderr` is the standard
error of the difference in means. This statistic measures how many standard errors away from the hypothesized mean the
observed means are.

By setting `lower.tail = FALSE` in the call to `pt` for `alternative = "greater"`, the function calculates the
probability of observing a t-statistic as extreme or more so in the direction of the alternative hypothesis (i.e.,
further from zero towards positive infinity), thus verifying that one computed mean is indeed greater than the
theoretical population mean under the null hypothesis.

## how can I achieve such an analysis in java? use as many standard libraries as possible such as org.apache.commons.math3.inference and if necessary, add custom code

To perform a t-test in Java using the Apache Commons Math library, specifically focusing on the scenario where you want
to test if one mean is greater than the other (similar to `alternative = "greater"` in R), you can follow these steps.
This guide assumes you have a basic understanding of Java and the Apache Commons Math library.

### Step 1: Add Apache Commons Math Dependency

First, ensure you have the Apache Commons Math library added to your project. If you're using Maven, add the following
dependency to your `pom.xml`:

```xml

<dependency>
    <groupId>org.apache.commons</groupId>
    <artifactId>commons-math3</artifactId>
    <version>3.6.1</version>
</dependency>
```

Adjust the version as necessary to match the latest release.

### Step 2: Perform the T-Test

You can use the `TTest` class from the `org.apache.commons.math3.stat.inference` package to perform the t-test. Here's
how you can do it:

```java
import org.apache.commons.math3.stat.inference.TTest;

public class TTestExample {
    public static void main(String[] args) {
        // Example data sets
        double[] sampleOne = {134, 146, 104, 119, 124, 161, 107, 83, 113, 129, 97, 123};
        double[] sampleTwo = {70, 118, 101, 85, 107, 132, 94};

        // Create a TTest instance
        TTest tTest = new TTest();

        // Perform the t-test
        double tStatistic = tTest.t(sampleOne, sampleTwo);

        System.out.println("T Statistic: " + tStatistic);
    }
}
```

This code snippet demonstrates how to perform a simple t-test comparing two samples. The `t` method of the `TTest` class
calculates the t-statistic, which measures the difference between the means of the two samples relative to the
variability of the data.

### Step 3: Interpret the Results

The t-statistic tells you how far apart the means of the two samples are, relative to the variability of the data. A
larger absolute value of the t-statistic indicates stronger evidence against the null hypothesis (that the means are
equal).

However, the Apache Commons Math library's `TTest` class does not directly support specifying the alternative hypothesis
as "greater". To analyze the results in the context of `alternative = "greater"`, you would typically compare the
t-statistic to the critical value from the t-distribution table corresponding to your chosen significance level (e.g.,
0.05 for a 95% confidence level). If the absolute value of the t-statistic exceeds this critical value, you would reject
the null hypothesis in favor of the alternative that one mean is greater than the other.

### Additional Considerations

- Ensure your data meets the assumptions of the t-test (e.g., normality, independence, homogeneity of variances).
- Adjust the significance level according to your needs.
- For more complex scenarios, such as unequal variances or paired samples, consult the documentation of the Apache
  Commons Math library or consider implementing custom logic.

This guide provides a starting point for performing a t-test in Java using the Apache Commons Math library. Depending on
your specific requirements, you may need to adapt or extend this example.

Citations:
- [1] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html?org/apache/commons/math3/stat/inference/TTest.html
- [2] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/inference/package-summary.html
- [3] https://stackoverflow.com/questions/47002374/apache-math-calculate-the-confidence-interval-in-a-t-test-which-has-two-sample
- [4] https://gist.github.com/d6e42f5ea538d7c4fd8f
- [5] https://forum.inductiveautomation.com/t/testing-the-org-apache-commons-math3-analysis-interpolation/74601
- [6] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/inference/OneWayAnova.html
- [7] https://www.baeldung.com/apache-commons-math
- [8] https://www.biostars.org/p/176092/
- [9] https://commons.apache.org/math/userguide/stat.html
- [10] https://discuss.gradle.org/t/problems-building-java-project-with-dependencies/20384

## be more specific, how can I use a t-distribution table in java to compare the t-statistic to check whether it exceeds the critical value. in particular, how can i get a critical value in java

To use a t-distribution table in Java to compare the t-statistic and check whether it exceeds the critical value,
especially for a one-sided test where you're interested in whether one mean is greater than the other, you can utilize
the Apache Commons Math library. This involves calculating the critical value for the desired significance level and
degrees of freedom, and then comparing it with your t-statistic.

Here's how you can achieve this:

### Step 1: Calculate the Critical Value

First, you need to calculate the critical value from the t-distribution table. Since you're doing a one-sided test,
you'll be looking up the critical value for the upper tail of the distribution. The critical value corresponds to the
significance level (e.g., 0.05 for a 95% confidence level).

You can use the `inverseCumulativeProbability` method from the `TDistribution` class in Apache Commons Math to find this
critical value. This method finds the t-value for a given probability in the distribution.

```java
import org.apache.commons.math3.distribution.TDistribution;

public class TTestAnalysis {
    public static void main(String[] args) {
        // Degrees of Freedom
        double degreesOfFreedom = 10; // Example: 10 degrees of freedom

        // Significance Level
        double significanceLevel = 0.05; // Example: 5%

        // Create a TDistribution object
        TDistribution tDist = new TDistribution(degreesOfFreedom);

        // Calculate the critical value for the upper tail
        double criticalValue = tDist.inverseCumulativeProbability(significanceLevel);

        System.out.println("Critical Value: " + criticalValue);
    }
}
```

### Step 2: Compare the T-Statistic with the Critical Value

After calculating the critical value, you can compare it with your t-statistic to see if it exceeds the critical value.
This comparison helps you determine whether to reject the null hypothesis in favor of the alternative that one mean is
greater than the other.

Assuming you have already calculated your t-statistic (let's say it's stored in a variable named `tStatistic`), you can
compare it as follows:

```java
// Assuming tStatistic is your calculated t-statistic
double tStatistic=...; // Your calculated t-statistic

// Check if the t-statistic exceeds the critical value
        if(tStatistic>criticalValue){
        System.out.println("Reject the null hypothesis. Evidence suggests one mean is greater.");
        }else{
        System.out.println("Do not reject the null hypothesis. There is not enough evidence to suggest one mean is greater.");
        }
```

### Summary

This approach allows you to perform a one-sided t-test in Java using the Apache Commons Math library, specifically
focusing on determining whether one mean is greater than the other. By calculating the critical value for the upper tail
of the t-distribution and comparing it with your t-statistic, you can make a decision based on the strength of evidence
against the null hypothesis.

Citations:
- [1] https://stackoverflow.com/questions/43859288/statistics-inquiry
- [2] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/TDistribution.html
- [3] https://commons.apache.org/math/userguide/stat.html
- [4] https://commons.apache.org/math/userguide/distribution.html
- [5] https://www.youtube.com/watch?v=IOn2nRBso00
- [6] https://ashwani-singh-nitk.medium.com/apache-commons-math-descriptive-statistics-55795af94565
- [7] https://github.com/apache/commons-statistics/blob/master/commons-statistics-distribution/src/main/java/org/apache/commons/statistics/distribution/TDistribution.java
- [8] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html?org/apache/commons/math3/util/FastMath.html
- [9] https://www.youtube.com/watch?v=4MlZxW117B4

## how do I choose degreesOfFreedom in TDistribution? what are the guidelines? I have a sample size of 1000

Choosing the degrees of freedom for a t-distribution in Java, particularly when working with the Apache Commons Math
library, is crucial for accurately interpreting the results of statistical tests. The degrees of freedom (df) in a
t-distribution are closely related to the sample size and the type of test being conducted. Here's how to choose them
based on your sample size of 1000:

### Guidelines for Choosing Degrees of Freedom

- **Independent Samples**: For independent samples, the degrees of freedom are usually calculated as the total sample
  size minus the number of groups. If you're comparing two groups, the degrees of freedom would be \(n - 2\), where
  \(n\) is the total sample size. For a sample size of 1000, comparing two groups, the degrees of freedom would be
  \(1000 - 2 = 998\).

- **Paired Samples**: For paired samples (e.g., pre-post measurements on the same subjects), the degrees of freedom are
  calculated slightly differently. They are the number of pairs minus 1. If you have 500 pairs, the degrees of freedom
  would be \(500 - 1 = 499\).

- **One-Sample Test**: For a one-sample t-test, the degrees of freedom are also calculated as the sample size minus 1.
  With a sample size of 1000, the degrees of freedom would be \(1000 - 1 = 999\).

### Implementing in Java Using Apache Commons Math

Given your sample size of 1000 and assuming you're conducting a two-sample t-test comparing two groups, here's how you
could calculate the degrees of freedom and use them with the `TDistribution` class:

```java
import org.apache.commons.math3.distribution.TDistribution;

public class TTestWithLargeSampleSize {
    public static void main(String[] args) {
        // Total sample size
        int totalSampleSize = 1000;

        // Number of groups
        int numberOfGroups = 2;

        // Calculate degrees of freedom
        double degreesOfFreedom = totalSampleSize - numberOfGroups + 1;

        // Create a TDistribution object
        TDistribution tDist = new TDistribution(degreesOfFreedom);

        // Now you can use tDist for various purposes, e.g., calculating critical values
        double criticalValue = tDist.inverseCumulativeProbability(0.975); // For a 95% confidence interval
        System.out.println("Critical Value: " + criticalValue);
    }
}
```

This code snippet demonstrates how to calculate the degrees of freedom for a two-group comparison from a sample size of
1000 and then use those degrees of freedom to create a `TDistribution` object.
From there, you can proceed to calculate critical values, probabilities, etc., as needed for your statistical analysis.

Remember, the choice of degrees of freedom is fundamental to the validity of your statistical conclusions, so it's
important to select it correctly based on the structure of your data and the type of test you're performing.

Citations:
- [1] https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/TDistribution.html
- [2] https://commons.apache.org/math/userguide/stat.html
- [3] https://github.com/apache/commons-statistics/blob/master/commons-statistics-distribution/src/main/java/org/apache/commons/statistics/distribution/TDistribution.java
- [4] https://stackoverflow.com/questions/47002374/apache-math-calculate-the-confidence-interval-in-a-t-test-which-has-two-sample
- [5] http://repository.transtep.com/repository/thirdparty/commons-math-1.2/docs/apidocs/index.html?org/apache/commons/math/distribution/TDistributionImpl.html
- [6] https://javadoc.scijava.org/Apache-Commons-Math/org/apache/commons/math3/random/EmpiricalDistribution.html
- [7] https://www.youtube.com/watch?v=MifbwoCumX0
- [8] https://www.baeldung.com/apache-commons-math
- [9] https://academic.oup.com/isd/article/6/1/4/6514767