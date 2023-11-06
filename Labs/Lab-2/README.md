# Coding Lab 2: Labor Market Discrimination 

[[PDF of Instructions](https://github.com/Mixtape-Sessions/Empirical-Bayes/raw/main/Lab/Lab-2/lab2.pdf)]

This coding lab will walk you through an empirical Bayes analysis of employer heterogeneity in labor market discrimination using data from the field experiment of Kline, Rose and Walters (2022). This experiment submitted 83,643 fictitious applications to 11,114 real job vacancies within 108 Fortune 500 firms. Each application was randomly assigned a distinctively-white or distinctively-Black name, stratified by job so that each vacancy received 4 white and 4 Black applications (though a few vacancies closed before all 8 applications could be sent). The key outcome is whether the application received a callback from the employer within 30 days.

1.  Import the data set `krw_data.csv` from the course website into a statistical software package of your choice (I recommend Stata or R). This is an application-level data set including job and firm identifiers, an indicator for a distinctively-white name, and a callback indicator. Summarize the variables in this data set.

2.  Regress an indicator for a callback on an indicator for a white name in the pooled sample of all applications. What is the average effect of a white name on callbacks? Compare robust and job-clustered standard errors for the race coefficient. How do these standard errors differ, and why?

3.  Compute the white/Black difference in callback rates separately for each job vacancy in the data set. Let $\hat{\Delta}_{jf}$ denote the contact gap for job $j$ within firm $f$. Take the average of $\hat{\Delta}_{jf}$ for each firm, resulting in 108 firm-specific estimates $\hat{\theta}_{f}$ . This is an unbiased estimate of the average effect of race at firm $f$, labeled $\theta_{f}$.

4.  Compute a standard error for each $\hat{\theta}_{f}$ as
    $s_{f}=\sqrt{\tfrac{1}{n_{f}(n_{f}-1)}\sum_{j=1}^{n_{f}}(\hat{\Delta}_{jf}-\hat{\theta}_{f})^{2}}$,
    where $n_{f}$ is the number of jobs for firm $f$. Collapse the data
    down to a firm-level data set with 108 observations on
    $\hat{\theta}_{f}$ and $s_{f}$.

5.  Suppose we view the $\theta_{f}$'s as random draws from a mixing
    distribution $G$. Using your unbiased estimates and standard errors,
    compute a bias-corrected estimate of the variance of $G$, labeled
    $\hat{\sigma}_{\theta}^{2}$. How does the standard deviation
    $\hat{\sigma}_{\theta}$ compare to the standard deviation of
    unbiased estimates $\hat{\theta}_{f}$? In economic terms, is
    $\hat{\sigma}_{\theta}$ big or small?

6.  Form linear shrinkage estimates of the effect of race at each firm. Plot histograms of unbiased and linear shrinkage estimates.

7.  The last part of this lab asks you to compute a non-parametric deconvolution estimate of $G$, the distribution of discrimination across firms. This is an advanced exercise, so do not worry if you find it difficult.

(a) Convert the estimates for each firm to a $z$-score,
        $z_{f}=\hat{\theta}_{f}/s_{f}$. Assume $z_{f}\sim N(\mu_{f},1)$,
        where $\mu_{f}=\theta_{f}/s_{f}$.

(b) Compute a log-spline deconvolution estimate of the distribution of $\mu_{f}$ across firms. [Hint: This can be done in R with the `deconvolveR` package.]

(c) Compute a kernel density estimate of the distribution of log standard errors, $\log s_{f}$. [Hint: This can be done in R with the `density` command in the `stats` package.]

(d) If $\mu_{f}$ and $\log s_{f}$ are independent, the density function for $\theta_{f}=\mu_{f}\exp(\log s_{f})$ is given by:

$$
  g_{\theta}(\theta) = \int g_{\mu}\left(\theta\exp(-t)\right)f(t)\exp(-t) dt, 
$$

where $g_{\mu}$ is the density function for $\mu_{f}$ and $f$ is the density function for $\log s_{f}$. Use this expression together with your results from parts (b) and (c) to compute an estimate of the distribution of $\theta_{f}$ across firms. Overlay this distribution on the histogram of unbiased estimates from part (6).

**Extra credit:** Use your log-spline estimates to compute non-parametric posterior mean estimates for each $\theta_{f}$. Plot these against the linear shrinkage estimates from part (6). What do you make of any differences between these estimates?
