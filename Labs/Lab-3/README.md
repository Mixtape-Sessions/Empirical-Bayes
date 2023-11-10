# Coding Lab 3: Classifying Firm-level Discrimination

[[PDF of Instructions](https://github.com/Mixtape-Sessions/Empirical-Bayes/raw/main/Lab/Lab-3/lab3.pdf)]

This coding lab uses the data from the employment discrimination experiment of Kline, Rose and Walters (2022) to conduct an empirical Bayes multiple testing analysis. The objective of the analysis is to determine which firms can be classified as discriminating against Black applicants while controlling the False Discovery Rate (FDR).

1.  Replicate steps 1-4 from coding lab 2 to obtain a data set of 108 firm-specific discrimination estimates $\hat{\theta}_{f}$ along with standard errors $s_{f}$.

2.  For each firm, form a test statistic $z_{f} = \hat{\theta}_{f} / s_{f}$, and use this statistic to compute the $p$-value $p_{f}$ from a one-tailed test of $H_{0}: \theta_{f} = 0$ vs. $H_{A}: \theta_{f} > 0$.

3.  Plot a histogram of the firm-specific $p$-values. What do you notice about this distribution?

4.  Let $\pi_{0}=\Pr[\theta_{f}=0]=\int1[\theta=0]dG(\theta)$ denote the share of firms in the population that are not discriminating. Use your $p$-values to compute an upper bound $\hat{\pi}_{0}$ on $\pi_{0}$. [Hint: compute the average height of the $p$-value density above some threshold $\lambda$. You can start by setting $\lambda=0.5$, though any $\lambda$ will do.] Interpret your estimated bound.

5.  Use the bound from part (8) to compute $q$-values as $q_{f} = [p_{f} \hat{\pi}_{0}] / \hat{F}_{p}(p_{f})$, where $\hat{F}_{p}$ is the estimated CDF of $p$-values. Recall that rejecting all null hypotheses with $q$-values below $\bar{q}$ controls FDR to at most $100\bar{q}\%$. Limiting FDR to 5%, how many firms can you classify as discriminating against Black applicants?

Solutions are in Lab-2 folder!
