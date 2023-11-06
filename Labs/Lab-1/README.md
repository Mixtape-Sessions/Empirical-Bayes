# Coding Lab 1: School Value-Added

[[PDF of Instructions](https://github.com/Mixtape-Sessions/Empirical-Bayes/raw/main/Lab/Lab-1/lab1.pdf)]

This coding lab will walk you through an example of empirical Bayes estimation of school value-added based on simulated data. Consider a population of students, each attending one of $J$ schools. The variable $D_{i} \in \{1,....,J\}$ indicates the school attended by student $i$. We are interested in the following constant-effects causal value-added model:
$$
  Y_{i}={\displaystyle \sum_{j=1}^{J}\theta_{j}D_{ij}+\beta X_{i}+\epsilon_{i}}
$$,
where $Y_{i}$ is a test score outcome for student $i$, $\theta_{j}$ is the causal effect of school $j$ and $D_{ij}=1\{D_{i}=j\}$ indicates attendance at $j$, $X_{i}$ is an observed control variable (e.g. a lagged test score), and $\epsilon_{i}$ represents unobserved determinants of students' potential outcomes and is assumed to satisfy $E[X_{i} \epsilon_{i}] = E[D_{ij} \epsilon_{i}] = 0\ \forall j$. The school-level parameters $\theta_{j}$ are assumed to be drawn from a normal mixing distribution:
$$
  \theta_{j}\sim N(\mu_{\theta},\sigma_{\theta}^{2}).
$$

<br/>
Our goal is to learn about the parameters of the mixing distribution and form accurate estimates of the individual school value-added parameters $\theta_{j}$.

1. Import the data set `vam_example_data.csv` from the course website into a statistical software package of your choice (I recommend Stata or R). This data set includes observations on $Y_{i},$ $D_{i}$, and $X_{i}$ simulated from the model above for 2,500 students, each attending one of $J=50$ schools. Since the data were simulated from a known data generating process, it also includes the true value-added of each student's school, given by $\theta_{d(i)}=\sum_{j}\theta_{j}D_{ij}$, which would not be known in a real-world application. Summarize the variables in this data set.

2. Create the school dummy variables $D_{ij}$. Fit two value-added models by ordinary least squares (OLS):

  (a) An *uncontrolled* OLS regression of $Y_{i}$ on the set of $D_{ij}$'s with no constant.

  (b) A *controlled* OLS regression of $Y_{i}$ on the $D_{ij}$'s, controlling for $X_{i}$ (again with no constant).

For each of these two models, collect the list of estimated value-added coefficients $\hat{\theta}_{j}
$ along with their robust standard errors $s_{j}$ for each school. Then collapse the data down to a  school-level data set with 50 observations on true value-added $\theta_{j}$, value-added estimates $\hat{\theta}_{j}$, and standard errors $s_{j}$.

3. For each of the two value-added models, use the OLS estimates and standard errors to form estimates $\hat{\mu}_{\theta}$ and $\hat{\sigma}_{\theta}^{2}$ of the mean and variance of value-added. How do the estimated variances differ for the controlled and uncontrolled models? What do you conclude from this comparison? 

4. Focus on estimates of the controlled model for the remainder of the question. Using these estimates, form linear shrinkage posteriors

$$
  \hat{\theta}_{j}^{*}=\left(\tfrac{\hat{\sigma}_{\theta}^{2}}{\hat{\sigma}_{\theta}^{2}+s_{j}^{2}}\right)\hat{\theta}_{j}+\left(\tfrac{s_{j}^{2}}{\hat{\sigma}_{\theta}^{2}+s_{j}^{2}}\right)\hat{\mu}_{\theta}.
$$

Summarize your results by making a plot that includes histograms of the OLS and shrunk estimates overlayed with the estimated mixing distribution (i.e. $N(\hat{\mu}_{\theta},\hat{\sigma}_{\theta}^{2})$). Compare standard deviations of the OLS estimates $\hat{\theta}_{j}$, the shrunk posteriors $\hat{\theta}_{j}^{*}$, the estimated mixing distribution, and the true-value-added parameters $\theta_{j}$. Which standard deviation is biggest, and which is smallest? Provide some intuition.

5. For each school, compute the squared difference between true value-added $\theta_{j}$ and the OLS estimate $\hat{\theta}_{j}$, as well as the squared difference between true value-added and the linear shrinkage estimate $\hat{\theta}_{j}^{*}$. Which estimator has lower mean squared error (MSE) averaged across the 50 schools in the sample? Does shrinkage cause squared error to move in the same direction for all schools? Give some intuition for your results.

