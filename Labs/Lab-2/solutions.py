# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# %%
df = pd.read_csv("krw_data.csv")
print(df.shape)
df.sort_values(by=["firm", "job", "white"], inplace=True)
df.head()

# %%
print("Number of firms = ", len(df.firm.unique()))
print("Number of jobs = ", len(df.job.unique()))


# %%
df.describe()

# %% [markdown]
# Question 2 - Regress an indicator for a callback on an indicator for a white name in the pooled sample of all
# applications. What is the average effect of a white name on callbacks? Compare robust and jobclustered
# standard errors for the race coefficient. How do these standard errors differ, and why?

# %%
# Using robust standard errors
model = smf.ols("callback ~ white", data=df)
fitted_model = model.fit(cov_type="HC1")
fitted_model.summary()

# %%
# Using jobclustered standard errors for the race coefficient
model = smf.ols("callback ~ white", data=df)
fitted_model = model.fit(cov_type="cluster", cov_kwds={"groups": df["job"]})
fitted_model.summary()


# %% [markdown]
# Clustered standard errors are lower than robust standard errors.

# %% [markdown]
# Question 3: Compute the white/Black difference in callback rates separately for each job vacancy in the data set.
# Let ˆΔjf denote the contact gap for job j within firm f. Take the average of ˆΔjf for each firm, resulting
# in 108 firm-specific estimates ˆθf . This is an unbiased estimate of the average effect of race at firm f,
# labeled θf .
#
#
# Question 4: Compute a standard error for each ˆθf as sf, where nf is the number
# of jobs for firm f. Collapse the data down to a firm-level data set with 108 observations on ˆθf and sf .

# %%
call_back_df = pd.pivot_table(
    df,
    values="callback",
    index=["firm", "job"],
    columns=["white"],
    aggfunc=[np.mean, "count"],
)
# call_back_df.reset_index(inplace=True)
call_back_df.columns = [
    "black_call_back_rate",
    "white_call_back_rate",
    "black_count",
    "black_count",
]
call_back_df.reset_index(inplace=True)

call_back_df["delta_hat_jf"] = (
    call_back_df["white_call_back_rate"] - call_back_df["black_call_back_rate"]
)

call_back_df.dropna(subset=["delta_hat_jf"], inplace=True)

print(call_back_df.columns)
call_back_df.head()

# %%
dictionary_of_sample_variance = (
    call_back_df.groupby("firm")["delta_hat_jf"].agg(np.var).to_dict()
)
dictionary_of_sample_variance[1]

# %% [markdown]
# EB Step 1: Compute \theta_hat_f and s_f^2
#
#
#

# %%
firm_level_theta_hat = call_back_df.groupby("firm").agg(
    {"delta_hat_jf": "sum", "job": "count"}
)
firm_level_theta_hat["theta_hat"] = (
    firm_level_theta_hat["delta_hat_jf"] / firm_level_theta_hat["job"]
)
firm_level_theta_hat.reset_index(inplace=True)

firm_level_theta_hat["sample_variance"] = (
    firm_level_theta_hat.apply(
        lambda x: dictionary_of_sample_variance[x["firm"]], axis=1
    )
    / firm_level_theta_hat["job"]
)
firm_level_theta_hat["sample_standard_deviation"] = (
    firm_level_theta_hat["sample_variance"] ** 0.5
)
print("Number of firms = ", firm_level_theta_hat.shape[0])
firm_level_theta_hat.head(8)

# %%
firm_level_theta_hat.theta_hat.hist(bins=50)
plt.show()

# %% [markdown]
# Question 5: Suppose we view the θf ’s as random draws from a mixing distribution G. Using your unbiased estimates
# and standard errors, compute a bias-corrected estimate of the variance of G, labeled ˆσ2
# θ . How does
# the standard deviation ˆσθ compare to the standard deviation of unbiased estimates ˆθf ? In economic
# terms, is ˆσθ big or small?


# %%
def mixing_fn(theta, sample_se):
    mu_theta = np.mean(theta)
    print("mu_theta = ", mu_theta)
    var_theta = np.sum(
        (theta - mu_theta) ** 2 - (sample_se**2) * ((len(theta) - 1) / len(theta))
    ) / len(theta)
    sigma_theta = np.sqrt(var_theta)
    print("sigma_theta = ", sigma_theta)
    return (mu_theta, sigma_theta)


# %%
mu_theta_hat, sigma_theta_hat = mixing_fn(
    theta=firm_level_theta_hat.theta_hat.values,
    sample_se=firm_level_theta_hat["sample_standard_deviation"].values,
)
mu_theta_hat - 1.96 * sigma_theta_hat, mu_theta_hat + 1.96 * sigma_theta_hat

# %%
firm_level_theta_hat.head(8)

# %% [markdown]
# Question 6: Form linear shrinkage estimates of the effect of race at each firm. Plot histograms of unbiased and
# linear shrinkage estimates.


# %%
def compute_theta_star(mu_hat, theta_hat, sigma_muhat, sample_se):
    snr = (sigma_muhat**2) / (
        sigma_muhat**2 + sample_se**2
    )  # signal to noise ratio
    theta_hat_star = snr * theta_hat + (1 - snr) * mu_hat
    return theta_hat_star


# %%
firm_level_theta_hat["theta_hat_star"] = compute_theta_star(
    mu_theta_hat,
    firm_level_theta_hat.theta_hat,
    sigma_theta_hat,
    firm_level_theta_hat.sample_standard_deviation,
)
firm_level_theta_hat.head()

# %%
plt.hist(
    [firm_level_theta_hat.theta_hat.values, firm_level_theta_hat.theta_hat_star.values],
    bins=30,
    label=["theta_hat", "theta_hat_star"],
)
plt.legend(loc="upper right")
plt.show()

# %% [markdown]
# Question 7: The last part of this lab asks you to compute a non-parametric deconvolution estimate of G, the
# distribution of discrimination across firms. This is an advanced exercise, so do not worry if you find it
# difficult.

# %%
