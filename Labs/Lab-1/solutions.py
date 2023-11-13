# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf

# %%
df = pd.read_csv("vam_example_data.csv")
print(df.shape)
df.head()

# %%
df.sort_values(by="D", inplace=True)
df.head()

# %%
np.sort(df.D.unique())

# %%
df.D = df.D.astype(str)
df.dtypes

# %%
df.describe()

# %%
df.D.astype(int).hist()

# %%
df.head(2)


# %%
def reg(formula, data):
    model = smf.ols(formula=formula, data=data)
    fitted_model = model.fit(cov_type="HC1")
    data1 = {
        "coef": fitted_model.params,
        "std err": fitted_model.bse,
        "t": fitted_model.tvalues,
        "P>|t|": fitted_model.pvalues,
        "[0.025": fitted_model.conf_int()[0],
        "0.975]": fitted_model.conf_int()[1],
    }
    return pd.DataFrame(data1)


# %%
ols_uncontrolled = reg(formula="Y ~ D + 0", data=df)
print(ols_uncontrolled.shape)
print(ols_uncontrolled.head())

# %%
form = "Y ~ D + X + 0"
ols_controlled = reg(formula=form, data=df)
print(ols_controlled.shape)
print(ols_controlled.head())
print(ols_controlled.index)

# %%
df["dummy"] = df.D.apply(lambda x: "D[" + x + "]")

# %%
df_group = df.groupby("dummy").agg({"theta_D": np.mean})
df_group.head()

# %%
ols_uncontrolled = pd.concat([ols_uncontrolled, df_group], axis=1)
ols_uncontrolled.head()

# %%
df[df.D == "1"].head()


# %%
def mixing_fn(theta, sample_se):
    mu_theta = np.mean(theta)
    print("mu_theta = ", mu_theta)
    var_theta = np.sum((theta - mu_theta) ** 2 - sample_se**2) / len(theta)
    sigma_theta = np.sqrt(var_theta)
    print("sigma_theta = ", sigma_theta)
    return (mu_theta, sigma_theta)


# %%
mu_theta_hat_uncontrolled, sigma_theta_hat_uncontrolled = mixing_fn(
    theta=ols_uncontrolled.coef.values, sample_se=ols_uncontrolled["std err"].values
)
(
    mu_theta_hat_uncontrolled - 1.96 * sigma_theta_hat_uncontrolled,
    mu_theta_hat_uncontrolled + 1.96 * sigma_theta_hat_uncontrolled,
)

# %%
np.mean(ols_uncontrolled["theta_D"])

# %%
ols_controlled = ols_controlled.drop("X")
ols_controlled = pd.concat([ols_controlled, df_group], axis=1)
ols_controlled.head()

# %%
mu_theta_hat_controlled, sigma_theta_hat_controlled = mixing_fn(
    theta=ols_controlled.coef.values, sample_se=ols_controlled["std err"].values
)
(
    mu_theta_hat_controlled - 1.96 * sigma_theta_hat_controlled,
    mu_theta_hat_controlled + 1.96 * sigma_theta_hat_controlled,
)


# %%
np.mean(ols_uncontrolled["theta_D"])


# %%
def compute_theta_star(mu_hat, theta_hat, sigma_muhat, sample_se):
    snr = sigma_muhat**2 / (sigma_muhat**2 + sample_se**2)
    theta_hat_star = snr * theta_hat + (1 - snr) * mu_hat
    return theta_hat_star


# %%
ols_controlled["theta_hat_star"] = compute_theta_star(
    mu_theta_hat_controlled,
    ols_controlled.coef,
    sigma_theta_hat_controlled,
    ols_controlled["std err"],
)

# %%
# both the estimates contain mean of the mu_theta
# variance dropped quite a bit for sigma_theta in second case

# %%
ols_controlled["squared_diff_theta_hat"] = (
    ols_controlled["theta_D"] - ols_controlled["coef"]
) ** 2
ols_controlled["squared_diff_theta_star"] = (
    ols_controlled["theta_D"] - ols_controlled["theta_hat_star"]
) ** 2

# %%
mu_theta_hat_controlled

# %%
np.mean(ols_controlled["std err"])

# %%
pd.options.display.float_format = "{:.5f}".format
ols_controlled.head(10)

# %%
# check the MSE of theta_hat and theta_hat_star, with respect to theta_D. The mean of the squared_diff columns in the below table gives the MSE
ols_controlled.describe()

# %%
plt.scatter(
    ols_controlled.coef,
    ols_controlled.theta_hat_star,
    label="Estimates",
    marker="o",
    color="blue",
)
plt.plot(
    [ols_controlled.coef, ols_controlled.theta_hat_star],
    [ols_controlled.coef, ols_controlled.theta_hat_star],
    linestyle="--",
    color="red",
)
plt.axhline(
    y=mu_theta_hat_controlled, linestyle=":", color="red", label="mu_theta", linewidth=1
)
plt.xlabel("theta hat")
plt.ylabel("theta hat star")
plt.legend()
plt.title("theta hat vs theta hat star")
plt.show()

# %%
sns.kdeplot(x=ols_controlled["theta_D"], label="theta_D")
sns.kdeplot(x=ols_controlled["theta_hat_star"], label="theta_hat_star")
sns.kdeplot(x=ols_controlled["coef"], label="theta_hat")
plt.axvline(
    np.mean(ols_controlled["theta_D"]), color="blue", linestyle="dashed", linewidth=0.8
)
plt.axvline(
    np.mean(ols_controlled["theta_hat_star"]),
    color="red",
    linestyle="dashed",
    linewidth=0.5,
)
plt.axvline(
    np.mean(ols_controlled["coef"]), color="green", linestyle="dashed", linewidth=0.8
)

plt.legend()
plt.show()

# %%
