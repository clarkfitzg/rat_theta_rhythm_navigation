source("functions.R")

# Compare ACF before and after injection

PreInjection = read_one_spike_train_file("IP18.12-9-16.PreInjection.Standard.Cylinder.2016-12-09_14-31-21.SE9.nse.converted.SpikeSortedForAutocorrelation")

PostInjection = read_one_spike_train_file("IP18.12-9-16.PostInjection.Standard.Cylinder.2016-12-09_19-18-01.SE15.nse.converted.SpikeSortedForAutocorrelation")

# The pre injection ACF is around 0.002, while the post injection ACF is around 0.01, around five times larger.
# What could that mean?
par(mfrow = c(2, 1))
plot_acf(PreInjection, lag.max = 800, main = "Pre Injection")
plot_acf(PostInjection, lag.max = 800, main = "Post Injection")

# R's default confidence interval lines are at 0.0026 for the pre injection ACF.
# This means we cannot really detect any autocorrelation at all in the PreInjection data.
a = acf(PreInjection$spike)
acf_clim(a)


