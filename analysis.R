source("functions.R")

# Compare ACF before and after injection

PreInjection = read_one_spike_train_file("IP18.12-9-16.PreInjection.Standard.Cylinder.2016-12-09_14-31-21.SE9.nse.converted.SpikeSortedForAutocorrelation")
PostInjection = read_one_spike_train_file("IP18.12-9-16.PostInjection.Standard.Cylinder.2016-12-09_19-18-01.SE15.nse.converted.SpikeSortedForAutocorrelation")
post21 = read_one_spike_train_file("IP21.3-29-17.PostInjection.Standard.Cylinder.2017-03-29_20-02-34.SE16.nse.converted.SpikeSortedForAutocorrelation")


# The pre injection ACF is around 0.002, while the post injection ACF is around 0.01, around five times larger.
# What could that mean?
pdf("pre_post_acf.pdf", width = 12, height = 12)
par(mfrow = c(3, 1))
ac_pre = plot_acf(PreInjection, main = "IP18.12-9-16.PreInjection")
ac_post = plot_acf(PostInjection, main = "IP18.12-9-16.PostInjection")
ac_post21 = plot_acf(post21, main = "IP21.3-29-17.PostInjection")
dev.off()

theta_index(ac_pre)
theta_index(ac_post)
theta_index(ac_post21)


if(FALSE)
{

# Trying to fit the Theta index to the autocorrelogram
o_pre = optim_theta_index(ac_pre, control = list(trace = 10, reltol = 1e-10, maxit = 1e4))

te_pre = theta_index(o_pre)


o_post = optim_theta_index(ac_post
, control = list(trace = 10, reltol = 1e-10, maxit = 1e4)
)

# How robust are these to initial parameters?

o_post2 = optim_theta_index(ac_post, par = c(a = 1, b = 1, c = 0.1, omega = default_omega, tau1 = 10, tau2 = 1)
, control = list(trace = 10, reltol = 1e-10, maxit = 1e4)
)

# Some of the parameters differ by orders of magnitude, based on the starting values for the optimization.
# This means calculating theta_index will be unreliable.
rel_diff(o_post$par, o_post2$par)

theta_index(o_post)
theta_index(o_post2)

# These values are quite different compared to the relative tolerance used in the numerical minimzation of 1e-10, which probably means I'm not truly finding the global minimum to this function.
rel_diff(o_post$value, o_post2$value)


# R's default confidence interval lines are at 0.0026 for the pre injection ACF.
# This means we cannot really detect any autocorrelation at all in the PreInjection data- it's not large enough to distinguish it from noise.
a = acf(PreInjection$spike)
acf_clim(a)

f = fft(PreInjection$spike)
plot.frequency.spectrum(f, xlimits = c(0, 40))

# This article suggests that standard Fourier spectrum works fine with enough data.
# https://www.tandfonline.com/doi/abs/10.1076/0929-1016(200010)31%3A4%3B1-2%3BFT481
spectrum(PreInjection$spike)


# The drug has cut the number of spikes down, about in half.

mean(PreInjection$spike)
# [1] 0.00357698
mean(PostInjection$spike)
# [1] 0.001930448

# With these numbers the highest possible frequency can be:
1000 * mean(PreInjection$spike)
# 3.58 hz

# Let's plot the indices of the spikes.
plot(which(as.logical(PreInjection$spike)))
# post21 has pronounced stair pattern, which means that the cell often does not fire for long periods of time, between 20k to 25k intervals, which is 20 to 25 seconds.

# PostInjection also has this stair pattern, with around 20k to 55k intervals.

# If there is a rhythm, then this plot will be approximately linear.
# If the plot is not linear, then there cannot be a rhythm.

s = spectrum(PreInjection$spike, spans = 201)


}
