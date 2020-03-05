source("functions.R")

# Compare ACF before and after injection

PreInjection = read_one_spike_train_file("IP18.12-9-16.PreInjection.Standard.Cylinder.2016-12-09_14-31-21.SE9.nse.converted.SpikeSortedForAutocorrelation")
PostInjection = read_one_spike_train_file("IP18.12-9-16.PostInjection.Standard.Cylinder.2016-12-09_19-18-01.SE15.nse.converted.SpikeSortedForAutocorrelation")
post21 = read_one_spike_train_file("IP21.3-29-17.PostInjection.Standard.Cylinder.2017-03-29_20-02-34.SE16.nse.converted.SpikeSortedForAutocorrelation")


# The pre injection ACF is around 0.002, while the post injection ACF is around 0.01, around five times larger.
# What could that mean?
pdf("pre_post_acf.pdf", width = 12, height = 12)
par(mfrow = c(3, 1))
plot_acf(PreInjection, main = "IP18.12-9-16.PreInjection")
plot_acf(PostInjection, main = "IP18.12-9-16.PostInjection")
plot_acf(post21, main = "IP21.3-29-17.PostInjection")
dev.off()

if(FALSE)
{

# R's default confidence interval lines are at 0.0026 for the pre injection ACF.
# This means we cannot really detect any autocorrelation at all in the PreInjection data- it's not large enough to distinguish it from noise.
a = acf(PreInjection$spike)
acf_clim(a)

f = fft(PreInjection$spike)
plot.frequency.spectrum(f, xlimits = c(0, 40))

spectrum(PreInjection$spike)
}
