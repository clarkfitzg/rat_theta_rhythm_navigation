# What does the spectral analysis look like if there actually *is* a signal?
source("functions.R")

set.seed(3480)

minutes = 8
obs_per_minute = 1000
hz = 8
signal_noise_ratio = 1

zeros = round(obs_per_minute / hz)

total_ones = hz * minutes

# This signal follows the period perfectly
x = c(rep(0, times = zeros), 1)
x = rep(x, times = total_ones)
n = length(x)

# The ones are so sparse I'm not going to worry about converting those to 0's
noise = sample.int(n, size = signal_noise_ratio * total_ones)
x[noise] = 1

a = acf(x, lag.max = 400, plot = FALSE)
a = a$acf[-1]
plot(a, type = "l")

