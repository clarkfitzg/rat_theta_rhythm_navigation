# What does the ACF and periodogram look like if there actually *is* a signal?

# It might be the case that when we add some irregularity to the signal, the ACF just doesn't really work anymore.

set.seed(3480)

minutes = 8
obs_per_second = 1000
hz = 6

total_obs = minutes * 60 * obs_per_second
# Hmm... 

# Add noise in two ways:
# 1. irregularity to the signal, so instead of being one signal firing exactly a specific frequency, it fires only at approximately that frequency every cycle, plus or minus.
hz_sd = 1 
# 2. random firing.
signal_noise_ratio = 1

n_expected_spikes = minutes * 60 * hz
random_hz = rnorm(2 * n_expected_spikes, mean = hz, sd = hz_sd)
spike_locations = obs_per_second / random_hz
spike_locations = round(cumsum(spike_locations))
spike_locations = spike_locations[spike_locations < total_obs]

signal = rep(0, times = total_obs)
signal[spike_locations] = 1

x = signal

# The ones are so sparse I'm not going to worry about converting those to 0's
if(FALSE){

noise = sample.int(n, size = signal_noise_ratio * total_ones)
x[noise] = 1
}

a = acf(x, lag.max = 400, plot = FALSE)
a = a$acf[-1]
plot(a, type = "l")

