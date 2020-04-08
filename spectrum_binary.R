# What does the ACF and periodogram look like if there actually *is* a signal?

# It might be the case that when we add some irregularity to the signal, the ACF just doesn't really work anymore.

set.seed(3480)

minutes = 8
obs_per_second = 1000
hz = 8
add_random_noise = TRUE
frac_signal_to_keep = 0.9

total_obs = minutes * 60 * obs_per_second
# Hmm... 

# Add noise in the following ways:
# 1. irregularity to the signal, so instead of being one signal firing exactly a specific frequency, it fires only at approximately that frequency every cycle, plus or minus.
hz_sd = 1 
# 2. random firing.
signal_noise_ratio = 1
# 3. drop some of the firing.

n_expected_spikes = minutes * 60 * hz
random_hz = rnorm(2 * n_expected_spikes, mean = hz, sd = hz_sd)
spike_locations = obs_per_second / random_hz
spike_locations = round(cumsum(spike_locations))
spike_locations = spike_locations[spike_locations < total_obs]
spikes_to_keep = round(frac_signal_to_keep * length(spike_locations))
spike_locations = sample(spike_locations, size = spikes_to_keep)

signal = rep(0, times = total_obs)
signal[spike_locations] = 1

x = signal

if(add_random_noise){
    p_noise = signal_noise_ratio * hz / obs_per_second
    noise_locations = sample(c(TRUE, FALSE), size = total_obs, replace = TRUE, prob = c(p_noise, 1 - p_noise))
    x = abs(x - noise_locations)
}

if(FALSE){

a = acf(x, lag.max = 800, plot = FALSE)
a = a$acf[-1]

plot(which(as.logical(x)))

spectrum(x, spans = 7, main = "7")

spectrum(x, spans = 51, main = "51")

png("spectrum_and_acf_with_signal.png", height = 2000, width = 1200)
par(mfrow = c(2, 1))
plot(a, type = "l")
spectrum(x, spans = 201, main = "201")
dev.off()

s = spectrum(x, spans = 501, main = "501")

spectrum(x, spans = 1001, main = "1001")

# We can definitely observe this signal in the plot of the spectrum.
recovered = obs_per_second * s$freq[which.max(s$spec)]
# And we can approximately recover the frequency of the original signal.

}
