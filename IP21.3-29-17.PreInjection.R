d = read.csv("IP21.3-29-17.PreInjection.Standard.Cylinder.2017-03-29_17-13-40.SE16.nse.converted.SpikeSortedForAutocorrelation.csv"
, header = FALSE
, col.names = c("index", "time", "fire"))

n = dim(d)[1]

# If these observations are every milliseconds then the data spans 8 minutes:
minutes = n / (60 * 1000)

# Verify that measurements are at regular timestamps:
table(diff(d$time))
# Yes.

a = acf(d$fire, lag.max = 800)

a2 = a$acf[-1]

# We see some negative correlation for the first 5 periods, consistent with the neuron not firing twice within 5 ms.
# After that it's a decaying positive correlation.
pdf("acf800.pdf")
plot(seq(-799, 800), c(rev(a2), a2)
    , type = "l"
    , xlab = "lag"
    , ylab = "autocorrelation"
)
dev.off()

plot(a2[1:60], type = "l")

# Not showing much difference.
plot(a2[1:400], type = "l")
