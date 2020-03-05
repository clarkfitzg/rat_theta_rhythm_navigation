# Read a single spike train file from a zip file that contains many such files.
read_one_spike_train_file = function(fname, zipfile = "data.zip"
    , col.names = c("index", "time", "spike"), ...)
{
    f = unz("data.zip", fname)
    d = read.table(f, header = FALSE, col.names = col.names, ...)
    
    time_deltas = table(diff(d$time))
    if(length(time_deltas) > 1)
        stop("Time intervals are not all the same.")

    n = dim(d)[1]
    # Assuming observations are at every millisecond.
    minutes = n / (60 * 1000)
    message(sprintf("Recording period was %g minutes.", minutes))

    d
}


plot_acf = function(d, lag.max = 100, ...){
    # Drop the first observation that has ACF 1.
    a = acf(d$spike, lag.max = lag.max, plot = FALSE)
    a = a$acf[-1]
    plot(a
        , type = "l"
        , xlab = "lag"
        , ylab = "autocorrelation"
        , ...
    )
}


# Compare ACF before and after injection

PreInjection = read_one_spike_train_file("IP18.12-9-16.PreInjection.Standard.Cylinder.2016-12-09_14-31-21.SE9.nse.converted.SpikeSortedForAutocorrelation")

PostInjection = read_one_spike_train_file("IP18.12-9-16.PostInjection.Standard.Cylinder.2016-12-09_19-18-01.SE15.nse.converted.SpikeSortedForAutocorrelation")


# The pre injection ACF is around 0.002, while the post injection ACF is around 0.01, around five times larger.
# What could that mean?
par(mfrow = c(2, 1))
plot_acf(PreInjection, main = "Pre Injection")
plot_acf(PostInjection, main = "Post Injection")
