# Last modified:
# Mon Jul 13 14:47:15 PDT 2020

# Jeff Calton and Clark Fitzgerald, CSUS
#
# This script analyzes neural spike data to determine if the neural activity
# shows any rhythmic modulation at a particular frequency, 6-12 Hz.
#
#
# Call it from the command line as follows:
#
#   $ Rscript theta_index.R data_directory
#
#
# For each file in the data directory, the script produces the following:
#
#   1. Computes the theta index, and appends it to the file theta_index.csv in the current directory.
#       This is based on (Tsanov et. al. 2011 p. 9492, eq. 1)
#   2. Saves a pdf plot containing the acf and the spectrum in the current directory.
#       The file name is based on the name of the data file.


# Parameters
############################################################

# The frequency we're most interested in.
TARGET_HZ = 8
# The ranges of frequencies that we'll consider for fitting.
HZ_RANGE = TARGET_HZ * c(0.1, 10)
MAX_PLOT_HZ = 50

# The largest lag for the autocorrelation.
LAG = 800

# Degree of smoothness of the periodogram, see examples in spectrum()
# Larger is more smooth.
SPANS = 201

# The number of observations per second
OBS_PER_SECOND = 1000

# The autocorrelation for small lags is low or negative because neurons can't fire right after one another.
# Drop these.
DROP_FIRST = 20

# Height and width in pixels of the png plots produced.
PNG_HEIGHT = 2000
PNG_WIDTH = 1200
FONT_SIZE = 30

# Assume that any file in the data directory with this suffix is a data file.
DATAFILE_SUFFIX = "SpikeSortedForAutocorrelation"
PLOTFILE_SUFFIX = "_plot.png"
ACF_SUFFIX = "_acf.csv"
PERIODOGRAM_SUFFIX = "_periodogram.csv"

# CSV file to save the actual theta index values
THETA_INDEX_FILE = "theta_index.csv"


# Functions
############################################################


process_one_file = function(fname, ...)
{
    message("processing ", fname)
    s = read_one_spike_train_file(fname, ...)
    a = lag_counts_from_spike_file(s)
    ti = theta_index(a)

    bfname = basename(fname)

    # Plotting
    png(paste0(bfname, PLOTFILE_SUFFIX), height = PNG_HEIGHT, width = PNG_WIDTH
        , pointsize = FONT_SIZE)
    par(mfrow = c(2, 1))
    plot_lag_counts(a = a
        , main = paste0("Counts", bfname)
        , sub = sprintf("Theta index = %g4", ti)
        )
    sp = spectrum(s[, "spike"], main = paste0("Periodogram with spans = ", SPANS), spans = SPANS)
    dev.off()

    # Save ACF and periodogram into CSV files
    adata = data.frame(lag = seq_along(a), counts = a)
    write.csv(adata, file = paste0(bfname, ACF_SUFFIX), row.names = FALSE)

    frequency = OBS_PER_SECOND * sp[["freq"]]
    pdata = data.frame(frequency, spectrum = sp[["spec"]])
    pdata = pdata[frequency < MAX_PLOT_HZ, ]
    write.csv(pdata, file = paste0(bfname, PERIODOGRAM_SUFFIX), row.names = FALSE)

    ti = data.frame(theta_index = ti, file = bfname)
    write.table(ti, file = THETA_INDEX_FILE, append = TRUE
                , row.names = FALSE, col.names = FALSE, sep = ",")
}


lag_counts_from_spike_file = function(d, lag.max = LAG)
{
    s = d[, "spike"]
    a = acf(s, lag.max = lag.max, plot = FALSE, type = "covariance", demean = FALSE)
    # Scale it into counts
    length(s) * a[["acf"]][-1]
}


read_one_spike_train_file = function(fname,
    col.names = c("index", "time", "spike"), ...)
{
    d = read.table(fname, header = FALSE, col.names = col.names, ...)
    
    time_deltas = table(diff(d[, "time"]))
    if(length(time_deltas) > 1)
        stop("Time intervals are not all the same.")

    n = dim(d)[1]
    # Assuming observations are at every millisecond.
    minutes = n / (60 * 1000)
    message(sprintf("Recording period was %g minutes.", minutes))

    d
}


plot_lag_counts = function(d, a = lag_counts_from_spike_file(d), ...)
{
    # Drop the first observation that has ACF 1.
	x = seq_along(a)
	x = c(rev(-x), 0, x)
	y = c(rev(a), 0, a)
    plot(x, y
        , type = "l"
        , xlab = "time (ms)"
        , ylab = "count"
        , ...
    )
    a
}


# The theta index should be small if there's no periodic signal.
# @param ac autocorrelation
# @param drop_first The first few observations have negative or lower autocorrelation, because cells can't fire consecutively. Drop these so they don't interfere with the fitting
theta_index = function(ac, tau1_range = c(0, 1e3), tau2_range = c(0, 1e5)
        , obs_per_sec = OBS_PER_SECOND
        , omega_range = HZ_RANGE * 2 * pi / obs_per_sec
        , drop_first = DROP_FIRST
){
    ac = ac[-seq(drop_first)]
    time = seq_along(ac)

    fa = function(omega){
        aterm = sin(pi/2 - omega*time)
        fit = lm(ac ~ aterm)
        sum(residuals(fit)^2)
    }

    fb = function(tau1){
        bterm = exp(-time / tau1)
        fit = lm(ac ~ bterm)
        sum(residuals(fit)^2)
    }

    fc = function(tau2){
        cterm = exp(-time^2 / tau2)
        fit = lm(ac ~ cterm)
        sum(residuals(fit)^2)
    }

    omega = optimize(fa, omega_range)[["minimum"]]
    tau1 = optimize(fb, tau1_range)[["minimum"]]
    tau2 = optimize(fc, tau2_range)[["minimum"]]

    bterm = exp(-time / tau1)
    aterm = sin(pi/2 - omega*time) * bterm
    cterm = exp(-time^2 / tau2)

    # Take out the intercept to be consistent with Tsanov.
    fit = lm(ac ~ 0 + aterm + bterm + cterm)

    # The aterm is not significant, even when in a model all by itself.
    # In other words, we can't tell that a is not 0.
    #fit = lm(ac ~ aterm)

    print(summary(fit))
    cf = coef(fit)
    # Take the absolute value because there's nothing preventing the a term from being negative, and the b term may also be negative depending on what c does.
    abs(cf["aterm"] / cf["bterm"])
}


main = function(dirname = ".")
{
    data_directory = dirname[1]
    files = list.files(data_directory, pattern = paste0("*", DATAFILE_SUFFIX, "$"))
    files = file.path(data_directory, files)

    if(!file.exists(THETA_INDEX_FILE)){
        writeLines("theta_index,file", THETA_INDEX_FILE)
    }

    lapply(files, function(x) tryCatch(process_one_file(x), error = function(e) message("ERROR processing ", x)))
}


# Run it!
#if(!interactive()){
if(TRUE){
    dirname = commandArgs(trailingOnly = TRUE)
    if(length(dirname) == 0){
        main()
    } else {
        main(dirname)
    }
}
#}
