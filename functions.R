rel_diff = function(reference, other) abs((reference - other) / reference)


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


plot_acf = function(d, lag.max = 400, ...){
    # Drop the first observation that has ACF 1.
    a = acf(d$spike, lag.max = lag.max, plot = FALSE)
    a = a$acf[-1]
	x = seq_along(a)
	x = c(rev(-x), 0, x)
	y = c(rev(a), 0, a)
    plot(x, y
        , type = "l"
        , xlab = "lag"
        , ylab = "autocorrelation"
        , ...
    )
    a
}


# Omega determines phase.
# Jeff said 6-10 Hz is what we're looking for.
target_hz = 8
obs_per_sec = 1000
default_omega = target_hz * 2 * pi * obs_per_sec

optim_theta_index = function(ac
    , par = c(a = 1, b = 1, c = 1, omega = default_omega, tau1 = 1, tau2 = 1)
    , ...)
{
    t = seq_along(ac)
    fn = function(par){
        a = par["a"]
        b = par["b"]
        c = par["c"]
        omega = par["omega"]
        tau1 = par["tau1"]
        tau2 = par["tau2"]

        sin_term = a * sin(pi/2 - omega * t)
        exp1 = exp(-abs(t) / tau1)
        exp2 = c * exp(-t^2 / tau2)
        ac_star = (a * sin_term + b) * exp1 + exp2
        delta = (ac - ac_star)^2
        sum(delta)
    }
    out = optim(par, fn, ...)
    out$hertz = out[["par"]][["omega"]] / (2 * pi * obs_per_sec)
    out
}


theta_index = function(optim_results)
{
    p = optim_results[["par"]]
    p["a"] / p["b"]
}


# What is distinguishable beyond noise?
# This is just the plot.acf method
acf_clim = function (x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL,
    ylim = NULL, main = NULL, ci.col = "blue", ci.type = c("white",
        "ma"), max.mfrow = 6, ask = Npgs > 1 && dev.interactive(),
    mar = if (nser > 2) c(3, 2, 2, 0.8) else par("mar"), oma = if (nser >
        2) c(1, 1.2, 1, 1) else par("oma"), mgp = if (nser >
        2) c(1.5, 0.6, 0) else par("mgp"), xpd = par("xpd"),
    cex.main = if (nser > 2) 1 else par("cex.main"), verbose = getOption("verbose"),
    ...)
{
    ci.type <- match.arg(ci.type)
    if ((nser <- ncol(x$lag)) < 1L)
        stop("x$lag must have at least 1 column")
    if (is.null(ylab))
        ylab <- switch(x$type, correlation = "ACF", covariance = "ACF (cov)",
            partial = "Partial ACF")
    if (is.null(snames <- x$snames))
        snames <- paste("Series ", if (nser == 1L)
            x$series
        else 1L:nser)
    with.ci <- ci > 0 && x$type != "covariance"
    with.ci.ma <- with.ci && ci.type == "ma" && x$type == "correlation"
    if (with.ci.ma && x$lag[1L, 1L, 1L] != 0L) {
        warning("can use ci.type=\"ma\" only if first lag is 0")
        with.ci.ma <- FALSE
    }
    clim0 <- if (with.ci)
        qnorm((1 + ci)/2)/sqrt(x$n.used)
    else c(0, 0)
    Npgs <- 1L
    nr <- nser
    if (nser > 1L) {
        sn.abbr <- if (nser > 2L)
            abbreviate(snames)
        else snames
        if (nser > max.mfrow) {
            Npgs <- ceiling(nser/max.mfrow)
            nr <- ceiling(nser/Npgs)
        }
        opar <- par(mfrow = rep(nr, 2L), mar = mar, oma = oma,
            mgp = mgp, ask = ask, xpd = xpd, cex.main = cex.main)
        on.exit(par(opar))
        if (verbose) {
            message("par(*) : ", appendLF = FALSE, domain = NA)
            str(par("mfrow", "cex", "cex.main", "cex.axis", "cex.lab",
                "cex.sub"))
        }
    }
    if (is.null(ylim)) {
        ylim <- range(x$acf[, 1L:nser, 1L:nser], na.rm = TRUE)
        if (with.ci)
            ylim <- range(c(-clim0, clim0, ylim))
        if (with.ci.ma) {
            for (i in 1L:nser) {
                clim <- clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1,
                  i, i]^2)))
                ylim <- range(c(-clim, clim, ylim))
            }
        }
    }
    for (I in 1L:Npgs) for (J in 1L:Npgs) {
        dev.hold()
        iind <- (I - 1) * nr + 1L:nr
        jind <- (J - 1) * nr + 1L:nr
        if (verbose)
            message(gettextf("Page [%d,%d]: i =%s; j =%s", I,
                J, paste(iind, collapse = ","), paste(jind, collapse = ",")),
                domain = NA)
        for (i in iind) for (j in jind) if (max(i, j) > nser) {
            frame()
            box(col = "light gray")
        }
        else {
            clim <- if (with.ci.ma && i == j)
                clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1, i, j]^2)))
            else clim0
            plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab,
                ylab = if (j == 1)
                  ylab
                else "", ylim = ylim, ...)
            abline(h = 0)
            if (with.ci && ci.type == "white")
                abline(h = c(clim, -clim), col = ci.col, lty = 2)
            else if (with.ci.ma && i == j) {
                clim <- clim[-length(clim)]
                lines(x$lag[-1, i, j], clim, col = ci.col, lty = 2)
                lines(x$lag[-1, i, j], -clim, col = ci.col, lty = 2)
            }
            title(if (!is.null(main))
                main
            else if (i == j)
                snames[i]
            else paste(sn.abbr[i], "&", sn.abbr[j]), line = if (nser >
                2)
                1
            else 2)
        }
        if (Npgs > 1) {
            mtext(paste("[", I, ",", J, "]"), side = 1, line = -0.2,
                adj = 1, col = "dark gray", cex = 1, outer = TRUE)
        }
        dev.flush()
    }
    clim
}


# From http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html
plot.frequency.spectrum <- function(X.k, xlimits=c(0,500)) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))

  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}
