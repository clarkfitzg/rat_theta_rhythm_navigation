

rel_diff = function(reference, other) abs((reference - other) / reference)


# Read a single spike train file from a zip file that contains many such files.
read_one_spike_train_file = function(fname, zipfile = "data.zip"
    , col.names = c("index", "time", "spike"), ...)
{
    f = unz(zipfile, fname)
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


acf_from_spike_file = function(d, lag.max = 800)
{
    a = acf(d$spike, lag.max = lag.max, plot = FALSE)
    a$acf[-1]
}


file_name_and_theta_index = function(fname, ...)
{
    s = read_one_spike_train_file(fname, ...)
    a = acf_from_spike_file(s)
    ti = theta_index(a)
    data.frame(file = fname, theta_index = ti)
}


plot_acf = function(d, a = acf_from_spike_file(d), ...)
{
    # Drop the first observation that has ACF 1.
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


plot_acf_pair = function(file_pre, file_post, zipfile = "data.zip", plot_func = plot_acf, ...)
{

    d_pre = read_one_spike_train_file(file_pre, zipfile)
    d_post = read_one_spike_train_file(file_post, zipfile)

par(mfrow = c(2, 1))
    plot_func(d_pre, main = "Pre Injection", sub = file_pre, ...)
    plot_func(d_post, main = "Post Injection", sub = file_post, ...)
    list(pre = d_pre, post = d_post)
}


plot_spike_locs = function(d, ...)
{
    plot(which(as.logical(d$spike)), ...)
}




# The theta index should be small if there's no periodic signal.
# @param ac autocorrelation
# @param drop_first The first few observations have negative or lower autocorrelation, because cells can't fire consecutively. Drop these so they don't interfere with the fitting
theta_index = function(ac, tau1_range = c(0, 1e3), tau2_range = c(0, 1e5)
        , obs_per_sec = 1000, target_hz = 8, omega_range = c(0.1, 10) * target_hz * 2 * pi / obs_per_sec
        , drop_first = 20
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

    omega = optimize(fa, omega_range)$minimum
    tau1 = optimize(fb, tau1_range)$minimum
    tau2 = optimize(fc, tau2_range)$minimum

    bterm = exp(-time / tau1)
    aterm = sin(pi/2 - omega*time) * bterm
    cterm = exp(-time^2 / tau2)

    # Take out the intercept to be consistent with Tsanov.
    fit = lm(ac ~ 0 + aterm + bterm + cterm)
    cf = coef(fit)
    # Take the absolute value because there's nothing preventing the a term from being negative, and the b term may also be negative depending on what c does.
    abs(cf["aterm"] / cf["bterm"])
}


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


theta_index_optim = function(optim_results)
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
