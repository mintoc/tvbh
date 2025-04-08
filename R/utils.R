##----------------------------
## Functions to fit
## time-varying Beverton-Holt via Extended Kalman Filter, TMB, BSSM
## CM: 23/4/2023 
##----------------------------
library(TMB)

## eKF code

Zf <- function(x, fec, s0){
    ##x <- x
    q <- exp(x[1])
    p <- exp(x[2])
    yhat <- -log(p/q * (exp(q) - 1) + exp(q)/(fec * s0))
    return(yhat)
}

## from
## lnr_expr <- expression(-log(exp(lnp)/exp(lnq) * (exp(exp(lnq)) - 1) + exp(exp(lnq))/(fec * s)))
## D(lnr_expr, "lnp")
## D(lnr_expr, "lnq")

ddlnp <- function(lnp, lnq, fec, s0){
    lnp <- lnp
    lnq <- lnq
    res <- -(exp(lnp)/exp(lnq) * (exp(exp(lnq)) - 1)/(exp(lnp)/exp(lnq) * 
                                                      (exp(exp(lnq)) - 1) + exp(exp(lnq))/(fec * s0)))
    rm(list = c("lnp", "lnq"))
    return(res)
}

ddlnq <- function(lnp, lnq, fec, s0){
    lnp <- lnp
    lnq <- lnq
    res <- -((exp(lnp)/exp(lnq) * (exp(exp(lnq)) * exp(lnq)) - exp(lnp) * 
              exp(lnq)/exp(lnq)^2 * (exp(exp(lnq)) - 1) + exp(exp(lnq)) * 
              exp(lnq)/(fec * s0))/(exp(lnp)/exp(lnq) * (exp(exp(lnq)) - 1) +
                                    exp(exp(lnq))/(fec * s0)))
    rm(list = c("lnp", "lnq"))
    return(res)
}

Zprimef <- function(x, fec, s0){
    lnq <- x[1]
    lnp <- x[2]
    ddlnq0 <- ddlnq(lnp, lnq, fec, s0)
    ddlnp0 <- ddlnp(lnp, lnq, fec, s0)
    return(matrix(c(as.numeric(ddlnq0), as.numeric(ddlnp0)), nrow = 1))
}

bh_ekf <- function(mod, fec, data){
    y <- c(NA_real_, log(data$Recruits)) ## easier for initial state and subsequent indexing
    s <- c(NA_real_, data$SSB) ## aligned with r
    n <- length(y)
    H <- mod$H
    T <- mod$T; R <- mod$R; Q <- mod$Q
    a0 <- mod$a0
    P0 <- mod$P0
    ## dimension of state vector
    m <- length(a0)
    ## predicted state a at time t given data to time t-1 and
    ## updated state at time t given data to time t
    attm1 <- at <- matrix(0, nc = n, nr = m)
    F <- matrix(0, nc = n, nr = 1)
    Pttm1 <- Pt <- array(0, dim = c(m, m, n))
    ## initial state vector
    at[,1] <- a0
    Pt[,,1] <- P0
    ## innovation residuals
    yresid <- rep(NA_real_, length(y))
    ## log-likelihood
    ll <- rep(0,n)
    for(i in 2:n){
        ## prediction
        attm1[,i] <- T %*% at[,i-1]###
        Pttm1[,,i] <- T %*% Pt[,,i-1] %*% t(T) + Q
        ## update
        Z1 <- Zprimef(attm1[,i], fec, s[i])
        ##ytilde <- y[i] - Zf(attm1[,i]) + Zprimef(attm1[,i]) %*% attm1[,i]
        ## innovation residual
        yresid[i] <- y[i] - Zf(attm1[,i], fec, s[i])
        ## innovation/residual covariance
        F[i] <- Z1 %*% Pttm1[,,i] %*% t(Z1) + H
        ## Kalman gain
        K <- Pttm1[,,i] %*% t(Z1) %*% solve(F[i])
        ## 
        ## updated state
        at[,i] <- attm1[,i] +  K %*% yresid[i]
        Pt[,,i] <- Pttm1[,,i] - K %*% Z1 %*% Pttm1[,,i]
        ## F is scalar so no need for determinant here
        ll[i] <- -1/2 * log(F[i]) - 1/2 * t(yresid[i]) %*% solve(F[i]) %*% (yresid[i])
    }
    ## negative log-likelihood
    nll <- -sum(ll)
    res <- list(at = at, Pt = Pt, attm1 = attm1, F = F, Pttm1 = Pttm1, yresid = yresid, nll = nll)
    return(res)
}

## negative log likelihood
nll_ekf <- function(theta, fec, data, mod){
    mod1 <- build(theta, mod)
    filt <- bh_ekf(mod = mod1, fec, data)
    nll <- filt$nll
    rm(list = c("filt", "mod1"))
    return(nll)
}

## FIXED INTERVAL SMOOTHER
my_ks <- function(mod, at, Pt, Pttm1){
    ## section 3.6.2 of Harvey
    n <- length(at[1,])
    T <- mod$T;
    ## dimension of state vector
    m <- length(at[,1])
    ## predicted state a at time t given data to time t-1 and
    ## updated state at time t given data to time t
    atT <- matrix(0, nc = n, nr = m)
    PtT <- array(0, dim = c(m, m, n))
    ## last smoothed state is filtered state
    atT[,n] <- at[,n]
    PtT[,,n] <- Pt[,,n]
    for(i in (n-1):1){
        Pstar <- Pt[,,i] %*% t(T) %*% solve(Pttm1[,,i+1])
        atT[,i] <- at[,i] + Pstar %*% (atT[,i+1] - T %*% at[,i])
        PtT[,,i] <- Pt[,,i] + Pstar %*% (PtT[,,i+1] - Pttm1[,,i+1]) %*% t(Pstar)
    }
    res <- list(atT = atT, PtT = PtT)
    return(res)
}

## TMB CODE
compile("tvbh.cpp")
dyn.load(dynlib("tvbh"))

fit_bh_tv <- function(data, fec, tau){
    s <- data$SSB
    r <- data$Recruits
    year <- data$Year
    n <- length(year)
    starts <- get_start(data, fec = fec, tau = tau)
    ##print(starts["fec0"])
    ##fec <- as.numeric(10 * 1/b_bh)
    tau <- 1
    ##
    obj <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec),
        parameters = list(
            lnq = rep(starts["lnq0"], n),
            lnsdq = log(starts["sd0"]),
            lnp = rep(starts["lnp0"], n),
            lnsdp = log(starts["sd0"]),
            lnsde = log(starts["sd0"])
        ),
        ##map = list(lnsdq = factor(NA), lnsdrs = factor(NA)),
        random = c("lnq", "lnp"),
        DLL = "tvbh",
        silent = TRUE)
    ##
    opt <- try(  
        optim(par = obj$par,
              fn = obj$fn,
              gr = obj$gr,
              method = "BFGS",
              hessian = TRUE),
        silent = TRUE)
    if(class(opt) != "try-error"){
        rep <- sdreport(obj)
        srep <- summary(rep)
        rm(rep)
        lnq <- as.data.frame(srep[rownames(srep) == "lnq", ])
        names(lnq)[names(lnq) == "Std. Error"] <- "SE"
        names(lnq)[names(lnq) == "Estimate"] <- "est"
        qhat <- data.frame(qhat_tmb = exp(lnq$est))
        qhat$qlwr_tmb <- exp(lnq$est - 2 * lnq$SE)
        qhat$qupr_tmb <- exp(lnq$est + 2 * lnq$SE)
        qhat$year <- year
        ## HERE!
        lnp <- as.data.frame(srep[rownames(srep) == "lnp", ])
        names(lnp)[names(lnp) == "Std. Error"] <- "SE"
        names(lnp)[names(lnp) == "Estimate"] <- "est"
        phat <- data.frame(phat_tmb = exp(lnp$est))
        phat$plwr_tmb <- exp(lnp$est - 2 * lnp$SE)
        phat$pupr_tmb <- exp(lnp$est + 2 * lnp$SE)
        phat$year <- year
    }else{
        qhat <- data.frame(qhat_tmb = NA_real_, qlwr_tmb = NA_real_, qupr_tmb = NA_real_, year = year)
        phat <- data.frame(phat_tmb = NA_real_, plwr_tmb = NA_real_, pupr_tmb = NA_real_, year = year)
        opt <- NA_real_
        srep <- NA_real_
    }
    return(list(opt = opt, qhat = qhat, phat = phat, srep = srep))
}

## q-varying only version
compile("tvbh_qvary.cpp")
dyn.load(dynlib("tvbh_qvary"))

fit_bh_tv_q_only <- function(data, fec, tau){
    s <- data$SSB
    r <- data$Recruits
    year <- data$Year
    n <- length(year)
    starts <- get_start(data, fec = fec, tau)
    tau <- 1
    ##
    obj <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec),
        parameters = list(
            lnq = rep(starts["lnq0"], n),
            lnsdq = log(0.1),
            lnp = starts["lnp0"],
            lnsde = log(starts["sd0"])
        ),
        random = c("lnq"),
        DLL = "tvbh_q_only",
        silent = TRUE)
    ##
    opt <- try(optim(par = obj$par,
                     fn = obj$fn,
                     gr = obj$gr,
                     method = "BFGS",
                     hessian = TRUE),
               silent = TRUE)
    if(class(opt) != "try-error"){
        rep <- sdreport(obj)
        srep <- summary(rep)
        rm(rep)
        lnq <- as.data.frame(srep[rownames(srep) == "lnq", ])
        names(lnq)[names(lnq) == "Std. Error"] <- "SE"
        names(lnq)[names(lnq) == "Estimate"] <- "est"
        qhat <- data.frame(qhat_tmb = exp(lnq$est))
        qhat$qlwr_tmb <- exp(lnq$est - 2 * lnq$SE)
        qhat$qupr_tmb <- exp(lnq$est + 2 * lnq$SE)
        qhat$year <- year
        ## HERE!
        lnp <- as.data.frame(t(srep[rownames(srep) == "lnp", ]))
        names(lnp)[names(lnp) == "Std. Error"] <- "SE"
        names(lnp)[names(lnp) == "Estimate"] <- "est"
        phat <- data.frame(phat_tmb = exp(lnp$est))
        phat$plwr_tmb <- exp(lnp$est - 2 * lnp$SE)
        phat$pupr_tmb <- exp(lnp$est + 2 * lnp$SE)
        phat <- cbind(phat, year)
    }else{
        qhat <- data.frame(qhat_tmb = NA_real_, qlwr_tmb = NA_real_, qupr_tmb = NA_real_, year = year)
        phat <- data.frame(phat_tmb = NA_real_, plwr_tmb = NA_real_, pupr_tmb = NA_real_, year = year)
        opt <- NA_real_
        srep <- NA_real_
    }
    return(list(opt = opt, qhat = qhat, phat = phat, srep = srep))
}


## p-varying only version
compile("tvbh_pvary.cpp")
dyn.load(dynlib("tvbh_pvary"))

fit_bh_tv_p_only <- function(data, fec, tau){
    s <- data$SSB
    r <- data$Recruits
    year <- data$Year
    n <- length(year)
    starts <- get_start(data, fec = fec, tau = tau)
    tau <- 1
    ##
    obj <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec),
        parameters = list(
            lnq = starts["lnq0"],
            lnp = rep(starts["lnq0"], n),
            lnsdp = log(0.1),
            lnsde = log(starts["sd0"])
        ),
        random = c("lnp"),
        DLL = "tvbh_p_only",
        silent = TRUE)
    ##
    opt <- try(optim(par = obj$par,
                     fn = obj$fn,
                     gr = obj$gr,
                     method = "BFGS",
                     hessian = TRUE),
               silent = TRUE)
    if(class(opt) != "try-error"){
        rep <- sdreport(obj)
        srep <- summary(rep)
        rm(rep)
        lnq <- as.data.frame(t(srep[rownames(srep) == "lnq", ]))
        names(lnq)[names(lnq) == "Std. Error"] <- "SE"
        names(lnq)[names(lnq) == "Estimate"] <- "est"
        qhat <- data.frame(qhat_tmb = exp(lnq$est))
        qhat$qlwr_tmb <- exp(lnq$est - 2 * lnq$SE)
        qhat$qupr_tmb <- exp(lnq$est + 2 * lnq$SE)
        qhat <- cbind(qhat, year)
        ## HERE!
        lnp <- as.data.frame(srep[rownames(srep) == "lnp", ])
        names(lnp)[names(lnp) == "Std. Error"] <- "SE"
        names(lnp)[names(lnp) == "Estimate"] <- "est"
        phat <- data.frame(phat_tmb = exp(lnp$est))
        phat$plwr_tmb <- exp(lnp$est - 2 * lnp$SE)
        phat$pupr_tmb <- exp(lnp$est + 2 * lnp$SE)
        phat <- cbind(phat, year)
        phat$year <- year        
    }else{
        qhat <- data.frame(qhat_tmb = NA_real_, qlwr_tmb = NA_real_, qupr_tmb = NA_real_, year = year)
        phat <- data.frame(phat_tmb = NA_real_, plwr_tmb = NA_real_, pupr_tmb = NA_real_, year = year)
        opt <- NA_real_
        srep <- NA_real_
    }
    return(list(opt = opt, qhat = qhat, phat = phat, srep = srep))
}

## static version
compile("bh_static.cpp")
dyn.load(dynlib("bh_static"))


## STARTING VALUE CODE
get_start <- function(data, fec, tau){
    ## tau is the age at recruitment
    r <- data$Recruits
    s <- data$SSB
    y <- r/s
    f0 <- try(glm(y ~ s, family = Gamma), silent = TRUE)
    if(class(f0)[1] != "try-error"){
        a_bh <- coef(f0)[2]
        b_bh <- coef(f0)[1]
    }else{
        ## make they're negative to enter heuristic
        a_bh <- b_bh <- -1
    }
    ## ensure they are positive
    if(a_bh < 0 | b_bh < 0){
        a_start <- as.numeric(1/quantile(r, 0.5))
        Smax <- max(s, na.rm = TRUE)
        Smin <- min(s, na.rm = TRUE)
        contrast <- (Smax - Smin)/Smax
        ## say it reaches 90% of maximal recruitment at median ssb -
        b_start <- median(s) * (a_start/0.9 - a_start)
        a_bh <- a_start
        b_bh <- b_start
    }
    q_start <- log(fec * b_bh)/tau
    p_start <- a_bh * q_start / (exp(q_start * tau) - 1)
    ## starting value for standard deviations
    y <- log(r/s)
    ## split the variance
    f0 <- lm(y ~ s)
    sd0 <- sqrt((summary(f0)$sigma^2)/3) ## split three ways - think about this
    ##
    res <- c(a0 = as.numeric(a_bh),
             b0 = as.numeric(b_bh),
             lnq0 = as.numeric(log(q_start)),
             lnp0 = as.numeric(log(p_start)),
             fec0 = as.numeric(fec),
             sd0 = sd0)
    return(res)
}

## function to calculate the deviance for bssm model fit
get_deviance_bssm_mcmc <- function(model, fit, particles = 10){
    posterior <- fit$posterior
    theta <- fit$theta
    ll <- sapply(1:nrow(theta), function(z){
        model$theta[] <- theta[z,]
        ll_smooth <-  particle_smoother(model, particles = particles, method = "psi")$logLik
        mcmc_static$posterior[z] - dnorm(exp(par$log_H.sd0), 0, 2, 1) - ll
    })
    return(list(D = -2 * ll))
}

## function to calculate the deviance information criterion for bssm
DIC <- function(fit, static = TRUE){
    ## extract the marginal log likelihood
    log_posterior <- fit$posterior
    if(static){
        log_prior <- apply(fit$theta, 1, function(z){log_prior_pdf0(z)})
    }else{
        log_prior <- apply(fit$theta, 1, function(z){log_prior_pdf12(z)})
    }
    log_like_short <- log_posterior - log_prior
    log_like_long <- rep(log_like_short, times = fit$counts)
    ## weights
    wts <- rep(fit$weights[, 1], times = fit$counts)
    idx <- which(wts > 0)
    wts <- wts[idx]
    log_like_long <- log_like_long[idx]
    D <- -2 * log_like_long
    pD <- 1/2 * weighted_var(D, wts)
    dic_static <- weighted_mean(D, wts) + pD
    return(dic_static)
}
