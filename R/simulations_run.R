##--------------------------
## Simulated data and fits 
## CM: 9/9/24
##
##--------------------------

library(ggplot2); theme_set(theme_bw())
library(viridis)

##------------
## SIMULATION 
##------------

## steepness - forms sr typologies here
h <- c(0.4, 0.7, 0.9)

## 1/a is maximum recruitment
a <- 1/80

## S0 - value we want for unfished biomass
S0 <- 300

## define the SPRF0 to result in that unfished biomass
phi0 <- a * S0 / (1 - (1 - h)/(4 * h))

## define b 1/b is the slope at the origin
b <- (phi0 - h * phi0) / (4 *h)

## check S0 is correct
(phi0 - b)/a

## note only h, S0 and a define scenarios
scen <- data.frame(h, a, b, phi0, S0)

scen$tau <- 1

## expand grid for dataframes
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

## 50 iterations per scenario
niter <- 50
sim_df <- expand.grid.df(scen,
                         expand.grid(
                             fec_true = 100, ## fecundity
                             sde = c(0.2, 0.4),
                             ##sdlnsar = c(0, 0, 0.1, 0.2), ## sd of ln slope @ origin
                             ##sdlnR0 = c(0, 0, 0.1, 0.2), ## sd of ln max recruit
                             sdlnq = c(0, 0.1, 0.2),
                             sdlnp = c(0, 0.1, 0.2),
                             contrast = c("lower", "upper", "full"),
                             n = c(40, 80),
                             iter = 1:niter))

sim_df <- subset(sim_df, !(sdlnq > 0 & sdlnp > 0))

sim_df$ID <- 1:nrow(sim_df)

## list to hold the simulated series
sim_list <- list()

for(i in 1:nrow(sim_df)){
    print(i)
    pars <- sim_df[i,]
    ## fecundity
    fec <- pars$fec_true
    ## age at recruitment
    tau <- pars$tau
    ## main parameters
    a <- pars$a ## 1/a is maximal recruitment
    b <- pars$b ## 1/b is slope at origin
    ## mortality terms
    ## density-independent mortality
    q <- log(b * fec)/tau
    ## density-dependent mortality
    p <- a * q / (exp(q * tau) - 1)
    ## check
    ## p/q * (exp(q * tau) - 1) ## should be a
    ## exp(q * tau)/fec ## should be b
    ##
    n <- pars$n
    year <- 1:n    
    ##
    sde <- pars$sde
    ## process sds
    sdlnq <- pars$sdlnq
    sdlnp <- pars$sdlnp
    ##
    lnq0 <- log(q)    
    lnp0 <- log(p)
    ##
    ## process deviations
    eta_lnq <- rnorm(n, 0, sdlnq)
    eta_lnp <- rnorm(n, 0, sdlnp)
    ## time-varying pars
    lnqt <- lnq0 + cumsum(eta_lnq)
    lnpt <- lnp0 + cumsum(eta_lnp)
    ##
    qt <- exp(lnqt)
    pt <- exp(lnpt)
    ## ensure that the mean is q
    qt <- qt/mean(qt) * q
    pt <- pt/mean(pt) * p
    ##
    at <- pt/qt * (exp(qt * tau) - 1)
    bt <- exp(qt * tau)/fec
    ## max and min ssb from contrast
    contrast <- pars$contrast
    sprf0 <- pars$phi0    
    S0 <- (sprf0 - b)/a
    ##
    if(contrast == "lower"){        
        mins <- 0
        maxs <- S0/2
    }
    if(contrast == "upper"){
        mins <- S0/2
        maxs <- S0
    }
    if(contrast == "full"){
        mins <- 0
        maxs <- S0
    }
    ## uniform distribution on s
    s <- runif(n, mins, maxs)
    r <- rlnorm(n, meanlog = log(1/(at + bt/s)), sdlog = sde)
    ##
    dat <- data.frame(s, r, qt, pt, year = 1:n)
    rownames(pars) <- NULL
    dat <- cbind(dat, pars)
    ##timeseries <- rbind(timeseries, dat)
    sim_list[[i]] <- dat
    rm(dat)
}

timeseries <- do.call(rbind, sim_list)

library(viridis)

## visualise
pdf("../tex/figures/sr_sim_visualisation_new.pdf", height = 6, width = 9)
for(j in h){
    for(i in c(0.2, 0.4)){
        ##p <- ggplot(subset(timeseries, h == j & iter == 1 & sde == i & n == 40),
        p <- ggplot(subset(timeseries, h == j & iter == 2 & sde == i & n == 40),
                    aes(x = s, y = r, colour = year)) +
            geom_point() +
            ##geom_hline(data = qmeans, aes(yintercept = qmean)) +
            facet_grid(contrast ~ sdlnq + sdlnp, labeller = label_both) +
            scale_colour_viridis() +
            xlab("Spawner biomass") +
            ylab("Recruitment") +
            ggtitle(paste("steepness =", j, "; sde =", i)) +
            theme(legend.position = "bottom")
        print(p)
    }
}
dev.off()


##------------
## ESTIMATION
##------------
## load in the estimation code for eKF and TMB
source("utils.R")

## BSSM code
library(bssm)
library(dplyr)
library(diagis)

## create shared library for BSSM fit 
Rcpp::sourceCpp("bh_bssm_new.cpp")
pntrs <- create_xptrs()


## containers for results
m <- nrow(sim_df)

names_df <- expand.grid(c("q", "p"), c("hat_ekf", "lwr_ekf", "upr_ekf",
                                       "hat_tmb", "lwr_tmb", "upr_tmb",
                                       "hat_psiAPF", "lwr_psiAPF", "upr_psiAPF"))

names_df <- expand.grid(c("q", "p"), c("hat_ekf", "hat_tmb","hat_bssm"))

new_names <- apply(names_df, 1, paste, collapse = "")

new_names <- new_names[order(new_names)]

timeseries[, new_names] <- NA

## slots for information criteria
names <- expand.grid(method = c("ekf", "tmb"),
                     ic = c("aic", "bic"),
                     model = c("static", "qvary", "pvary"))

names <- rbind(names,
               expand.grid(method = "bssm", ic = "dic", model = c("static", "qvary", "pvary")))

sim_df[, apply(names, 1, paste, collapse = "_")] <- NA

## variables to remove at end of each iteration
vars2remove <- c("dat", "fec0", "starts", "mod", "build", "f0", "f1", "f2", "ic_idx", "best", "smooth0", "smooth1", "smooth2", "smooth", "initial_theta", "known_params", "known_tv_params", "model", "mcmc_static", "dic", "mcmc_qvary", "states1", "mcmc_pvary", "states2", "states", "s", "r", "year", "n", "tau", "obj0", "opt0", "srep", "srep_df0", "obj1", "opt1", "srep_df1", "obj2", "opt2", "srep_df2")

## track progress in a file
cat("", file = "progress.txt")

time0 <- Sys.time()

for(j in 1:m){
    print(j)
    ##
    idx <- which(timeseries$ID == j)
    dat <- timeseries[idx,]
    dat$SSB <- dat$s
    dat$Recruits <- dat$r
    dat$Year <- dat$year
    ## pretend fec known
    fec0 <- unique(dat$fec_true)
    tau <- dat$tau[1]
    starts <- get_start(dat, fec0, tau = tau)
    ##~~~~~~~~~~~~~
    ## eKF fits
    ##~~~~~~~~~~~~~
    ## with(dat, plot(SSB, r, xlim = c(0, max(SSB)), ylim = c(0, max(r))))
    ## curve(1/(starts["a0"] + starts["b0"]/x), add = TRUE)
    ##--------------
    ## STATIC MODEL
    ##--------------
    mod <- list(H = matrix(starts["sd0"]^2),
                T = diag(2),
                ##Q = diag(c(starts["sd0"]^2, 0)),
                Q = diag(rep(starts["sd0"]^2, 2)),
                a0 = c(starts["lnq0"], starts["lnp0"]),
                P0 = diag(2) * 1e7)
    ## build function
    ## no process variation
    build <- function(theta, mod){
        mod$Q[1,1] <- 0 
        mod$Q[2,2] <- 0 
        mod$H[1,1] <- exp(theta[1])^2
        return(mod)
    }
    ## fit
    f0 <- try(optim(
        par = rep(log(starts["sd0"]), 1),
        fn = nll_ekf,
        fec = fec0,
        mod = mod,
        data = dat,
        control = list(maxit = 1e3),
        method = "BFGS"),
        silent = TRUE)
    if(class(f0) != "try-error"){
        if(f0$convergence == 0){
            sim_df$ekf_aic_static[j] <- 2 * f0$value + 2 * 3
            sim_df$ekf_bic_static[j] <- 2 * f0$value + log(nrow(dat)) * 3
            ## smoothed states
            mod_hat <- build(f0$par, mod)
            filt <- bh_ekf(mod = mod_hat, fec0, data = dat)
            smooth0 <- my_ks(mod_hat, at = filt$at, Pt = filt$Pt, Pttm1 = filt$Pttm1)
        }
    }
    ##--------
    ## q-vary 
    ##--------
    build <- function(theta, mod){
        mod$Q[1,1] <- exp(theta[1])^2 ## so thetas are log sds
        mod$Q[2,2] <- 0 
        mod$H[1,1] <- exp(theta[2])^2 
        return(mod)
    }
    ##
    f1 <- try(optim(
        par = rep(log(starts["sd0"]), 2),
        fn = nll_ekf,
        fec = fec0,
        mod = mod,
        data = dat,
        control = list(maxit = 1e3),
        method = "BFGS"),
        silent = TRUE)
    ##
    if(class(f1) != "try-error" ){
        ## sometimes matrix is singular so calculating all pieces first
        ## smoothed states
        mod_hat <- build(f1$par, mod)
        filt <- bh_ekf(mod = mod_hat, fec0, data = dat)
        smooth1 <- try(my_ks(mod_hat, at = filt$at, Pt = filt$Pt, Pttm1 = filt$Pttm1), silent = TRUE)
        if(f1$convergence == 0 & class(smooth1) != "try-error"){
            sim_df$ekf_aic_qvary[j] <- 2 * f1$value + 2 * 4
            sim_df$ekf_bic_qvary[j] <- 2 * f1$value + log(nrow(dat)) * 4
        }
    }
    ##--------
    ## p-vary 
    ##--------
    build <- function(theta, mod){
        mod$Q[1,1] <- 0 
        mod$Q[2,2] <- exp(theta[1])^2 
        mod$H[1,1] <- exp(theta[2])^2 
        return(mod)
    }
    ##
    f2 <- try(optim(
        par = rep(log(starts["sd0"]), 2),
        fn = nll_ekf,
        fec = fec0,
        mod = mod,
        data = dat,
        control = list(maxit = 1e3),
        method = "BFGS"),
        silent = TRUE)
    ##
    if(class(f2) != "try-error" ){
        ## smoothed states
        mod_hat <- build(f2$par, mod)
        filt <- bh_ekf(mod = mod_hat, fec0, data = dat)
        smooth2 <- try(my_ks(mod_hat, at = filt$at, Pt = filt$Pt, Pttm1 = filt$Pttm1), silent = TRUE)
        ##
        if(f2$convergence == 0 & class(smooth2) != "try-error"){
            sim_df$ekf_aic_pvary[j] <- 2 * f2$value + 2 * 4
            sim_df$ekf_bic_pvary[j] <- 2 * f2$value + log(nrow(dat)) * 4
        }
    }
    ## get best ekf model
    ic_idx <- grep(paste0("^", "ekf_aic"), names(sim_df))
    best <- apply(sim_df[j, ic_idx], 1, function(z){
        if(all(is.na(z))){
            3
        }else{
            which.min(z)-1 ## as smooths start at 0
        }
    })
    ## and the best fitting states
    if(best == 3){
        ## all are unconverged - left as NA's
    }else{
        smooth <- get(paste0("smooth", best))
        if(class(smooth) == "try-error"){
            dat$qhat_ekf <- NA
            dat$phat_ekf <- NA
        }else{
            dat$qhat_ekf <- exp(smooth$atT[1,-1])
            dat$phat_ekf <- exp(smooth$atT[2,-1])
        }
    }
    ##~~~~~~~
    ## BSSM
    ##~~~~~~~
    ##-------------
    ## BSSM STATIC
    ##-------------
    ## starting at a maximum of 1.5 for sampler
    if(starts["sd0"] < 1.5){
        initial_theta <- c(log_H = log(starts["sd0"]))
    }else{
        initial_theta <- c(log_H = log(1.5))
    }
    known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                      P11 = 1e0, P22 = 1e0,
                      tau = tau,
                      fec = fec0)
    ##
    known_tv_params <- matrix(dat$SSB)
    ##
    model <- ssm_nlg(y = log(dat$r),
                     a1=pntrs$a1,
                     P1 = pntrs$P1, 
                     Z = pntrs$Z_fn,
                     H = pntrs$H_fn,
                     T = pntrs$T_fn,
                         R = pntrs$R0_fn, 
                     Z_gn = pntrs$Z_gn,
                     T_gn = pntrs$T_gn,
                     theta = initial_theta,
                     log_prior_pdf = pntrs$log_prior_pdf0_unif,
                     ## using uniform priors as half normal occassionally
                     ## goes to very low values on log scale that result in core
                     ## dump - happened with iteration 96 on qvary model
                     known_params = known_params,
                     known_tv_params = known_tv_params,
                     n_states = 2,
                     n_etas = 2,
                     state_names = c("logq", "logp"))
    ##
    mcmc_static <- try(run_mcmc(model, particles = 50, iter = 1e4, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(1) * 0.1), silent = TRUE)
    ##
    if(class(mcmc_static) != "try-error"){
        dic <- try(DIC(mcmc_static, static = TRUE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            sim_df$bssm_dic_static[j] <- dic
            states0 <- summary(mcmc_static, variable = "states")
        }
    }
    ##-------------
    ## BSSM q vary
    ##-------------
    if(starts["sd0"] < 1.5){
        initial_theta <- c(log_H = log(starts["sd0"]), log_R = log(starts["sd0"]))
    }else{
        initial_theta <- log(c(1.5, 0.5))
    }
    ##
    known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                      P11 = 1e0, P22 = 1e0,
                      tau = tau,
                      fec = fec0)
    ##
    model <- ssm_nlg(y = log(dat$r),
                     a1 = pntrs$a1,
                     P1 = pntrs$P1, 
                     Z = pntrs$Z_fn,
                     H = pntrs$H_fn,
                     T = pntrs$T_fn,
                     R = pntrs$R1_fn, 
                     Z_gn = pntrs$Z_gn,
                     T_gn = pntrs$T_gn,
                     theta = initial_theta,
                     log_prior_pdf = pntrs$log_prior_pdf12_unif,
                     known_params = known_params,
                     known_tv_params = known_tv_params,
                     n_states = 2,
                     n_etas = 2, 
                     state_names = c("logq", "logp"))
    ##
    mcmc_qvary <- try(run_mcmc(model, particles = 50, iter = 1e4, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(2) * 0.1), silent = TRUE)
    if(class(mcmc_qvary) != "try-error"){
        dic <- try(DIC(mcmc_qvary, static = FALSE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            sim_df$bssm_dic_qvary[j] <- dic
            states1 <- summary(mcmc_qvary, variable = "states")            
        }
    }
    ##-------------
    ## BSSM p vary
    ##-------------
    if(starts["sd0"] < 1.5){
        initial_theta <- c(log_H = log(starts["sd0"]), log_R = log(starts["sd0"]))
    }else{
        initial_theta <- log(c(1.5, 0.5))
    }
    known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                      ##P11 = 1e2, P22 = 1e2,
                      P11 = 1e0, P22 = 1e0,
                      tau = tau,
                      fec = fec0)
    ##
    model <- ssm_nlg(y = log(dat$r),
                     a1=pntrs$a1,
                     P1 = pntrs$P1, 
                     Z = pntrs$Z_fn,
                     H = pntrs$H_fn,
                     T = pntrs$T_fn,
                     R = pntrs$R2_fn, 
                     Z_gn = pntrs$Z_gn,
                     T_gn = pntrs$T_gn,
                     theta = initial_theta,
                     log_prior_pdf = pntrs$log_prior_pdf12_unif,
                     known_params = known_params,
                     known_tv_params = known_tv_params,
                     n_states = 2,
                     n_etas = 2,
                     state_names = c("logq", "logp"))
    ##
    mcmc_pvary <- try(run_mcmc(model, particles = 50, iter = 1e4, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(2) * 0.1), silent = TRUE)
    ##
    if(class(mcmc_pvary) != "try-error"){
        dic <- try(DIC(mcmc_pvary, static = FALSE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            sim_df$bssm_dic_pvary[j] <- dic
            states2 <- summary(mcmc_pvary, variable = "states")
        }
    }
    ## get best bssm model
    ic_idx <- grep(paste0("^", "bssm_dic"), names(sim_df))
    best <- apply(sim_df[j, ic_idx], 1, function(z){
        if(all(is.na(z))){
            3
        }else{
            which.min(z)-1 ## as smooths start at 0
        }
    })
    ## and the best fitting states
    if(best == 3){
        ## all are unconverged
    }else{
        states <- get(paste0("states", best))
        states <- subset(states, time != (dat$n[1] + 1))
        ## might be better of exponentiating this and then taking the mean
        ## do this later
        dat$qhat_bssm <- exp(subset(states, variable == "logq")$Mean)
        dat$phat_bssm <- exp(subset(states, variable == "logp")$Mean)
    }
    ##~~~~~~~
    ## TMB
    ##~~~~~~~
    s <- dat$SSB
    r <- dat$Recruits
    year <- dat$Year
    n <- length(year)
    ##------------
    ## TMB static
    ##------------
    obj0 <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec0),
        parameters = list(
            lnq = starts["lnq0"],
            lnp = starts["lnp0"],
            lnsde = log(starts["sd0"])
        ),
        DLL = "bh_static",
        silent = TRUE)
    ##
    opt0 <- try(  
        optim(par = obj0$par,
              fn = obj0$fn,
              gr = obj0$gr,
              method = "BFGS",
              control = list(maxit = 1e3)),
        silent = TRUE)
    if(class(opt0) != "try-error"){
        if(opt0$convergence == 0){
            sim_df$tmb_aic_static[j] <- 2 * opt0$value + 2 * 3
            sim_df$tmb_bic_static[j] <- 2 * opt0$value + log(nrow(dat)) * 3
            ## states
            srep <- try(summary(sdreport(obj0)), silent = TRUE)
            if(class(srep)[1] != "try-error"){
                srep <- srep[c("lnq", "lnp"), ][rep(c(1, 2), each = dat$n[1]),]
                srep_df0 <- as.data.frame(srep)
                srep_df0$variable <- rownames(srep)
            }else{
                srep_df0 <- expand.grid("Estimate" = NA, "Std. Error" = NA, variable = rep(c("lnq", "lnp"), each = nrow(dat)))
            }
        }
    }
    ##------------
    ## TMB qvary
    ##------------
    obj1 <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec0),
        parameters = list(
            lnq = rep(starts["lnq0"], n),
            lnsdq = log(starts["sd0"]),
            lnp = starts["lnp0"],
            lnsde = log(starts["sd0"])
        ),
        random = c("lnq"),
        DLL = "tvbh_qvary",
        silent = TRUE)
    ##
    opt1 <- try(
        optim(par = obj1$par,
              fn = obj1$fn,
              gr = obj1$gr,
              method = "BFGS",
              control = list(maxit = 1e3)),
        silent = TRUE)
    if(class(opt1) != "try-error"){
        if(opt1$convergence == 0){
            sim_df$tmb_aic_qvary[j] <- 2 * opt1$value + 2 * 4
            sim_df$tmb_bic_qvary[j] <- 2 * opt1$value + log(nrow(dat)) * 4
            ## states
            srep <- try(summary(sdreport(obj1)), silent = TRUE)
            if(class(srep)[1] != "try-error"){
                srepq <- srep[rownames(srep) == "lnq", ]
                srepp <- t(srep[rownames(srep) == "lnp", ])[rep(1, dat$n[1]),]
                ##
                srepq_df <- as.data.frame(srepq)
                srepq_df$variable <- rownames(srepq)
                ##
                srepp_df <- as.data.frame(srepp)
                srepp_df$variable <- "lnp"
                srep_df1 <- rbind(srepq_df, srepp_df)
            }else{
                srep_df1 <- expand.grid("Estimate" = NA, "Std. Error" = NA, variable = rep(c("lnq", "lnp"), each = nrow(dat)))
            }
        }
    }
    ##-----------
    ## TMB pvary 
    ##-----------
    obj2 <- MakeADFun(
        data = list(ssb = s,
                    y = log(r),
                    tau = tau,
                    fec = fec0),
        parameters = list(
            lnq = starts["lnq0"],
            lnp = rep(starts["lnp0"], n),
            lnsdp = log(starts["sd0"]),
            lnsde = log(starts["sd0"])
        ),
        random = c("lnp"),
        DLL = "tvbh_pvary",
        silent = TRUE)
        ##
    opt2 <- try(  
        optim(par = obj2$par,
              fn = obj2$fn,
              gr = obj2$gr,
              method = "BFGS",
              control = list(maxit = 1e3)),
        silent = TRUE)
    if(class(opt2) != "try-error"){
        if(opt2$convergence == 0){
            sim_df$tmb_aic_pvary[j] <- 2 * opt2$value + 2 * 4
            sim_df$tmb_bic_pvary[j] <- 2 * opt2$value + log(nrow(dat)) * 4
            ## states
            srep <- try(summary(sdreport(obj2)), silent = TRUE)
            if(class(srep)[1] != "try-error"){
                srepp <- srep[rownames(srep) == "lnp", ]
                srepq <- t(srep[rownames(srep) == "lnq", ])[rep(1, dat$n[1]),]
                ##
                srepp_df <- as.data.frame(srepp)
                srepp_df$variable <- rownames(srepp)
                ##
                srepq_df <- as.data.frame(srepq)
                srepq_df$variable <- "lnq"
                srep_df2 <- rbind(srepq_df, srepp_df)
            }else{
                srep_df2 <- expand.grid("Estimate" = NA, "Std. Error" = NA, variable = rep(c("lnq", "lnp"), each = nrow(dat)))
            }
        }
    }
    ## get best bssm model
    ic_idx <- grep(paste0("^", "tmb_aic"), names(sim_df))
    best <- apply(sim_df[j, ic_idx], 1, function(z){
        if(all(is.na(z))){
            3
        }else{
            which.min(z)-1 ## as smooths start at 0
        }
    })
    ## and the best fitting states
    if(best == 3){
        ## all are unconverged
    }else{
        srep_df <- get(paste0("srep_df", best))
        dat$qhat_tmb <- exp(subset(srep_df, variable == "lnq")$Estimate)
        dat$phat_tmb <- exp(subset(srep_df, variable == "lnp")$Estimate)
    }
    ##---------------------
    ## fill in the results
    ##---------------------
    for(nn in new_names){
        timeseries[idx, nn] <- dat[, nn]
    }
    ##clean up temporary objects in loop
    rm(list = vars2remove)
    gc()
    ## progress
    cat(paste0(j, "/", m, " = ", round(j/m * 100, 2),"% \n"), file = "progress.txt", append = TRUE)
}

time1 <- Sys.time()

time1 - time0

save(list = c("sim_df", "timeseries"), file = "inference.RData")
