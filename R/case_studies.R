##----------------
## Real data fits
## CM: 7/12/23
## Notes:
## 08-Aug-2024 - used contrast to select examples
## 20-Dec-2024 - used model selection
## 7-Jan-25 included grid of starting values
##----------------

library(ggplot2); theme_set(theme_bw())
library(viridis)

## read in the data
timeseries <- read.csv("case_studies_data.csv")

regions <- unique(timeseries$Region)
regions <- regions[order(regions)]

## get contrast for each stock
stockid <- unique(timeseries$ID)

contrast_df <- data.frame(stockid, contrast = NA)

for(i in 1:nrow(contrast_df)){
    sub <- subset(timeseries, ID == contrast_df$stockid[i])
    Smax <- max(sub$SSB)
    Smin <- min(sub$SSB)
    contrast_df$contrast[i] <- (Smax - Smin)/Smax
}

plot(contrast_df$contrast, ylim = c(0, 1))

## plot out the SR series and contrast for each stock
pdf("../tex/figures/sr_real_data.pdf", height = 6, width = 7)
for(i in 1:nrow(contrast_df)){
    id <- contrast_df$stockid[i]
    sub <- subset(timeseries, ID == id)
    p <- ggplot(sub, aes(x = SSB, y = Recruits, colour = Year)) +
        geom_point() +
        scale_colour_viridis() +
        expand_limits(x = 0, y = 0) +
        theme(legend.position = "bottom") +
        ggtitle(paste(id, round(contrast_df$contrast[i], 2), collapse = " "))
    print(p)    
}
dev.off()

## used this to select three stocks per contrast
quantile(contrast_df$contrast, p = c(0.1, 0.5, 0.9))

## chosen stocks to fit
stocks2fit_df <- read.csv("../data/choose_stocks.csv")

stocks2fit <- stocks2fit_df$stockid

timeseries <- subset(timeseries, ID %in% stocks2fit)

##------------
## ESTIMATION 
##------------
source("utils.R")
## for confidence intervals
beta <- matrix(c(1, 0, 1, -2, 1, 2), nr = 2)

library(bssm)
Rcpp::sourceCpp("bh_bssm_new.cpp")
pntrs <- create_xptrs()

library(dplyr)
library(diagis)

## make estimate slots for ekf and tmb (these are run first)
names_df <- expand.grid(c("q", "p"),
                        c("hat_ekf", "lwr_ekf", "upr_ekf",
                          "hat_tmb", "lwr_tmb", "upr_tmb"))

new_names <- apply(names_df, 1, paste, collapse = "")

new_names <- new_names[order(new_names)]

timeseries[, new_names] <- NA

## make information criteria table
real_df <- data.frame(ID = stocks2fit)

m <- nrow(real_df)

## slots for information criteria
names <- expand.grid(method = c("ekf", "tmb"),
                     ic = "aic",
                     model = c("static", "qvary", "pvary"))

names <- rbind(names,
               expand.grid(method = "bssm", ic = "dic", model = c("static", "qvary", "pvary")))

real_df[, apply(names, 1, paste, collapse = "_")] <- NA

## variables to remove at end of each iteration
vars2remove <- c("dat", "fec0", "starts", "mod", "build", "f0", "f1", "f2", "ic_idx", "best", "smooth0", "smooth1", "smooth2", "smooth", "initial_theta", "s", "r", "year", "n", "tau", "obj0", "opt0", "srep", "srep_df0", "obj1", "opt1", "srep_df1", "obj2", "opt2", "srep_df2")

time0 <- Sys.time()

## grid of starting values in time-invariant fits
sde_vec <- c(0.1, 0.2, 0.4, 0.6)

## grid of starting values for standard deviations in time-varying fits
start_grid <- expand.grid(sdproc = c(0.01, 0.05, 0.1, 0.2),
                          sde = c(0.1, 0.2, 0.4, 0.6))


## quick look at the staring values
par(mfrow = c(2, 3))
for(j in 1:m){
    print(j)
    ##
    idx <- which(timeseries$ID == real_df$ID[j])
    dat <- timeseries[idx,]
    fec0 <- 100
    tau <- 1
    ##~~~~~~~~~~
    ## eKF fits 
    ##~~~~~~~~~~
    starts <- get_start(dat, fec = fec0, tau = tau)
    with(dat, plot(SSB, Recruits, xlim = c(0, max(SSB)), ylim = c(0, max(Recruits)), bty = "l"))
    curve(1/(starts["a0"] + starts["b0"]/x), add = TRUE)
    legend("topleft", legend = real_df$ID[j], bty = "n")
}

## run eKF and TMB
for(j in 1:m){
    print(j)
    ##
    idx <- which(timeseries$ID == real_df$ID[j])
    dat <- timeseries[idx,]
    ## set fecundity
    fec0 <- 100
    tau <- 1
    starts <- get_start(dat, fec = fec0, tau = tau)
    ## with(dat, plot(SSB, Recruits, xlim = c(0, max(SSB)), ylim = c(0, max(Recruits))))
    ## curve(1/(starts["a0"] + starts["b0"]/x), add = TRUE)    
    ##~~~~~~~~~~~~~
    ## eKF fits
    ##~~~~~~~~~~~~~
    ##--------------
    ## STATIC MODEL
    ##--------------
    ## setup the model
    mod <- list(H = matrix(0),
                T = diag(2),
                Q = diag(2) * 0,
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
    ##
    f0_list <- list()
    for(k in 1:length(sde_vec)){
        f0_list[[k]] <- try(optim(
            par = log(sde_vec[k]),
            fn = nll_ekf,
            fec = fec0,
            mod = mod,
            data = dat,
            control = list(maxit = 1e3),
            method = "BFGS"),
            silent = TRUE)
    }
    f0_ll <- unlist(lapply(f0_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            x$value
        }
    }))
    f0 <- f0_list[[which.min(f0_ll)]]
    if(class(f0) != "try-error"){
        if(f0$convergence == 0){
            real_df$ekf_aic_static[j] <- 2 * f0$value + 2 * 3
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
    f1_list <- list()
    for(k in 1:nrow(start_grid)){
        f1_list[[k]] <- try(optim(
            par = log(start_grid[k,]),
            fn = nll_ekf,
            fec = fec0,
            mod = mod,
            data = dat,
            control = list(maxit = 1e3),
            method = "BFGS"),
            silent = TRUE)
    }
    f1_ll <- unlist(lapply(f1_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            x$value
        }
    }))
    f1 <- f1_list[[which.min(f1_ll)]]
    ##
    if(class(f1) != "try-error" ){
        ## sometimes matrix is singular so calculating all pieces first
        ## smoothed states
        mod_hat <- build(f1$par, mod)
        filt <- bh_ekf(mod = mod_hat, fec0, data = dat)
        smooth1 <- try(my_ks(mod_hat, at = filt$at, Pt = filt$Pt, Pttm1 = filt$Pttm1), silent = TRUE)
        if(f1$convergence == 0 & class(smooth1) != "try-error"){
            real_df$ekf_aic_qvary[j] <- 2 * f1$value + 2 * 4
        }
    }
    ##--------
    ## p-vary 
    ##--------
    build <- function(theta, mod){
        mod$Q[1,1] <- 0 ## so thetas are log sds
        mod$Q[2,2] <- exp(theta[1])^2
        mod$H[1,1] <- exp(theta[2])^2 
        return(mod)
    }
    ##
    f2_list <- list()
    for(k in 1:nrow(start_grid)){
        f2_list[[k]] <- try(optim(
            par = log(start_grid[k,]),
            fn = nll_ekf,
            fec = fec0,
            mod = mod,
            data = dat,
            control = list(maxit = 1e3),
            method = "BFGS"),
            silent = TRUE)
    }
    f2_ll <- unlist(lapply(f2_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            x$value
        }
    }))
    f2 <- f2_list[[which.min(f2_ll)]]
    ##
    if(class(f2) != "try-error" ){
        ## smoothed states
        mod_hat <- build(f2$par, mod)
        filt <- bh_ekf(mod = mod_hat, fec0, data = dat)
        smooth2 <- try(my_ks(mod_hat, at = filt$at, Pt = filt$Pt, Pttm1 = filt$Pttm1), silent = TRUE)
        ##
        if(f2$convergence == 0 & class(smooth2) != "try-error"){
            real_df$ekf_aic_pvary[j] <- 2 * f2$value + 2 * 4
        }
    }
    ## get best ekf model
    ic_idx <- grep(paste0("^", "ekf_aic"), names(real_df))
    best <- apply(real_df[j, ic_idx], 1, function(z){
        if(all(is.na(z))){
            3
        }else{
            which.min(z)-1 ## as smooths start at 0
        }
    })
    ## and the best fitting states
    if(best == 3){
        ## all are unconverged, leave as NA's
    }else{
        smooth <- get(paste0("smooth", best))
        if(class(smooth) == "try-error"){
            dat$qhat_ekf <- NA
            dat$phat_ekf <- NA
        }else{
            lnq_atT <- smooth$atT[1,-1]
            lnq_PtT <- smooth$PtT[1,1,-1]
            lnq_se_atT <- sqrt(lnq_PtT)
            lnq_smooth_atT <- cbind(lnq_atT, lnq_se_atT) %*% beta
            ## qt
            lnp_atT <- smooth$atT[2,-1]
            lnp_PtT <- smooth$PtT[2,2,-1]
            lnp_se_atT <- sqrt(lnp_PtT)
            lnp_smooth_atT <- cbind(lnp_atT, lnp_se_atT) %*% beta
            ##
            qmat <- exp(lnq_smooth_atT)
            dat$qhat_ekf <- qmat[,1]
            dat$qlwr_ekf <- qmat[,2]
            dat$qupr_ekf <- qmat[,3]
            ##
            pmat <- exp(lnp_smooth_atT)
            dat$phat_ekf <- pmat[,1]
            dat$plwr_ekf <- pmat[,2]
            dat$pupr_ekf <- pmat[,3]
        }
    }
    ##~~~~~~~
    ## TMB
    ##~~~~~~~
    s <- dat$SSB
    r <- dat$Recruits
    year <- dat$Year
    n <- length(year)
    if(real_df$ID[j] == "cod.27.7a"){
        starts["lnp0"] <- -10 ## otherwise non-finite starting values 
    }
    ##------------
    ## TMB static
    ##------------
    sdreport0_list <- list()
    opt0_ll <- rep(NA, length(sde_vec))
    for(k in 1:length(sde_vec)){
        obj0 <- MakeADFun(
            data = list(ssb = s,
                        y = log(r),
                        tau = tau,
                        fec = fec0),
            parameters = list(
                lnq = starts["lnq0"],
                lnp = starts["lnp0"],
                lnsde = log(sde_vec[k])
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
            sdrep0 <- summary(sdreport(obj0))
            sdreport0_list[[k]] <- sdrep0
            if(any(is.na(sdrep0[, "Std. Error"]) | is.nan(sdrep0[, "Std. Error"]) | opt0$convergence != 0)){
                opt0_ll[k] <- Inf
            }else{
                opt0_ll[k] <- opt0$value
            }
        }else{
            sdreport0_list[[k]] <- NA
            opt0_ll[k] <- Inf
        }
        rm(list = c("obj0", "opt0"))
    }
    ## find one with highest likelihood
    if(any(is.finite(opt0_ll))){
        idx0 <- which.min(opt0_ll)
        real_df$tmb_aic_static[j] <- 2 * opt0_ll[idx0] + 2 * 3
        srep <- sdreport0_list[[idx0]]
        srep <- srep[c("lnq", "lnp"), ][rep(c(1, 2), each = nrow(dat)),]
        srep_df0 <- as.data.frame(srep)
        srep_df0$variable <- rownames(srep)
    }
    ##------------
    ## TMB qvary
    ##------------
    sdreport1_list <- list()
    opt1_ll <- rep(NA, length(sde_vec))
    for(k in 1:nrow(start_grid)){
        obj1 <- MakeADFun(
            data = list(ssb = s,
                        y = log(r),
                        tau = tau,
                        fec = fec0),
            parameters = list(
                lnq = rep(starts["lnq0"], n),
                lnsdq = log(start_grid$sdproc[k]),
                lnp = starts["lnp0"],
                lnsde = log(start_grid$sde[k])
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
            sdrep1 <- summary(sdreport(obj1))
            sdreport1_list[[k]] <- sdrep1
            if(any(is.na(sdrep1[, "Std. Error"]) | is.nan(sdrep1[, "Std. Error"]) | opt1$convergence != 0)){
                opt1_ll[k] <- Inf
            }else{
                opt1_ll[k] <- opt1$value
            }
        }else{
            sdreport1_list[[k]] <- NA
            opt1_ll[k] <- Inf
        }
        rm(list = c("obj1", "opt1"))
    }
    ## find one with highest likelihood
    if(any(is.finite(opt1_ll))){
        idx1 <- which.min(opt1_ll)
        real_df$tmb_aic_qvary[j] <- 2 * opt1_ll[idx1] + 2 * 4
        srep <- sdreport1_list[[idx1]]
        srepq <- srep[rownames(srep) == "lnq", ]
        srepp <- t(srep[rownames(srep) == "lnp", ])[rep(1, nrow(dat)),]
        ##
        srepq_df <- as.data.frame(srepq)
        srepq_df$variable <- rownames(srepq)
        ##
        srepp_df <- as.data.frame(srepp)
        srepp_df$variable <- "lnp"
        srep_df1 <- rbind(srepq_df, srepp_df)
    }
    ##-----------
    ## TMB pvary 
    ##-----------
    sdreport2_list <- list()
    opt2_ll <- rep(NA, length(sde_vec))
    for(k in 1:nrow(start_grid)){
        obj2 <- MakeADFun(
            data = list(ssb = s,
                        y = log(r),
                        tau = tau,
                        fec = fec0),
            parameters = list(
                lnq = starts["lnq0"],
                lnp = rep(starts["lnp0"], n),
                lnsdp = log(start_grid$sdproc[k]),
                lnsde = log(start_grid$sde[k])
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
            sdrep2 <- try(summary(sdreport(obj2)), silent = TRUE)
            if(!"try-error" %in% class(sdrep2)){
                sdreport2_list[[k]] <- sdrep2
                if(any(is.na(sdrep2[, "Std. Error"]) | is.nan(sdrep2[, "Std. Error"]) | opt2$convergence != 0)){
                    opt2_ll[k] <- Inf
                }else{
                    opt2_ll[k] <- opt2$value
                }
            }else{
                sdreport2_list[[k]] <- NA
                opt2_ll[k] <- Inf
            }
        }else{
               sdreport2_list[[k]] <- NA
                opt2_ll[k] <- Inf 
        }
        rm(list = c("obj2", "opt2"))
    }
    ## find one with highest likelihood
    if(any(is.finite(opt2_ll))){
        idx2 <- which.min(opt2_ll)
        real_df$tmb_aic_pvary[j] <- 2 * opt2_ll[idx2] + 2 * 4
        srep <- sdreport2_list[[idx2]]
        srepp <- srep[rownames(srep) == "lnp", ]
        srepq <- t(srep[rownames(srep) == "lnq", ])[rep(1, nrow(dat)),]
        ##
        srepp_df <- as.data.frame(srepp)
        srepp_df$variable <- rownames(srepp)
        ##
        srepq_df <- as.data.frame(srepq)
        srepq_df$variable <- "lnq"
        srep_df2 <- rbind(srepq_df, srepp_df)
    }   
    ## get best tmb model
    ic_idx <- grep(paste0("^", "tmb_aic"), names(real_df))
    best <- apply(real_df[j, ic_idx], 1, function(z){
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
        dat$qlwr_tmb <- exp(subset(srep_df, variable == "lnq")$Estimate - 1.96 * subset(srep_df, variable == "lnq")[, "Std. Error"])
        dat$qupr_tmb <- exp(subset(srep_df, variable == "lnq")$Estimate + 1.96 * subset(srep_df, variable == "lnq")[, "Std. Error"])
        ##
        dat$phat_tmb <- exp(subset(srep_df, variable == "lnp")$Estimate)
        dat$plwr_tmb <- exp(subset(srep_df, variable == "lnp")$Estimate - 1.96 * subset(srep_df, variable == "lnp")[, "Std. Error"])
        dat$pupr_tmb <- exp(subset(srep_df, variable == "lnp")$Estimate + 1.96 * subset(srep_df, variable == "lnp")[, "Std. Error"])        
    }
    ##---------
    ## fill in
    ##---------
    for(nn in new_names){
        timeseries[idx, nn] <- dat[, nn]
    }
    ##clean up objects in loop
    rm(list = vars2remove)
    gc()
    ## progress
    ##cat(paste0(j, "/", m, " = ", round(j/m * 100, 2),"% \n"), file = "progress.txt", append = TRUE)
}

##----------
## run BSSM
##----------
names_df <- expand.grid(c("q", "p"), c("hat_bssm", "lwr_bssm", "upr_bssm"))

new_names <- apply(names_df, 1, paste, collapse = "")

new_names <- new_names[order(new_names)]

timeseries[, new_names] <- NA

vars2remove <- c("dat", "fec0", "starts", "initial_theta", "known_params", "known_tv_params", "model", "mcmc_static", "dic", "mcmc_qvary")

niter_bssm <- 2e4

start_bssm <- expand.grid(sdproc = c(0.01, 0.1),
                          sde = c(0.2, 0.4))

sde_vec <- c(0.2, 0.4)

## wrap in a pdf to look at convergence
pdf("bssm_convergence.pdf", height = 9, width = 11)
for(j in 1:m){
    print(j) ##
    idx <- which(timeseries$ID == real_df$ID[j])
    dat <- timeseries[idx,] ##dat$SSB <- dat$s
    ## set fecundity
    fec0 <- 100
    tau <- 1
    starts <- get_start(dat, fec0, tau = tau)
    ## with(dat, plot(SSB, Recruits, xlim = c(0, max(SSB)), ylim = c(0, max(Recruits))))
    ## curve(1/(starts["a0"] + starts["b0"]/x), add = TRUE)
    ##-------------
    ## BSSM STATIC
    ##-------------
    ## run two chains
    known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                      ##P11 = 1e2, P22 = 1e2,
                      P11 = 1e0, P22 = 1e0,
                      tau = tau,
                      fec = fec0)
    ##
    known_tv_params <- matrix(dat$SSB)
    ##
    bssm_static <- list()
    ##
    for(k in 1:length(sde_vec)){
        initial_theta <- c(log_H = log(sde_vec[k]))
        ##
        model <- ssm_nlg(y = log(dat$Recruits),
                         a1=pntrs$a1,
                         P1 = pntrs$P1, 
                         Z = pntrs$Z_fn,
                         H = pntrs$H_fn,
                         T = pntrs$T_fn,
                         ##R = pntrs$R_fn,
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
        ## running for much longer
        mcmc_static <- try(run_mcmc(model, particles = 50, iter = niter_bssm, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(1) * 0.1), silent = TRUE)
        ##
        bssm_static[[k]] <- mcmc_static
    }
    sde1 <- bssm_static[[1]]$theta[rep(1:nrow(bssm_static[[1]]$theta), times = bssm_static[[1]]$count),]
    sde2 <- bssm_static[[2]]$theta[rep(1:nrow(bssm_static[[2]]$theta), times = bssm_static[[2]]$count),]
    lnsde_chain <- cbind(sde1, sde2) ##
    par(mfrow = c(1, 2))
    matplot(lnsde_chain, type = "l", lty = 1, col = 1:2, bty = "l", main = paste(real_df$ID[j], "BSSM static"))
    plot(density(sde1), bty = "l")
    lines(density(sde2), col = "red") ##
    ##
    dic_static_vec <- unlist(lapply(bssm_static, DIC, static = TRUE))
    mcmc_static <- bssm_static[[which.min(dic_static_vec)]]
    ##
    if(class(mcmc_static) != "try-error"){
        dic <- try(DIC(mcmc_static, static = TRUE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            real_df$bssm_dic_static[j] <- dic
            states0 <- summary(mcmc_static, variable = "states")
        }
    }
    ##-------------
    ## BSSM q vary
    ##-------------
    bssm_qvary <- list()
    for(k in 1:nrow(start_bssm)){
        initial_theta <- c(log(start_bssm$sde[k]), log(start_bssm$sdproc[k]))        
        ##
        known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                          P11 = 1e0, P22 = 1e0,
                          tau = tau,
                          fec = fec0)
        ##
        model <- ssm_nlg(y = log(dat$Recruits),
                         a1 = pntrs$a1,
                         P1 = pntrs$P1, 
                         Z = pntrs$Z_fn,
                         H = pntrs$H_fn,
                         T = pntrs$T_fn,
                         ##R = pntrs$R_fn,
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
        mcmc_qvary <- try(run_mcmc(model, particles = 50, iter = niter_bssm, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(2) * 0.1), silent = TRUE)
        bssm_qvary[[k]] <- mcmc_qvary
    }
    sde_list <- lapply(bssm_qvary, function(x){
                             x$theta[rep(1:nrow(x$theta), times = x$count), 1]
    })
    sde_mat <- do.call(cbind, sde_list)
    sdq_list <- lapply(bssm_qvary, function(x){
        x$theta[rep(1:nrow(x$theta), times = x$count), 2]
    })
    sdq_mat <- do.call(cbind, sdq_list)
    ##
    par(mfrow = c(2, 2))
    ## sde
    matplot(sde_mat, type = "l", lty = 1, col = 1:4, bty = "l", main = paste(real_df$ID[j], "BSSM qvary ln(sde)"))
    d0 <- density(sde_mat[,1])
    plot(d0, bty = "l", ylim = c(0, 1.2 * max(d0$y)), main = "ln(sde)")
    for(k in 2:4){
        lines(density(sde_mat[,k]), col = k)
    }
    ## sdq
    matplot(sdq_mat, type = "l", lty = 1, col = 1:4, bty = "l", main = paste(real_df$ID[j], "BSSM qvary ln(sdq)"))
    d0 <- density(sdq_mat[,1])
    plot(d0, bty = "l", ylim = c(0, 1.2 * max(d0$y)), main = "ln(sdq)")
    for(k in 2:4){
        lines(density(sdq_mat[,k]), col = k)
    }
    dic_qvary_vec <- unlist(lapply(bssm_qvary, function(x){
        dic <- try(DIC(x, static = FALSE), silent = TRUE)
        if(class(dic) == "try-error"){
            Inf
        }else{
            dic
        }
    }))
    mcmc_qvary <- bssm_qvary[[which.min(dic_qvary_vec)]]
    ##
    if(class(mcmc_qvary) != "try-error"){
        dic <- try(DIC(mcmc_qvary, static = FALSE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            real_df$bssm_dic_qvary[j] <- dic
            states1 <- summary(mcmc_qvary, variable = "states")            
        }
    }
    ##-------------
    ## BSSM p vary
    ##-------------
    bssm_pvary <- list()
    for(k in 1:nrow(start_bssm)){
        initial_theta <- c(log(start_bssm$sde[k]), log(start_bssm$sdproc[k]))
        known_params <- c(a11 = starts["lnq0"], a22 = starts["lnp0"],
                          P11 = 1e0, P22 = 1e0,
                          tau = tau,
                          fec = fec0)
        ##
        model <- ssm_nlg(y = log(dat$Recruits),
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
        mcmc_pvary <- try(run_mcmc(model, particles = 50, iter = niter_bssm, mcmc_type = "is2", sampling_method = "psi", verbose = FALSE, S = diag(2) * 0.1), silent = TRUE)
        bssm_pvary[[k]] <- mcmc_pvary
    }
    ## plot the chains
    sde_list <- lapply(bssm_pvary, function(x){
        x$theta[rep(1:nrow(x$theta), times = x$count), 1]
    })
    sde_mat <- do.call(cbind, sde_list)
    ##
    sdp_list <- lapply(bssm_pvary, function(x){
        x$theta[rep(1:nrow(x$theta), times = x$count), 2]
    })
    sdp_mat <- do.call(cbind, sdp_list)
    ##
    par(mfrow = c(2, 2))
    ## sde
    matplot(sde_mat, type = "l", lty = 1, col = 1:4, bty = "l", main = paste(real_df$ID[j], "BSSM pvary ln(sde)"))
    d0 <- density(sde_mat[,1])
    plot(d0, bty = "l", ylim = c(0, 1.2 * max(d0$y)), xlim = range(sde_mat), main = "ln(sde)")
    for(k in 2:4){
        lines(density(sde_mat[,k]), col = k)
    }
    ## sdp
    matplot(sdp_mat, type = "l", lty = 1, col = 1:4, bty = "l", main = paste(real_df$ID[j], "BSSM pvary ln(sdp)"))
    d0 <- density(sdp_mat[,1])
    plot(d0, bty = "l", ylim = c(0, 1.2 * max(d0$y)), xlim = range(sdp_mat), main = "ln(sdp)")
    for(k in 2:4){
        lines(density(sdp_mat[,k]), col = k)
    }
    dic_pvary_vec <- unlist(lapply(bssm_pvary, function(x){
        dic <- try(DIC(x, static = FALSE), silent = TRUE)
        if(class(dic) == "try-error"){
            Inf
        }else{
            dic
        }
    }))
    mcmc_pvary <- bssm_pvary[[which.min(dic_pvary_vec)]]
    ## 
    if(class(mcmc_pvary) != "try-error"){
        dic <- try(DIC(mcmc_pvary, static = FALSE), silent = TRUE)
        if(class(dic) != "try-error" & !is.nan(dic)){
            real_df$bssm_dic_pvary[j] <- dic
            states2 <- summary(mcmc_pvary, variable = "states")
        }
    }
    ## get best bssm model
    ic_idx <- grep(paste0("^", "bssm_dic"), names(real_df))
    best <- apply(real_df[j, ic_idx], 1, function(z){
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
        states <- subset(states, time != (nrow(dat) + 1))
        ## might be better of exponentiating this and then taking the mean
        ## do this later
        dat$qhat_bssm <- exp(subset(states, variable == "logq")$Mean)
        dat$qlwr_bssm <- exp(subset(states, variable == "logq")[, "2.5%"])
        dat$qupr_bssm <- exp(subset(states, variable == "logq")[, "97.5%"])
        ##
        dat$phat_bssm <- exp(subset(states, variable == "logp")$Mean)
        dat$plwr_bssm <- exp(subset(states, variable == "logp")[, "2.5%"])
        dat$pupr_bssm <- exp(subset(states, variable == "logp")[, "97.5%"])
    }
    ##---------
    ## fill in
    ##---------
    for(nn in new_names){
        timeseries[idx, nn] <- dat[, nn]
    }
    ##clean up objects in loop
    rm(list = vars2remove)
    gc()
}
dev.off()

## save(timeseries, file = "timeseries_fits_14_Jan_2025.RData")

## calculate difference from minimum information criterion per estimation method
real_df1 <- real_df

for(i in 1:nrow(real_df1)){
    for(j in c("ekf", "bssm", "tmb")){
        idx <- grep(j, names(real_df1))
        real_df1[i, idx] <- real_df1[i, idx] - min(real_df1[i, idx], na.rm = TRUE)
    }
}


real_df2 <- cbind(real_df1[, 1], round(real_df1[, c("ekf_aic_static", "ekf_aic_qvary", "ekf_aic_pvary",
                                                 "bssm_dic_static", "bssm_dic_qvary", "bssm_dic_pvary",
                                                 "tmb_aic_static", "tmb_aic_qvary", "tmb_aic_pvary")], 1))

## output information criteria table
library(xtable)

real_tab <- xtable(real_df2, digits = 1)

print.xtable(real_tab, file = "../tex/tables/real_information_criteria.tex", include.rownames = FALSE)


## plot the estimated series
mycols <- c("darkgrey", "orange2", "darkblue")

## plots of parameters only
vars2keep <- c("ID", "Year",
               "qhat_ekf", "qlwr_ekf", "qupr_ekf",
               "phat_ekf", "plwr_ekf", "pupr_ekf",
               "qhat_tmb", "qlwr_tmb", "qupr_tmb",
               "phat_tmb", "plwr_tmb", "pupr_tmb",
               "qhat_bssm", "qlwr_bssm", "qupr_bssm",
               "phat_bssm", "plwr_bssm", "pupr_bssm")


## parameter estimates
pdf("../tex/figures/Fig6_real_par_ests_Jan_2025.pdf", height = 7, width = 10)
range_mult <- c(0.5, 1.5)
par(mfrow = c(3, 4), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 2))
for(j in 1:length(stocks2fit)){
    tmp <- subset(timeseries, ID == stocks2fit[j])
    qmin <- range_mult[1] * min(tmp[, grep("^qhat", names(tmp))], na.rm = TRUE)
    qmax <- range_mult[2] * max(tmp[, grep("^qhat", names(tmp))], na.rm = TRUE)    
    qrange <- c(qmin, qmax)
    with(tmp, matplot(Year, cbind(qhat_ekf, qlwr_ekf, qupr_ekf), col = mycols[1], type = "l", lty = c(1, 2, 2), ylim = qrange))
    with(tmp, matlines(Year, cbind(qhat_bssm, qlwr_bssm, qupr_bssm), col = mycols[2], type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(qhat_tmb, qlwr_tmb, qupr_tmb), col = mycols[3], type = "l", lty = c(1, 2, 2)))
    legend("topleft", legend = paste(stocks2fit[j], "q"), bty = "n")
    ##
    pmin <- range_mult[1] * min(tmp[, grep("^phat", names(tmp))], na.rm = TRUE)    
    pmax <- range_mult[2] * max(tmp[, grep("^phat", names(tmp))], na.rm = TRUE)
    prange <- c(pmin, pmax)
    with(tmp, matplot(Year, cbind(phat_ekf, plwr_ekf, pupr_ekf), col = mycols[1], type = "l", lty = c(1, 2, 2), ylim = prange))
    with(tmp, matlines(Year, cbind(phat_bssm, plwr_bssm, pupr_bssm), col = mycols[2], type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(phat_tmb, plwr_tmb, pupr_tmb), col = mycols[3], type = "l", lty = c(1, 2, 2)))
    legend("topleft", legend = paste(stocks2fit[j], "p"), bty = "n")
}
mtext(side = 1, text = "Year", line = 0.5, outer = TRUE)
mtext(side = 2, text = "Parameter estimate", line = 0.5, outer = TRUE)
mtext(side = 4, text = "Low contrast", at = 0.85, line = 0.25, outer = TRUE)
mtext(side = 4, text = "Medium contrast", at = 0.5, line = 0.25, outer = TRUE)
mtext(side = 4, text = "High contrast", at = 0.17, line = 0.25, outer = TRUE)
dev.off()


## for the appendix - all fits
pdf("../tex/figures/real_fits_all_3D.pdf", height = 10, width = 8)
for(k in 0:1){
    par(mfrow = c(3, 3), oma = c(4, 4, 1, 3))
    for(j in (k*3 + 1):(k*3 + 3)){
        print(j)
        tmp <- subset(timeseries, ID == stocks2fit[j])
        fec0 <- 100
        cols <- viridis(nrow(tmp))
        maxS <- max(tmp$SSB)
        minS <- min(tmp$SSB)
        maxR <- max(tmp$Recruits)
        maxY <- max(tmp$Year)
        minY <- min(tmp$Year)
        ##
        for(method in c("ekf", "bssm", "tmb")){
            if(stocks2fit[j] == "cod.27.7a"){scale.y = 0.7}else{scale.y = 0.5}
            sr3d <- with(tmp, scatterplot3d(x = Year, y = SSB, z = Recruits, xlim = c(minY, maxY), ylim = c(0, maxS), zlim = c(0, maxR), pch = NA,  color = "darkgrey", box = FALSE, type = "h", angle = 60, grid = TRUE, cex.symbols = 2, scale.y = scale.y, mar = c(0, 0, 0, 0), xlab = "", ylab = "", zlab = ""))
            with(tmp, sr3d$points3d(x = Year, y = SSB, z = Recruits, col = cols, pch = 19, cex = 1.5))
            x <- seq(0, max(tmp$SSB), length = 1e3)
            sapply(1:nrow(tmp), function(zz){
                qt <- tmp[zz, paste0("qhat_", method)]
                pt <- tmp[zz, paste0("phat_", method)]
                at <- pt/qt * (exp(qt) - 1)
                bt <- exp(qt)/fec0
                z <- 1/(at + bt/x)
                y <- rep(tmp[zz, "Year"], length(x))
                sr3d$points3d(x = y, y = x, z = z, col = cols[zz], type = "l", lwd = 1)
            })
            if(method == "bssm") mtext(side = 3, text = stocks2fit[j], line = -3, cex = 1)
            if(method == "tmb"){
                coords <- sr3d$xyz.convert(maxY + 15, (maxS + minS)/3, 0)
                text(coords$x, coords$y, label = "SSB", srt = 0, cex = 1.5, xpd = NA)
            }        
            coords <- sr3d$xyz.convert(minY-5, (maxS + minS)/4, maxR)
            if(j %in% c(1, 4)) legend(coords$x, coords$y, legend = toupper(method), bty = "n", cex = 1.5)
        }
        mtext(side = 1, text = "Year", line = 1.5, outer = TRUE)
        mtext(side = 2, text = "Recruits", line = 1.5, outer = TRUE)    
    }
}
dev.off()


## BSSM fits for main ms
pdf("../tex/figures/real_fits_bssm_3D.pdf", height = 10, width = 8)
par(mfrow = c(3, 2), oma = c(4, 4, 1, 3))
for(j in 1:m){
    print(j)
    tmp <- subset(timeseries, ID == stocks2fit[j])
    fec0 <- 100
    cols <- viridis(nrow(tmp))
    maxS <- max(tmp$SSB)
    minS <- min(tmp$SSB)
    maxR <- max(tmp$Recruits)
    maxY <- max(tmp$Year)
    minY <- min(tmp$Year)
    ##
    method <- "bssm"
    if(stocks2fit[j] == "cod.27.7a"){scale.y = 0.7}else{scale.y = 1}
    sr3d <- with(tmp, scatterplot3d(x = Year, y = SSB, z = Recruits, xlim = c(minY, maxY), ylim = c(0, maxS), zlim = c(0, maxR), pch = NA,  color = "darkgrey", box = FALSE, type = "h", angle = 60, grid = TRUE, cex.symbols = 2, scale.y = scale.y, mar = c(0, 0, 0, 0), xlab = "", ylab = "", zlab = ""))
    with(tmp, sr3d$points3d(x = Year, y = SSB, z = Recruits, col = cols, pch = 19, cex = 1.5))
    x <- seq(0, max(tmp$SSB), length = 1e3)
    sapply(1:nrow(tmp), function(zz){
        qt <- tmp[zz, paste0("qhat_", method)]
        pt <- tmp[zz, paste0("phat_", method)]
        at <- pt/qt * (exp(qt) - 1)
        bt <- exp(qt)/fec0
        z <- 1/(at + bt/x)
        y <- rep(tmp[zz, "Year"], length(x))
        sr3d$points3d(x = y, y = x, z = z, col = cols[zz], type = "l", lwd = 1)
    })
    mtext(side = 3, text = stocks2fit[j], line = -4)
    if(j == 2) mtext(side = 4, text = "Low contrast", line = 1)
    if(j == 4){
        mtext(side = 4, text = "Medium contrast", line = 1)
        coords <- sr3d$xyz.convert(maxY + 15, (maxS + minS)/2.5, 0)
        text(coords$x, coords$y, label = "SSB", srt = 0, cex = 1.5, xpd = NA)
    }
    if(j == 6) mtext(side = 4, text = "High contrast", line = 1)
    ##box()
}
mtext(side = 1, text = "Year", line = 1.5, outer = TRUE)
mtext(side = 2, text = "Recruits", line = 1.5, outer = TRUE)
dev.off()

##-------------------
## RICKER COMPARISON 
##-------------------
## comparison with the Ricker
names_df <- expand.grid(c("q", "p"), c("hat_ricker", "lwr_ricker", "upr_ricker"))

new_names <- apply(names_df, 1, paste, collapse = "")

new_names <- new_names[order(new_names)]

timeseries[, new_names] <- NA

ids <- unique(timeseries$ID)

m <- length(ids)

library(dlm)

sde_vec <- c(0.1, 0.2, 0.4, 0.6)

ricker_aic_df <- data.frame(ID = ids, kf_static = NA, kf_qvary = NA, kf_pvary = NA)

for(j in 1:m){
    print(j)
    ##
    idx <- which(timeseries$ID == ids[j])
    dat <- timeseries[idx,]
    if(!all(diff(dat$Year) == 1)){
        stop("Years not in order")
    }
    fec0 <- 100
    ##---------
    ## KF fits
    ##---------
    X <- as.matrix(dat$SSB)
    y <- with(dat, log(Recruits/SSB))
    ##-----------
    ## KF static
    ##-----------
    build0 <- function(u) {
        dlm(FF = matrix(c(1,1), ncol = 2),
            V = exp(u[1]),
            GG = diag(rep(1,2)), 
            W = 0 * diag(2),
            JFF = matrix(c(0,1), ncol = 2),
            X = X,
            m0 = c(0, 0),
            C0 = diag(rep(1e+7,2))
            )
    }
    ##
    kf0_list <- list()
    for(k in 1:length(sde_vec)){
        kf0_list[[k]] <- try(
            dlmMLE(y, parm = log(sde_vec[k]), build0, hessian=TRUE),
            silent = TRUE)
    }
    kf0_ll <- unlist(lapply(kf0_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            x$value
        }
    }))
    outMLE <- kf0_list[[which.min(kf0_ll)]]    
    if(class(outMLE) != "try-error"){
        if(outMLE$convergence == 0){
            ricker_aic_df$kf_static[j] <- 2 * outMLE$value + 2 * 3
            ##real_df$ekf_bic_static[j] <- 2 * f0$value + log(nrow(dat)) * 3
            ## smoothed states
            mod_hat0 <- build0(outMLE$par)
        }
    }
    ##-----------
    ## KF qvary
    ##-----------
    build1 <- function(u) {
        dlm(FF = matrix(c(1,1), ncol = 2),
            V = exp(u[1]),
            GG = diag(rep(1,2)), 
            W = diag(c(exp(u[2]), 0)),
            JFF = matrix(c(0,1), ncol = 2),
            ##JGG = matrix(c(0,0,0,0), nrow = 2),
            X = X,
            m0 = c(0, 0),
            C0 = diag(rep(1e+7,2))
            )
    }
    ##
    kf1_list <- list()
    for(k in 1:nrow(start_grid)){
        kf1_list[[k]] <- try(
            dlmMLE(y, parm = c(log(start_grid$sde[k]), log(start_grid$sdproc[k])), build1, hessian=TRUE),
            silent = TRUE)
    }
    kf1_ll <- unlist(lapply(kf1_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            if(x$convergence != 0){
                Inf
            }else{
                x$value
            }
        }
    }))
    outMLE <- kf1_list[[which.min(kf1_ll)]]
    if(class(outMLE) != "try-error"){
        if(outMLE$convergence == 0){
            ricker_aic_df$kf_qvary[j] <- 2 * outMLE$value + 2 * 4
            ##real_df$ekf_bic_static[j] <- 2 * f0$value + log(nrow(dat)) * 3
            ## smoothed states
            mod_hat1 <- build1(outMLE$par)
        }
    }
    ##-----------
    ## KF pvary
    ##-----------
    build2 <- function(u) {
        dlm(FF = matrix(c(1,1), ncol = 2),
            V = exp(u[1]),
            GG = diag(rep(1,2)), 
            W = diag(c(0, exp(u[2]))),
            JFF = matrix(c(0,1), ncol = 2),
            ##JGG = matrix(c(0,0,0,0), nrow = 2),
            X = X,
            m0 = c(0, 0),
            C0 = diag(rep(1e+7,2))
            )
    }
    ##
    kf2_list <- list()
    for(k in 1:nrow(start_grid)){
        kf2_list[[k]] <- try(
            dlmMLE(y, parm = c(log(start_grid$sde[k]), log(start_grid$sdproc[k])), build2, hessian=TRUE),
            silent = TRUE)
    }
    kf2_ll <- unlist(lapply(kf2_list, function(x){
        if(class(x) == "try-error"){
            Inf
        }else{
            x$value
        }
    }))
    outMLE <- kf2_list[[which.min(kf2_ll)]]
    if(class(outMLE) != "try-error"){
        if(outMLE$convergence == 0){
            ricker_aic_df$kf_pvary[j] <- 2 * outMLE$value + 2 * 4
            ## smoothed states
            mod_hat2 <- build2(outMLE$par)
        }
    }
    ##
    mod_hat <- get(paste0("mod_hat", which.min(ricker_aic_df[j, -1]) - 1))
    par.est <- sqrt(exp(outMLE$par)) # compare with input values of 0.2 and 0.1
    LL.est <- outMLE$value
    conv.est <- outMLE$convergence
    ## filter and smooth
    filt <- dlmFilter(y, mod_hat)
    smooth <- dlmSmooth(filt)
    ##
    ahat <- dropFirst(smooth$s[,1])
    bhat <- dropFirst(smooth$s[,2])
    ## uncertainty
    state_var <- with(smooth, dlmSvd2var(U.S, D.S))
    ase <- dropFirst(unlist(lapply(state_var, function(z){sqrt(z[1,1])})))
    bse <- dropFirst(unlist(lapply(state_var, function(z){sqrt(z[2,2])})))
    amat <- cbind(ahat, ase) %*% beta
    bmat <- cbind(bhat, bse) %*% beta
    qmat <- -amat + log(fec0)
    pmat <- -bmat
    dat$qhat_ricker <- qmat[,1]
    dat$qlwr_ricker <- qmat[,2]
    dat$qupr_ricker <- qmat[,3]
    ##
    dat$phat_ricker <- pmat[,1]
    dat$plwr_ricker <- pmat[,2]
    dat$pupr_ricker <- pmat[,3]
    ## recruitment predictions
    for(nn in new_names){
        timeseries[idx, nn] <- dat[, nn]
    }
}

vars2keep <- c(vars2keep, new_names)
range_mult <- c(0.5, 1.5)
pdf("../tex/figures/Fig7_ricker_comparison.pdf", height = 7, width = 10)
par(mfrow = c(3, 4), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 2))
for(j in 1:length(stocks2fit)){
    tmp <- subset(timeseries, ID == stocks2fit[j])
    tmp <- tmp[, vars2keep]
    ## q
    qmat <- tmp[, grep("^q", names(tmp))]
    qmin <- range_mult[1] * min(tmp[, grep("^qhat", names(tmp))], na.rm = TRUE)
    qmax <- range_mult[2] * max(tmp[, grep("^qhat", names(tmp))], na.rm = TRUE)    
    qrange <- c(qmin, qmax)
    with(tmp, matplot(Year, cbind(qhat_ekf, qlwr_ekf, qupr_ekf), col = "darkgrey", type = "l", lty = c(1, 2, 2), ylim = qrange))
    with(tmp, matlines(Year, cbind(qhat_bssm, qlwr_bssm, qupr_bssm), col = "orange2", type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(qhat_tmb, qlwr_tmb, qupr_tmb), col = "darkblue", type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(qhat_ricker, qlwr_ricker, qupr_ricker), col = "forestgreen", type = "l", lty = c(1, 2, 2)))
    legend("topleft", legend = paste(stocks2fit[j], "q"), bty = "n")
    ##
    pmat <- tmp[, grep("^p", names(tmp))]
    pmat <- tmp[, grep("^p", names(tmp))]
    pmin <- range_mult[1] * min(tmp[, grep("^phat", names(tmp))], na.rm = TRUE)
    pmax <- range_mult[2] * max(tmp[, grep("^phat", names(tmp))], na.rm = TRUE)    
    prange <- c(pmin, pmax)
    ##
    with(tmp, matplot(Year, cbind(phat_ekf, plwr_ekf, pupr_ekf), col = "darkgrey", type = "l", lty = c(1, 2, 2), ylim = prange))
    with(tmp, matlines(Year, cbind(phat_bssm, plwr_bssm, pupr_bssm), col = "orange2", type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(phat_tmb, plwr_tmb, pupr_tmb), col = "darkblue", type = "l", lty = c(1, 2, 2)))
    with(tmp, matlines(Year, cbind(phat_ricker, plwr_ricker, pupr_ricker), col = "forestgreen", type = "l", lty = c(1, 2, 2)))
    legend("topleft", legend = paste(stocks2fit[j], "p"), bty = "n")
}
mtext(side = 1, text = "Year", line = 0.5, outer = TRUE)
mtext(side = 2, text = "Parameter estimate", line = 0.5, outer = TRUE)
mtext(side = 4, text = "Low contrast", at = 0.85, line = 0.25, outer = TRUE)
mtext(side = 4, text = "Medium contrast", at = 0.5, line = 0.25, outer = TRUE)
mtext(side = 4, text = "High contrast", at = 0.17, line = 0.25, outer = TRUE)
dev.off()
