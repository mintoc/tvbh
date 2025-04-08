##---------------------------------------------
## Analyse simulation testing of time-varying BH model
## CM: 27/6/23
## run simulations_run.R first
## base the analyses on those for which all three methods converged
## so analyse convergence first
## then model selection
## then performance metrics
##---------------------------------------------

library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(reshape)

library(TMB)
source("utils.R")

library(rpart)
library(rpart.plot)
library(DALEX)

load("./inference.RData")

m <- length(unique(timeseries$ID))

## overall convergence
methods <- c("ekf", "tmb", "bssm")

idvars <- c("h", "a", "b", "phi0", "S0", "tau", "n", "fec_true", "sde", 
            "sdlnq", "sdlnp", "contrast", "iter", "ID")

##---------------------
## ANALYSE CONVERGENCE
##---------------------

sim_df$tmb_static <- ifelse(is.na(sim_df$tmb_aic_static), 1, 0)
sim_df$tmb_qvary <- ifelse(is.na(sim_df$tmb_aic_qvary), 1, 0)
sim_df$tmb_pvary <- ifelse(is.na(sim_df$tmb_aic_pvary), 1, 0)

sim_df$ekf_static <- ifelse(is.na(sim_df$ekf_aic_static), 1, 0)
sim_df$ekf_qvary <- ifelse(is.na(sim_df$ekf_aic_qvary), 1, 0)
sim_df$ekf_pvary <- ifelse(is.na(sim_df$ekf_aic_pvary), 1, 0)

sim_df$bssm_static <- ifelse(is.na(sim_df$bssm_dic_static), 1, 0)
sim_df$bssm_qvary <- ifelse(is.na(sim_df$bssm_dic_qvary), 1, 0)
sim_df$bssm_pvary <- ifelse(is.na(sim_df$bssm_dic_pvary), 1, 0)

method_vars <- apply(expand.grid(c("ekf", "tmb", "bssm"), c("static", "qvary", "pvary")), 1, paste, collapse = "_")

vars2keep <- c(idvars, method_vars)

lconverge <- melt(sim_df[, vars2keep], id.vars = idvars)

lconverge$method <- unlist(lapply(strsplit(as.character(lconverge$variable), "_"), "[", 1))
lconverge$model <- unlist(lapply(strsplit(as.character(lconverge$variable), "_"), "[", 2))

lconverge$h <- factor(lconverge$h)
lconverge$n <- factor(lconverge$n)
lconverge$sde <- factor(lconverge$sde)
lconverge$sdlnq <- factor(lconverge$sdlnq)
lconverge$sdlnp <- factor(lconverge$sdlnp)

oconv <- aggregate(value ~ method, sum, data = lconverge)

## factors incluencing convergence
rtree <-  rpart(value ~ method + factor(h) + factor(n) + factor(sde) +
                    factor(sdlnq) + factor(sdlnp) +
                    + contrast + model,
                data = lconverge, control = rpart.control(cp = 0.005))

## perumatation based variable importance
library(DALEX)

model_vars <- c("h", "n", "sde", "sdlnq", "sdlnp", "contrast", "method", "model")

explainer_rt <- explain(model = rtree, 
                        data = lconverge[, model_vars], 
                        y = lconverge$value, 
                        label = "Convergence tree")

vip_converge <- model_parts(explainer = explainer_rt, 
                            loss_function = loss_root_mean_square,
                            B = 100)

yrange <- range(rtree$frame$yval)

breaks<-seq(0.9 * yrange[1], 1.1 * yrange[2], length = 101)
greys <- colorRampPalette(c("white", "grey"))
cols <- greys(100)

col_index<-cut(rtree$frame$yval, breaks)

## plot
pdf("../tex/figures/FigS1_convergence_tree.pdf", height = 7, width = 8)
print(plot(vip_converge))
prp(rtree, box.col = cols[col_index], faclen=0, varlen=0, do.par=FALSE, left=T, type=4, clip.right.labs=FALSE)
dev.off()

## IDs that converged
ids2keep <- sim_df$ID[rowSums(sim_df[, method_vars]) == 0]


##---------------------
## TRUE MODEL RECOVERY
##---------------------
## for each method and information criteria plot the confusion matrix and get the model choice accuracy
## those that converged for all models and methods
sim_df_conv <- sim_df[sim_df$ID %in% ids2keep, ]

sim_df_conv$true <- NA
sim_df_conv$true[sim_df_conv$sdlnq == 0 & sim_df_conv$sdlnp == 0] <- "true_static"
sim_df_conv$true[sim_df_conv$sdlnq == 0 & sim_df_conv$sdlnp == 0.1] <- "true_pvary-low"
sim_df_conv$true[sim_df_conv$sdlnq == 0 & sim_df_conv$sdlnp == 0.2] <- "true_pvary-high"
sim_df_conv$true[sim_df_conv$sdlnq == 0.1 & sim_df_conv$sdlnp == 0] <- "true_qvary-low"
sim_df_conv$true[sim_df_conv$sdlnq == 0.2 & sim_df_conv$sdlnp == 0] <- "true_qvary-high"


## show correct cells
correct <- data.frame(xmin = c(0.5, 1.5, 2.5, 3.5, 4.5),
                      xmax = c(1.5, 2.5, 3.5, 4.5, 5.5),
                      ymin = c(0.5, 1.5, 1.5, 2.5, 2.5),
                      ymax = c(1.5, 2.5, 2.5, 3.5, 3.5),
                      true = "static", ic = "static", Freq = 0)

## run first with all variables to find out what is important
accuracy_df <- NULL

for(method in c("ekf", "tmb", "bssm")){
    if(method %in% c("ekf", "tmb")){
        ics <- c("aic", "bic")
    }else{
        ics <- "dic"
    }
    for(ic in ics){
        ## confusion
        ic_idx <- grep(paste0("^", method, "_", ic), names(sim_df_conv))
        ic_names <- c(names(sim_df_conv)[ic_idx], "Unconverged")
        idx <- apply(sim_df_conv[, ic_idx], 1, function(z){
            if(all(is.na(z))){
                4
            }else{
                which.min(z)
            }
        })
        sim_df_conv$method_ic <- ic_names[idx]
        ## conf <- with(sim_df_conv, table(n, sde, contrast, true, aic))
        conf <- with(sim_df_conv, table(n, sde, contrast, true, method_ic, sdlnq, sdlnp))
        conf_df <- as.data.frame(conf)
        conf_df$contrast <- factor(as.character(conf_df$contrast), levels = c("upper", "lower", "full"))
        names(conf_df)[names(conf_df) == "n"] <- "nobs"
        ##
        true_df <- aggregate(Freq ~ true + sde + sdlnq + sdlnp + contrast + nobs, FUN = sum, data = conf_df)
        true_df <- subset(true_df, Freq > 0)
        ##true_df <- aggregate(Freq ~ true + sde + contrast + nobs + sdlnsar + sdlnR0, FUN = sum, data = conf_df)
        names(true_df)[names(true_df) == "Freq"] <- "Total"
        conf_df <- merge(true_df, conf_df)
        conf_df$percent <- with(conf_df, Freq/Total * 100)
        conf_df$label <- with(conf_df, paste0(round(Freq/Total, 2) * 100, "%\n(", Freq, ")"))
        conf_df$label[conf_df$Freq == 0] <- NA
        conf_df$true <- gsub("true_", "", conf_df$true)
        conf_df$ic <- gsub(paste(method, ic, "", sep = "_"), "", conf_df$method_ic)
        ##conf_df <- droplevels(conf_df)
        conf_df$true <- factor(conf_df$true, levels = c("static", "qvary-low", "qvary-high", "pvary-low", "pvary-high"))
        conf_df$ic <- factor(conf_df$ic, levels = c("static", "qvary", "pvary", "qpvary", "Unconverged"))
        ## plot
        ## 40 years
        p40 <- ggplot(subset(conf_df, nobs == 40), aes(x = true, y = ic, fill = percent)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "red", na.value = "white") +
            geom_text(aes(label = label)) +
            geom_rect(data = correct, fill=NA, colour = "black", aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + 
            theme(legend.position = "none") +
            ##facet_grid(sde ~ contrast, labeller = "label_both") +
            facet_grid(contrast ~ sde, labeller = "label_both") +
            xlab("True model") +
            ylab(paste(toupper(method), toupper(ic), "selected model")) +
            ggtitle(paste("40 year timeseries", toupper(method), toupper(ic)))
        ## 80 years
        p80 <- ggplot(subset(conf_df, nobs == 80), aes(x = true, y = ic, fill = percent)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "red", na.value = "white") +
            geom_text(aes(label = label)) +
            geom_rect(data = correct, fill=NA, colour = "black", aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + 
            theme(legend.position = "none") +
            ##facet_grid(sde ~ contrast, labeller = "label_both") +
            facet_grid(contrast ~ sde, labeller = "label_both") +
            xlab("True model") +
            ylab(paste(toupper(method), toupper(ic), "selected model")) +
            ggtitle(paste("80 year timeseries", toupper(method), toupper(ic)))
        ## plot to file
        pdf(paste0(paste("../tex/figures/model_selection", method, ic, sep = "_"), ".pdf"), height = 12, width = 8)
        print(p40)
        print(p80)
        dev.off()
        ## accuracy
        ## EKF AIC model accuracy
        sim_df_conv$method_ic <- gsub(paste(method, ic, "", sep = "_"), "", sim_df_conv$method_ic)
        sim_df_conv$correct <- 0
        sim_df_conv$correct[with(sim_df_conv, true == "true_static" & method_ic == "static")] <- 1
        sim_df_conv$correct[with(sim_df_conv, true %in% c("true_qvary-low", "true_qvary-high") & method_ic == "qvary")] <- 1
        sim_df_conv$correct[with(sim_df_conv, true %in% c("true_pvary-low", "true_pvary-high") & method_ic == "pvary")] <- 1
        ##aggregate(correct ~ contrast + sde + n, FUN = function(z){sum(z)/length(z)}, data = sim_df_conv)
        tmp_accuracy <- aggregate(correct ~ sdlnq + sdlnp + sde + contrast + n + h, FUN = function(z){sum(z)/length(z)}, data = sim_df_conv)
        tmp_accuracy$method <- method
        tmp_accuracy$ic <- ic
        accuracy_df <- rbind(accuracy_df, tmp_accuracy)
    }
}

aic_dic <- subset(accuracy_df, ic %in% c("aic", "dic"))

rtree <-  rpart(correct ~ method + factor(h) + factor(n) + factor(sde) +
                    factor(sdlnq) + factor(sdlnp) +
                    + contrast,
                data = aic_dic, control = rpart.control(cp = 0.01))

## perumatation based variable importance
library(DALEX)

model_vars <- c("h", "n", "sde", "sdlnq", "sdlnp", "contrast", "method", "ic")

explainer_rt <- explain(model = rtree, 
                        data = aic_dic[, model_vars], 
                        y = aic_dic$correct, 
                        label = "True model recovery tree")

vip_converge <- model_parts(explainer = explainer_rt, 
                            loss_function = loss_root_mean_square,
                            B = 100)

plot(vip_converge)

yrange <- range(rtree$frame$yval)

breaks<-seq(0.9 * yrange[1], 1.1 * yrange[2], length = 101)
greys <- colorRampPalette(c("white", "grey"))
cols <- greys(100)

col_index<-cut(rtree$frame$yval, breaks)

pdf("../tex/figures/FigS2_true_model_recovery.pdf", height = 7, width = 8)
print(plot(vip_converge))
prp(rtree, box.col = cols[col_index], faclen=0, varlen=0, do.par=FALSE, left=T, type=4, clip.right.labs=FALSE)
dev.off()

## now re-run with contrast, sdlnp, sdlnq, method and n - most important variables

accuracy_df <- NULL

for(method in c("ekf", "tmb", "bssm")){
    if(method %in% c("ekf", "tmb")){
        ics <- c("aic", "bic")
    }else{
        ics <- "dic"
    }
    for(ic in ics){
        ## confusion
        ic_idx <- grep(paste0("^", method, "_", ic), names(sim_df_conv))
        ic_names <- c(names(sim_df_conv)[ic_idx], "Unconverged")
        idx <- apply(sim_df_conv[, ic_idx], 1, function(z){
            if(all(is.na(z))){
                4
            }else{
                which.min(z)
            }
        })
        sim_df_conv$method_ic <- ic_names[idx]
        ## conf <- with(sim_df_conv, table(n, sde, contrast, true, aic))
        conf <- with(sim_df_conv, table(n, sde, contrast, true, method_ic, sdlnq, sdlnp))
        conf_df <- as.data.frame(conf)
        conf_df$contrast <- factor(as.character(conf_df$contrast), levels = c("upper", "lower", "full"))
        names(conf_df)[names(conf_df) == "n"] <- "nobs"
        ##
        true_df <- aggregate(Freq ~ true + sde + sdlnq + sdlnp + contrast + nobs, FUN = sum, data = conf_df)
        true_df <- subset(true_df, Freq > 0)
        names(true_df)[names(true_df) == "Freq"] <- "Total"
        conf_df <- merge(true_df, conf_df)
        conf_df$percent <- with(conf_df, Freq/Total * 100)
        conf_df$true <- gsub("true_", "", conf_df$true)
        conf_df$ic <- gsub(paste(method, ic, "", sep = "_"), "", conf_df$method_ic)
        ##conf_df <- droplevels(conf_df)
        conf_df$true <- factor(conf_df$true, levels = c("static", "qvary-low", "qvary-high", "pvary-low", "pvary-high"))
        conf_df$ic <- factor(conf_df$ic, levels = c("static", "qvary", "pvary", "qpvary", "Unconverged"))
        ## accuracy
        ## EKF AIC model accuracy
        sim_df_conv$method_ic <- gsub(paste(method, ic, "", sep = "_"), "", sim_df_conv$method_ic)
        sim_df_conv$correct <- 0
        sim_df_conv$correct[with(sim_df_conv, true == "true_static" & method_ic == "static")] <- 1
        sim_df_conv$correct[with(sim_df_conv, true %in% c("true_qvary-low", "true_qvary-high") & method_ic == "qvary")] <- 1
        sim_df_conv$correct[with(sim_df_conv, true %in% c("true_pvary-low", "true_pvary-high") & method_ic == "pvary")] <- 1
        tmp_accuracy <- aggregate(correct ~ sdlnq + sdlnp + contrast + n, FUN = function(z){sum(z)/length(z)}, data = sim_df_conv)
        tmp_accuracy$method <- method
        tmp_accuracy$ic <- ic
        accuracy_df <- rbind(accuracy_df, tmp_accuracy)
    }
}

aic_dic <- subset(accuracy_df, ic %in% c("aic", "dic"))

aic_dic$contrast <- factor(as.character(aic_dic$contrast), levels = c("upper", "lower", "full"))

aic_dic$method <- toupper(aic_dic$method)
aic_dic$method[aic_dic$method == "EKF"] <- "eKF"
aic_dic$ic <- toupper(aic_dic$ic)

aic_dic$method_ic <- with(aic_dic, paste(method, ic, sep = "-"))

aic_dic$method_ic <- factor(aic_dic$method_ic, levels = c("eKF-AIC", "BSSM-DIC", "TMB-AIC"))

pdf("../tex/figures/true_model_recovery_summary.pdf", height = 8, width = 6)
ggplot(aic_dic, aes(x = contrast, y = correct, colour = factor(method_ic), group = method_ic)) +
    geom_line() +
    ##geom_point(size = 3) +
    scale_colour_manual("Method-\nInformation criterion", values = c("darkgrey", "orange2", "darkblue")) + 
    facet_grid(sdlnp + sdlnq  ~ n, labeller = label_both) +
    ##scale_y_continuous(lim = c(0, 1)) +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0.5, lty = 3) +
    scale_y_continuous(lim = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2))  +
    xlab("Contrast") +
    ylab("Proportion of times correct model chosen")
dev.off()



##--------------
## EXAMPLE FITS 
##--------------
timeseries <- subset(timeseries, ID %in% ids2keep)

set.seed(101)
ex_ids <- labels <-NULL

sds <- c(0.1, 0.2)
ns <- c(40, 80)
contrasts <- c("upper", "full")

ii <- 1
for(pari in c("q", "p")){
    for(sdi in sds){
        for(ni in ns){
            for(ci in contrasts){
                if(pari == "q"){
                    id <- sample(subset(sim_df_conv, sdlnp == 0 & sdlnq == sdi & sde == 0.2 & n == ni & h == 0.7 & contrast == ci)$ID, 1)
                }else{
                    id <- sample(subset(sim_df_conv, sdlnp == sdi & sdlnq == 0 & sde == 0.2 & n == ni & h == 0.7 & contrast == ci)$ID, 1)
                }
                ex_ids[ii] <- id
                labels[ii] <- paste(pari, sdi, ni, ci, collapse = "_")
                ii <- ii + 1
            }
        }
    }
}

ex <- subset(timeseries, ID %in% ex_ids)

pdf("../tex/figures/Fig4_simulation_example.pdf", height = 9, width = 9)
par(mfrow = c(4, 4), oma = c(2, 2, 4, 1), mar = c(2, 2, 1, 1))
for(i in 1:length(ex_ids)){
    dat <- subset(timeseries, ID == ex_ids[i])
    ## q series
    if(length(grep("q", labels[i])) > 0){
        ##matplot(dat[, c("qt", "qhat_ekf", "qhat_bssm", "qhat_tmb")], type = "l", lty = c(2, 1, 1, 1), col = c("black", "darkgrey", "orange2", "darkblue"), bty = "l")
        ylim <- range(dat[, c("qt", "qhat_ekf", "qhat_bssm", "qhat_tmb")])
        plot(dat$qt, ylim = ylim, bty = "l", col = "slategrey")
        matlines(dat[, c("qhat_ekf", "qhat_bssm", "qhat_tmb")], type = "l", lty = 1, col = c("darkgrey", "orange2", "darkblue"), bty = "l", lwd = 1.5)
        ##points(dat$qt)
        ##legend("topright", legend = labels[i], bty = "n")
    }else{
        ylim <- range(dat[, c("pt", "phat_ekf", "phat_bssm", "phat_tmb")])
        plot(dat$pt, ylim = ylim, bty = "l", col = "slategrey")
        matlines(dat[, c("phat_ekf", "phat_bssm", "phat_tmb")], type = "l", lty = 1, col = c("darkgrey", "orange2", "darkblue"), bty = "l", lwd = 1.5)
        ##legend("topright", legend = labels[i], bty = "n")
    }
    if(i == 4){
        legend("topright", legend = c("True", "eKF", "BSSM", "TMB"), lty = c(NA, 1, 1, 1), pch = c(1, NA, NA, NA), col = c("slategrey", "darkgrey", "orange2", "darkblue"), bty = "n", cex = 1.5, lwd = 1.5)
    }    
    if(i %in% c(1, 3)) mtext(side = 3, text = "Upper", line = 1)
    if(i %in% c(2, 4)) mtext(side = 3, text = "Full", line = 1)
    if(i == 1) mtext(side = 2, line = 2, text = expression(hat(q) * " low variation"))
    if(i == 5) mtext(side = 2, line = 2, text = expression(hat(q) * " high variation"))
    if(i == 9) mtext(side = 2, line = 2, text = expression(hat(p) * " low variation"))
    if(i == 13) mtext(side = 2, line = 2, text = expression(hat(p) * " high variation"))
}
mtext(side = 1, line = 0.5, outer = TRUE, text = "Year")
mtext(side = 3, line = 2, outer = TRUE, at = 0.25, text = "40 year timeseries")
mtext(side = 3, line = 2, outer = TRUE, at = 0.75, text = "80 year timeseries")
dev.off()


##---------------------
## Performance metrics 
##---------------------
error <- function(pred, true){pred - true}

## relative bias
rbias <- function(pred, true){mean(error(pred, true))/mean(true)}

## relative mean absolute error
rmae <- function(pred, true){
    et <- error(pred, true)
    mean(abs(et)) / mean(true)
}

## relative root mean squared error
rrmse <- function(pred, true){
    et <- error(pred, true)
    rmse <- sqrt(mean(et^2))
    rmse / mean(true)
}

perf <- list()

## container for those that have missing, NaN or Inf values
na_id <- NULL

m <- length(unique(timeseries$ID))

pb <- txtProgressBar(min = 1, max = m, initial = 1, style = 3) 

tmp <- expand.grid(par = c("q", "p"),
                   ##measure = c("hat", "lwr", "upr"),
                   measure = c("hat"),
                   method = c("ekf", "tmb", "bssm"))
vars2check <- paste(paste0(tmp$par, tmp$measure), tmp$method, sep = "_")

for(j in 1:m){
    ##print(j)
    tmp <- subset(timeseries, ID == ids2keep[j])[-1,]  ## removing the first observation as sometimes missing for eKF
    ## q
    qtrue <- tmp$qt
    qperfj <- sapply(methods, function(z){
        qmean <- tmp[, paste0("qhat_", z)]
        ##
        bias <- rbias(qmean, qtrue)
        mae <- rmae(qmean, qtrue)
        rmse <- rrmse(qmean, qtrue)
        return(c(rbias = bias,
                 rmae = mae,
                 rrmse = rmse))
    })
    ## p
    ptrue <- tmp$pt
    pperfj <- sapply(methods, function(z){
        pmean <- tmp[, paste0("phat_", z)]
        ##
        bias <- rbias(pmean, ptrue)
        mae <- rmae(pmean, ptrue)
        rmse <- rrmse(pmean, ptrue)
        return(c(rbias = bias,
                 rmae = mae,
                 rrmse = rmse))
    })
    ## check to see if any non-convergence
    na_check <- any(is.na(tmp[, vars2check])) | any(apply(tmp[, vars2check], 1, is.infinite)) | any(tmp[, vars2check] > 1e2)
    if(na_check){
        na_id <- c(na_id, ids2keep[j])
    }
    qperfj <- t(qperfj)
    pperfj <- t(pperfj)
    ##
    qperfj_df <- as.data.frame(qperfj)
    qperfj_df$par <- "q"
    pperfj_df <- as.data.frame(pperfj)
    pperfj_df$par <- "p"
    ##
    perfj_df <- rbind(qperfj_df, pperfj_df)##, saoperfj_df, R0perfj_df)
    ##
    perfj_df$method <- rep(rownames(qperfj_df), 2)
    rownames(perfj_df) <- NULL
    perfj_df <-
        cbind(tmp[rep(1, length(methods) * 2),idvars], perfj_df)
    perf[[j]] <- perfj_df
    setTxtProgressBar(pb, j)
}

idx <- unlist(lapply(perf, is.data.frame))

perf <- do.call(rbind, perf)

perf$Notes <- NULL

perf$na_bin <- ifelse(perf$ID %in% na_id, 1, 0)

rtree <-  rpart(na_bin ~ factor(h) + factor(sdlnq) +
                    factor(sdlnp) + factor(sde) + contrast + method + n,
                data = perf, control = rpart.control(cp = 0.005))

prp(rtree, faclen=0, varlen=0, do.par=FALSE, left=T, type=4, clip.right.labs=FALSE)

## remove the additional non-convergers
perf_nona <- perf[!perf$ID %in% na_id, ]

perf_nona$method <- factor(perf_nona$method, levels = c("ekf", "tmb", "bssm"))

with(perf_nona, boxplot(rbias ~ method + par, notch = TRUE, ylim = c(-3, 3))); abline(h = 0, lty = 2)
with(perf_nona, boxplot(rmae ~ method + par, notch = TRUE, ylim = c(0, 2))); abline(h = 0, lty = 2)
with(perf_nona, boxplot(rrmse ~ method + par, notch = TRUE, ylim = c(0, 2))); abline(h = 0, lty = 2)

pms <- c("rbias", "rmae", "rrmse")

perf_nona$contrast <- as.factor(perf_nona$contrast)

pdf("../tex/figures/FigS2_performance_metrics.pdf", height = 7, width = 8)
for(pm in pms){
    f <- as.formula(paste(pm, "h + sdlnq + sdlnp + sde + contrast + method + n", sep = "~"))
    for(parj in c("q", "p")){ ## problem with R0
        dat <- subset(perf_nona, par == parj)
        dat$h <- factor(dat$h)
        dat$n <- factor(dat$n)
        dat$sde <- factor(dat$sde)
        dat$sdlnq <- factor(dat$sdlnq)
        dat$sdlnp <- factor(dat$sdlnp)
        ## sometimes very large values - set these to top 1%
        max_pm <- quantile(dat[, pm], prob = 0.99, na.rm = TRUE)
        dat[dat[, pm] > max_pm & !is.na(dat[, pm]), pm] <- max_pm
        rtree <-  rpart(f, data = dat, control = rpart.control(cp = 0))
        if(parj == "R0"){
            rprune <- prune(rtree, cp = 0.001)
        }else{
            rprune <- prune(rtree, cp = 0.005)
        }
        ## DALEX variable importance
        model_vars <- c("h", "n", "sde", "sdlnq", "sdlnp", "contrast", "method")
        explainer_rt <- explain(model = rprune, 
                                data = dat[, model_vars], 
                                y = dat[, pm], 
                                label = paste(pm, parj))
        ##
        vip_converge <- model_parts(explainer = explainer_rt, 
                                    loss_function = loss_root_mean_square,
                                    B = 100)
        ##
        print(plot(vip_converge))
        ##
        if(pm == "rbias"){
            breaks<-seq(-1.1 * max(abs(rprune$frame$yval)), 1.1 * max(abs(rprune$frame$yval)), length=101)    
        }else{
            breaks<-seq(0.9 * min(rprune$frame$yval), 1.1 * max(rprune$frame$yval), length=101)    
        }
        if(pm == "rbias"){
            mypal <- colorRampPalette(c("blue", "white", "red"))
        }else{
            mypal <- colorRampPalette(c("white", "red"))
        }
        cols <- mypal(100)
        col_index<-cut(rprune$frame$yval, breaks)
        ##
        prp(rprune, box.col = cols[col_index], faclen=0, varlen=0, do.par=FALSE, left=T, type=4, clip.right.labs=FALSE, main = paste(parj, pm))
    }
}
dev.off()

perf_nona$par <- factor(perf_nona$par, levels = c("q", "p"))

perf_nona$contrast <- factor(as.character(perf_nona$contrast), levels = c("upper", "lower", "full"))

perf_nona <- subset(perf_nona, par %in% c("q", "p"))

perf_nona$Method <- NA
perf_nona$Method[perf_nona$method == "ekf"] <- "eKF"
perf_nona$Method[perf_nona$method == "tmb"] <- "TMB"
perf_nona$Method[perf_nona$method == "bssm"] <- "BSSM"

perf_nona$Method <- factor(perf_nona$Method, levels = c("eKF", "BSSM", "TMB"))

pdf("../tex/figures/Fig5_rerror_boxplots.pdf", height = 10, width = 8)
ggplot(perf_nona, aes(x = factor(sdlnq), y = rbias, fill = Method)) +
    geom_boxplot(outliers = FALSE, size = 0.3) +
    facet_grid(par + h ~ contrast, scales = "free_y", labeller = "label_both") +
    ##scale_fill_manual(values = c("red", "mediumpurple", "darkorange")) +
    ##scale_fill_manual(values = c("mediumpurple","lightblue", "darkorange")) +
    scale_fill_manual(values = c("lightgrey", "white", "darkgrey")) +
    coord_cartesian(ylim = c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(legend.position = "bottom") +
    xlab("Standard deviation of ln(q_y) random walk process") +
    ylab("Relative error")
dev.off()

##--------------------
## PERFORMANCE TABLES
##--------------------
## large table to look up
library(reshape)

perf_nona$na_bin <- NULL
perf_nona$n.1 <- NULL
perf_nona$method <- NULL

perf_nona_long <- melt(perf_nona, id.vars = c(idvars, "par", "Method"))

## main variables
vars2keep <- c("h", "contrast", "Method", "par", "sde", "sdlnq", "variable", "value")

## q
perf_summary_q <- cast(subset(perf_nona_long[, vars2keep], par == "q"),
                       h + contrast + sde + sdlnq ~ variable + Method, function(zz){round(median(zz), 2)})

## remove duplicated id variables 
perf_summary_q$h0 <- perf_summary_q$h
perf_summary_q$contrast0 <- perf_summary_q$contrast
perf_summary_q$sde0 <- perf_summary_q$sde
for(i in 2:nrow(perf_summary_q)){
    if(perf_summary_q$contrast0[i] == perf_summary_q$contrast0[i-1]){
        perf_summary_q$contrast[i] <- ""
    }
    if(perf_summary_q$h0[i] == perf_summary_q$h0[i-1]){
        perf_summary_q$h[i] <- ""
    }
    if(perf_summary_q$sde0[i] == perf_summary_q$sde0[i-1]){
        perf_summary_q$sde[i] <- ""
    }    
}

perf_summary_q$h0 <- perf_summary_q$contrast0 <- perf_summary_q$sde0 <- NULL

library(xtable)

qtab <- xtable(perf_summary_q)

print.xtable(qtab, file = "../tex/tables/simulation_performance_q.tex", include.rownames = FALSE)

## p
perf_summary_p <- cast(subset(perf_nona_long[, vars2keep], par == "p"),
                     h + contrast + sde + sdlnq ~ variable + Method, function(zz){round(median(zz), 2)})

## remove duplicated id variables 
perf_summary_p$h0 <- perf_summary_p$h
perf_summary_p$contrast0 <- perf_summary_p$contrast
perf_summary_p$sde0 <- perf_summary_p$sde
for(i in 2:nrow(perf_summary_p)){
    if(perf_summary_p$contrast0[i] == perf_summary_p$contrast0[i-1]){
        perf_summary_p$contrast[i] <- ""
    }
    if(perf_summary_p$h0[i] == perf_summary_p$h0[i-1]){
        perf_summary_p$h[i] <- ""
    }
    if(perf_summary_p$sde0[i] == perf_summary_p$sde0[i-1]){
        perf_summary_p$sde[i] <- ""
    }    
}

perf_summary_p$h0 <- perf_summary_p$contrast0 <- perf_summary_p$sde0 <- NULL


ptab <- xtable(perf_summary_p)

print.xtable(ptab, file = "../tex/tables/simulation_performance_p.tex", include.rownames = FALSE)


## collapse across variances
## q pretty much unbiased except for upper
overall_summary_q <- cast(subset(perf_nona_long[, vars2keep], par == "q"),
                     contrast + h ~ variable + Method, function(zz){round(median(zz), 2)})

overall_q_long <- melt(as.data.frame(overall_summary_q[, 1:5]), id = c("contrast", "h"))
overall_q_long$par <- "q"

overall_summary_p <- cast(subset(perf_nona_long[, vars2keep], par == "p"),
                     contrast + h ~ variable + Method, function(zz){round(median(zz), 2)})

overall_p_long <- melt(as.data.frame(overall_summary_p[, 1:5]), id = c("contrast", "h"))
overall_p_long$par <- "p"

overall_long <- rbind(overall_q_long, overall_p_long)

overall_long$par <- factor(overall_long$par, levels = c("q", "p"))

pdf("../tex/figures/rerror_summary.pdf", height = 8, width = 6)
ggplot(overall_long, aes(x = h, y = value, fill = variable)) +
    geom_point(pch = 21, size = 3) +
    scale_fill_manual(values = c("lightgrey", "white", "darkgrey")) +
    facet_grid(contrast ~ par) +
    geom_hline(yintercept = 0, lty = 2) +
    theme(legend.position = "bottom")
dev.off()
k
