rm(list=ls())
###################################
## Load packages
###################################

## if you are connected to internet, you can load the latest version of the website from github with the devtools package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)

## if you have LIME downloaded on your computer, you can skip to loading it in R
library(LIME)

## compare to LB-SPR
library(LBSPR)

###################################
## Directories
###################################

## set your working directory
main_dir <- "C:\\Git_Projects\\LIME_demo"

## data directory
data_dir <- file.path(main_dir, "data")

## results directory
res_dir <- file.path(main_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

########################################
## Example life history: snapper
########################################

## life history list
## fixed values required: vbk, linf, lwa, lwb, M50
## fixed values that have defaults included: binwidth, CVlen, M, F1
## starting values that will be estimated parameters: S50, SigmaR
cr_lh <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.3, M=0.43, F1=0.34)

par(mfrow=c(2,2))
plot(cr_lh$L_a, ylim=c(0, max(cr_lh$L_a)*1.1), lwd=3, type="l", xlab="Age", ylab="Length")
plot(cr_lh$W_a, ylim=c(0, max(cr_lh$W_a)*1.1), lwd=3, type="l", xlab="Age", ylab="Weight")
plot(cr_lh$Mat_a, ylim=c(0, max(cr_lh$Mat_a)*1.1), lwd=3, type="l", xlab="Age", ylab="Proportion mature")
plot(cr_lh$S_a, ylim=c(0, max(cr_lh$S_a)*1.1), lwd=3, type="l", xlab="Age", ylab="Selectivity")

########################################
## Simulate a population in equilibrium
########################################
cr_lhsim <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.01, SigmaF=0.01, M=0.43, F1=0.34)

Nyears_sim <- 9
Nobs_length <- 1000

simdata <- sim_pop(lh=cr_lhsim, Nyears=Nyears_sim, Fdynamics="Constant", Rdynamics="Constant", Nyears_comp=Nyears_sim, comp_sample=rep(Nobs_length,Nyears_sim), nburn=50, seed=123, modname="simsnap")

par(mfrow=c(2,2))
plot(simdata$F_t, ylim=c(0, max(simdata$F_t)*1.2), type="l", lwd=3, xlab="Year", ylab="Fishing mortality")
plot(simdata$R_t, ylim=c(0, max(simdata$R_t)*1.2), type="l", lwd=3, xlab="Year", ylab="Recruitment")
plot(simdata$ML_t, ylim=c(0, max(simdata$ML_t)*1.2), type="l", lwd=3, xlab="Year", ylab="Mean length in the catch")
plot(simdata$SPR_t, ylim=c(0, max(simdata$SPR_t)*1.2), type="l", lwd=3, xlab="Year", ylab="Spawning potential ratio (SPR)")

########################################
## Format simulated data
########################################
input_simdata <- list("years"=1:Nyears_sim, "LF"=simdata$LF, "obs_per_year"=rep(Nobs_length,Nyears_sim))

########################################
## Run LIME model for simulated data
########################################
results_sim <- run_LIME(modpath=NULL, write=FALSE, lh=cr_lhsim, input_data=input_simdata, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, fix_f=0, simulation=FALSE, param_adjust="SigmaF", val_adjust=0.1, f_true=FALSE, fix_param=FALSE)

########################################
## Run LBSPR model for simulated data
########################################
LB_pars <- new("LB_pars")
LB_pars@Linf <- cr_lhsim$linf
LB_pars@CVLinf <- cr_lhsim$CVlen
LB_pars@L50 <- cr_lhsim$ML50 
LB_pars@L95 <- cr_lhsim$ML95
LB_pars@MK <- cr_lhsim$M/cr_lhsim$vbk         

LB_pars@SL50 <- cr_lhsim$SL50 
LB_pars@SL95 <- cr_lhsim$SL95
LB_pars@BinWidth <- cr_lhsim$binwidth
LB_pars@Walpha <- cr_lhsim$lwa
LB_pars@Wbeta <- cr_lhsim$lwb
LB_pars@Steepness <- min(0.99, cr_lhsim$h)
LB_pars@R0 <- cr_lhsim$R0

LB_pars@Species <- "snapper"

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- seq(cr_lhsim$binwidth/2, length=ncol(input_simdata$LF), by=cr_lhsim$binwidth)
LB_lengths@LData <- t(input_simdata$LF)
LB_lengths@Years <- 1:Nyears_sim
LB_lengths@NYears <- Nyears_sim

lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel"))
lbspr_res2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))


LBSPR_outs <- list()
LBSPR_outs$pLF <- lbspr_res@pLCatch
LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]

LBSPR_outs2 <- list()
LBSPR_outs2$pLF <- lbspr_res2@pLCatch
LBSPR_outs2$SL50 <- lbspr_res2@Ests[,"SL50"]
LBSPR_outs2$SL95 <- lbspr_res2@Ests[,"SL95"]
LBSPR_outs2$FM <- lbspr_res2@Ests[,"FM"]
LBSPR_outs2$SPR <- lbspr_res2@Ests[,"SPR"]

###########################################
## Compare results
###########################################
results_sim$df
rep <- results_sim$Report
sdrep <- results_sim$Sdreport
der <- results_sim$Derived
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$F_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$FM*cr_lh$M, col="red", lwd=2)
      lines(LBSPR_outs2$FM*cr_lh$M, col="orange", lwd=2)
      lines(simdata$F_t, col="black", lwd=4)

      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$R_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(rep(1,length(simdata$R_t)), col="red", lwd=2)
      lines(simdata$R_t, col="black", lwd=4)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0, max(rep$SPR_t)*3), xlab="Year", ylab="SPR", xaxs="i", yaxs="i")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,max(rep$SPR_t)*3, max(rep$SPR_t)*3), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$SPR, col="red", lwd=2)
      lines(LBSPR_outs2$SPR, col="orange", lwd=2)
      lines(simdata$SPR_t, col="black", lwd=4)
      
      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,3), xaxs="i", yaxs="i", xlab="SPR", ylab="F/F30")
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,3,3), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(simdata, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.na(F30)==FALSE){
      	points(x=simdata$SPR_t[length(simdata$SPR_t)], y=simdata$F_t[length(simdata$F_t)]/F30, col="black", pch=19, cex=2)
      	points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*cr_lhsim$M)/F30, col="red", pch=19, cex=2)
      	points(x=LBSPR_outs2$SPR[length(LBSPR_outs2$SPR)], y=(LBSPR_outs2$FM[length(LBSPR_outs2$FM)]*cr_lhsim$M)/F30, col="orange", pch=19, cex=2)

      }

########################################
## Simulate a population with variability
########################################
cr_lhsim <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.7, SigmaF=0.1, M=0.43, F1=0.34)

Nyears_sim <- 9
Nobs_length <- 1000

simdata <- sim_pop(lh=cr_lhsim, Nyears=Nyears_sim, Fdynamics="Constant", Rdynamics="Constant", Nyears_comp=Nyears_sim, comp_sample=rep(Nobs_length,Nyears_sim), nburn=50, seed=123, modname="simsnap")

par(mfrow=c(2,2))
plot(simdata$F_t, ylim=c(0, max(simdata$F_t)*1.2), type="l", lwd=3)
plot(simdata$R_t, ylim=c(0, max(simdata$R_t)*1.2), type="l", lwd=3)
plot(simdata$ML_t, ylim=c(0, max(simdata$ML_t)*1.2), type="l", lwd=3)
plot(simdata$SPR_t, ylim=c(0, max(simdata$SPR_t)*1.2), type="l", lwd=3)

########################################
## Format simulated data
########################################
input_simdata <- list("years"=1:Nyears_sim, "LF"=simdata$LF, "obs_per_year"=rep(Nobs_length,Nyears_sim))

########################################
## Run LIME model for simulated data
########################################
results_sim <- run_LIME(modpath=NULL, write=FALSE, lh=cr_lhsim, input_data=input_simdata, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, fix_f=0, simulation=FALSE, param_adjust="SigmaF", val_adjust=0.1, f_true=FALSE, fix_param=FALSE)

########################################
## Run LBSPR model for simulated data
########################################
LB_pars <- new("LB_pars")
LB_pars@Linf <- cr_lhsim$linf
LB_pars@CVLinf <- cr_lhsim$CVlen
LB_pars@L50 <- cr_lhsim$ML50 
LB_pars@L95 <- cr_lhsim$ML95
LB_pars@MK <- cr_lhsim$M/cr_lhsim$vbk         

LB_pars@SL50 <- cr_lhsim$SL50 
LB_pars@SL95 <- cr_lhsim$SL95
LB_pars@BinWidth <- cr_lhsim$binwidth
LB_pars@Walpha <- cr_lhsim$lwa
LB_pars@Wbeta <- cr_lhsim$lwb
LB_pars@Steepness <- min(0.99, cr_lhsim$h)
LB_pars@R0 <- cr_lhsim$R0

LB_pars@Species <- "snapper"

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- seq(cr_lhsim$binwidth/2, length=ncol(input_simdata$LF), by=cr_lhsim$binwidth)
LB_lengths@LData <- t(input_simdata$LF)
LB_lengths@Years <- 1:Nyears_sim
LB_lengths@NYears <- Nyears_sim


lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel"))
lbspr_res2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))


LBSPR_outs <- list()
LBSPR_outs$pLF <- lbspr_res@pLCatch
LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]

LBSPR_outs2 <- list()
LBSPR_outs2$pLF <- lbspr_res2@pLCatch
LBSPR_outs2$SL50 <- lbspr_res2@Ests[,"SL50"]
LBSPR_outs2$SL95 <- lbspr_res2@Ests[,"SL95"]
LBSPR_outs2$FM <- lbspr_res2@Ests[,"FM"]
LBSPR_outs2$SPR <- lbspr_res2@Ests[,"SPR"]

###########################################
## Compare results
###########################################
results_sim$df
rep <- results_sim$Report
sdrep <- results_sim$Sdreport
der <- results_sim$Derived
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$F_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$FM*cr_lh$M, col="red", lwd=2)
      lines(LBSPR_outs2$FM*cr_lh$M, col="orange", lwd=2)
      lines(simdata$F_t, col="black", lwd=4)

      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$R_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(rep(1,length(simdata$R_t)), col="red", lwd=2)
      lines(simdata$R_t, col="black", lwd=4)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0, max(rep$SPR_t)*3), xlab="Year", ylab="SPR", xaxs="i", yaxs="i")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,max(rep$SPR_t)*3, max(rep$SPR_t)*3), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$SPR, col="red", lwd=2)
      lines(LBSPR_outs2$SPR, col="orange", lwd=2)
      lines(simdata$SPR_t, col="black", lwd=4)
      
      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,3), xaxs="i", yaxs="i", xlab="SPR", ylab="F/F30")
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,3,3), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(simdata, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.na(F30)==FALSE){
      	points(x=simdata$SPR_t[length(simdata$SPR_t)], y=simdata$F_t[length(simdata$F_t)]/F30, col="black", pch=19, cex=2)
      	points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*cr_lhsim$M)/F30, col="red", pch=19, cex=2)
      	points(x=LBSPR_outs2$SPR[length(LBSPR_outs2$SPR)], y=(LBSPR_outs2$FM[length(LBSPR_outs2$FM)]*cr_lhsim$M)/F30, col="orange", pch=19, cex=2)

      }

#####################################################
## Simulate a population without variability and with trends
#####################################################
cr_lhsim <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.01, SigmaF=0.01, M=0.43, F1=0.34)

Nyears_sim <- 9
Nobs_length <- 1000

simdata <- sim_pop(lh=cr_lhsim, Nyears=Nyears_sim, Fdynamics="Increasing", Rdynamics="Pulsed", Nyears_comp=Nyears_sim, comp_sample=rep(Nobs_length,Nyears_sim), nburn=50, seed=123, modname="simsnap")

par(mfrow=c(2,2))
plot(simdata$F_t, ylim=c(0, max(simdata$F_t)*1.2), type="l", lwd=3)
plot(simdata$R_t, ylim=c(0, max(simdata$R_t)*1.2), type="l", lwd=3)
plot(simdata$ML_t, ylim=c(0, max(simdata$ML_t)*1.2), type="l", lwd=3)
plot(simdata$SPR_t, ylim=c(0, max(simdata$SPR_t)*1.2), type="l", lwd=3)

########################################
## Format simulated data
########################################
input_simdata <- list("years"=1:Nyears_sim, "LF"=simdata$LF, "obs_per_year"=rep(Nobs_length,Nyears_sim))

########################################
## Run LIME model for simulated data
########################################
results_sim <- run_LIME(modpath=NULL, write=FALSE, lh=cr_lhsim, input_data=input_simdata, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, fix_f=0, simulation=FALSE, param_adjust="SigmaF", val_adjust=0.1, f_true=FALSE, fix_param=FALSE)

########################################
## Run LBSPR model for simulated data
########################################
LB_pars <- new("LB_pars")
LB_pars@Linf <- cr_lhsim$linf
LB_pars@CVLinf <- cr_lhsim$CVlen
LB_pars@L50 <- cr_lhsim$ML50 
LB_pars@L95 <- cr_lhsim$ML95
LB_pars@MK <- cr_lhsim$M/cr_lhsim$vbk         

LB_pars@SL50 <- cr_lhsim$SL50 
LB_pars@SL95 <- cr_lhsim$SL95
LB_pars@BinWidth <- cr_lhsim$binwidth
LB_pars@Walpha <- cr_lhsim$lwa
LB_pars@Wbeta <- cr_lhsim$lwb
LB_pars@Steepness <- min(0.99, cr_lhsim$h)
LB_pars@R0 <- cr_lhsim$R0

LB_pars@Species <- "snapper"

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- seq(cr_lhsim$binwidth/2, length=ncol(input_simdata$LF), by=cr_lhsim$binwidth)
LB_lengths@LData <- t(input_simdata$LF)
LB_lengths@Years <- 1:Nyears_sim
LB_lengths@NYears <- Nyears_sim


lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel"))
lbspr_res2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))


LBSPR_outs <- list()
LBSPR_outs$pLF <- lbspr_res@pLCatch
LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]

LBSPR_outs2 <- list()
LBSPR_outs2$pLF <- lbspr_res2@pLCatch
LBSPR_outs2$SL50 <- lbspr_res2@Ests[,"SL50"]
LBSPR_outs2$SL95 <- lbspr_res2@Ests[,"SL95"]
LBSPR_outs2$FM <- lbspr_res2@Ests[,"FM"]
LBSPR_outs2$SPR <- lbspr_res2@Ests[,"SPR"]

###########################################
## Compare results
###########################################
results_sim$df
rep <- results_sim$Report
sdrep <- results_sim$Sdreport
der <- results_sim$Derived
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$F_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$FM*cr_lh$M, col="red", lwd=2)
      lines(LBSPR_outs2$FM*cr_lh$M, col="orange", lwd=2)
      lines(simdata$F_t, col="black", lwd=4)

      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$R_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(rep(1,length(simdata$R_t)), col="red", lwd=2)
      lines(simdata$R_t, col="black", lwd=4)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0, max(rep$SPR_t)*3), xlab="Year", ylab="SPR", xaxs="i", yaxs="i")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,max(rep$SPR_t)*3, max(rep$SPR_t)*3), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$SPR, col="red", lwd=2)
      lines(LBSPR_outs2$SPR, col="orange", lwd=2)
      lines(simdata$SPR_t, col="black", lwd=4)
      
      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,3), xaxs="i", yaxs="i", xlab="SPR", ylab="F/F30")
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,3,3), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(simdata, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.na(F30)==FALSE){
      	points(x=simdata$SPR_t[length(simdata$SPR_t)], y=simdata$F_t[length(simdata$F_t)]/F30, col="black", pch=19, cex=2)
      	points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*cr_lhsim$M)/F30, col="red", pch=19, cex=2)
      	points(x=LBSPR_outs2$SPR[length(LBSPR_outs2$SPR)], y=(LBSPR_outs2$FM[length(LBSPR_outs2$FM)]*cr_lhsim$M)/F30, col="orange", pch=19, cex=2)

      }

#####################################################
## Simulate a population with variability and trends
#####################################################
cr_lhsim <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.7, SigmaF=0.1, M=0.43, F1=0.34)

Nyears_sim <- 9
Nobs_length <- 1000

simdata <- sim_pop(lh=cr_lhsim, Nyears=Nyears_sim, Fdynamics="Increasing", Rdynamics="Pulsed", Nyears_comp=Nyears_sim, comp_sample=rep(Nobs_length,Nyears_sim), nburn=50, seed=123, modname="simsnap")

par(mfrow=c(2,2))
plot(simdata$F_t, ylim=c(0, max(simdata$F_t)*1.2), type="l", lwd=3)
plot(simdata$R_t, ylim=c(0, max(simdata$R_t)*1.2), type="l", lwd=3)
plot(simdata$ML_t, ylim=c(0, max(simdata$ML_t)*1.2), type="l", lwd=3)
plot(simdata$SPR_t, ylim=c(0, max(simdata$SPR_t)*1.2), type="l", lwd=3)

########################################
## Format simulated data
########################################
input_simdata <- list("years"=1:Nyears_sim, "LF"=simdata$LF, "obs_per_year"=rep(Nobs_length,Nyears_sim))

########################################
## Run LIME model for simulated data
########################################
results_sim <- run_LIME(modpath=NULL, write=FALSE, lh=cr_lhsim, input_data=input_simdata, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, fix_f=0, simulation=FALSE, param_adjust="SigmaF", val_adjust=0.1, f_true=FALSE, fix_param=FALSE)

########################################
## Run LBSPR model for simulated data
########################################
LB_pars <- new("LB_pars")
LB_pars@Linf <- cr_lhsim$linf
LB_pars@CVLinf <- cr_lhsim$CVlen
LB_pars@L50 <- cr_lhsim$ML50 
LB_pars@L95 <- cr_lhsim$ML95
LB_pars@MK <- cr_lhsim$M/cr_lhsim$vbk         

LB_pars@SL50 <- cr_lhsim$SL50 
LB_pars@SL95 <- cr_lhsim$SL95
LB_pars@BinWidth <- cr_lhsim$binwidth
LB_pars@Walpha <- cr_lhsim$lwa
LB_pars@Wbeta <- cr_lhsim$lwb
LB_pars@Steepness <- min(0.99, cr_lhsim$h)
LB_pars@R0 <- cr_lhsim$R0

LB_pars@Species <- "snapper"

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- seq(cr_lhsim$binwidth/2, length=ncol(input_simdata$LF), by=cr_lhsim$binwidth)
LB_lengths@LData <- t(input_simdata$LF)
LB_lengths@Years <- 1:Nyears_sim
LB_lengths@NYears <- Nyears_sim


lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel"))
lbspr_res2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))


LBSPR_outs <- list()
LBSPR_outs$pLF <- lbspr_res@pLCatch
LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]

LBSPR_outs2 <- list()
LBSPR_outs2$pLF <- lbspr_res2@pLCatch
LBSPR_outs2$SL50 <- lbspr_res2@Ests[,"SL50"]
LBSPR_outs2$SL95 <- lbspr_res2@Ests[,"SL95"]
LBSPR_outs2$FM <- lbspr_res2@Ests[,"FM"]
LBSPR_outs2$SPR <- lbspr_res2@Ests[,"SPR"]

###########################################
## Compare results
###########################################
results_sim$df
rep <- results_sim$Report
sdrep <- results_sim$Sdreport
der <- results_sim$Derived
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$F_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$FM*cr_lh$M, col="red", lwd=2)
      lines(LBSPR_outs2$FM*cr_lh$M, col="orange", lwd=2)
      lines(simdata$F_t, col="black", lwd=4)

      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$R_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(rep(1,length(simdata$R_t)), col="red", lwd=2)
      lines(simdata$R_t, col="black", lwd=4)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0, max(rep$SPR_t)*3), xlab="Year", ylab="SPR", xaxs="i", yaxs="i")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,max(rep$SPR_t)*3, max(rep$SPR_t)*3), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$SPR, col="red", lwd=2)
      lines(LBSPR_outs2$SPR, col="orange", lwd=2)
      lines(simdata$SPR_t, col="black", lwd=4)
      
      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,3), xaxs="i", yaxs="i", xlab="SPR", ylab="F/F30")
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,3,3), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(simdata, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.na(F30)==FALSE){
      	points(x=simdata$SPR_t[length(simdata$SPR_t)], y=simdata$F_t[length(simdata$F_t)]/F30, col="black", pch=19, cex=2)
      	points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*cr_lhsim$M)/F30, col="red", pch=19, cex=2)
      	points(x=LBSPR_outs2$SPR[length(LBSPR_outs2$SPR)], y=(LBSPR_outs2$FM[length(LBSPR_outs2$FM)]*cr_lhsim$M)/F30, col="orange", pch=19, cex=2)

      }

########################################
## Format real data
########################################
lf <- readRDS(file.path(data_dir, "Snapper_LenFreq.rds"))
totals <- readRDS(file.path(data_dir, "Snapper_TotalObs.rds"))
days <- readRDS(file.path(data_dir, "Snapper_Days.rds"))

Years_lf <- as.numeric(rownames(lf))
Years_lf_total <- min(Years_lf):max(Years_lf)
Nyears_model <- length(Years_lf_total) + 0
Years_model <- (Years_lf_total[1]-(Nyears_model - nrow(lf))):(max(Years_lf_total))
obs_per_year <- rep(0, length(Years_model))
	obs_per_year[which(Years_model %in% Years_lf)] <- days
	names(obs_per_year) <- Years_model

write.csv(lf, "Snapper_LenFreq9.csv")
write.csv(obs_per_year, "Snapper_ESS9_days.csv")

input_data <- list("years"=Years_model, "LF"=lf, "obs_per_year"=obs_per_year)

########################################
## Run LIME model
########################################
results <- run_LIME(modpath=NULL, write=FALSE, lh=cr_lh, input_data=input_data, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, fix_f=0, simulation=FALSE, param_adjust="SigmaF", val_adjust=0.1, f_true=FALSE, fix_param=FALSE)

########################################
## Run LBSPR model for simulated data
########################################
LB_pars <- new("LB_pars")
LB_pars@Linf <- cr_lh$linf
LB_pars@CVLinf <- cr_lh$CVlen
LB_pars@L50 <- cr_lh$ML50 
LB_pars@L95 <- cr_lh$ML95
LB_pars@MK <- cr_lh$M/cr_lh$vbk         

LB_pars@SL50 <- cr_lh$SL50 
LB_pars@SL95 <- cr_lh$SL95
LB_pars@BinWidth <- cr_lh$binwidth
LB_pars@Walpha <- cr_lh$lwa
LB_pars@Wbeta <- cr_lh$lwb
LB_pars@Steepness <- min(0.99, cr_lh$h)
LB_pars@R0 <- cr_lh$R0

LB_pars@Species <- "snapper"

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- seq(cr_lh$binwidth/2, length=ncol(input_data$LF), by=cr_lh$binwidth)
LB_lengths@LData <- t(input_data$LF)
LB_lengths@Years <- 1:Nyears_sim
LB_lengths@NYears <- Nyears_sim

lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel"))
lbspr_res2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))


LBSPR_outs <- list()
LBSPR_outs$pLF <- lbspr_res@pLCatch
LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]

LBSPR_outs2 <- list()
LBSPR_outs2$pLF <- lbspr_res2@pLCatch
LBSPR_outs2$SL50 <- lbspr_res2@Ests[,"SL50"]
LBSPR_outs2$SL95 <- lbspr_res2@Ests[,"SL95"]
LBSPR_outs2$FM <- lbspr_res2@Ests[,"FM"]
LBSPR_outs2$SPR <- lbspr_res2@Ests[,"SPR"]

###########################################
## Compare results
###########################################
results$df
rep <- results$Report
sdrep <- results$Sdreport
der <- results$Derived
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$F_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$FM*cr_lh$M, col="red", lwd=2)
      lines(LBSPR_outs2$FM*cr_lh$M, col="orange", lwd=2)

      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(rep$R_t)*3), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(rep(1,length(simdata$R_t)), col="red", lwd=2)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0, max(rep$SPR_t)*3), xlab="Year", ylab="SPR", xaxs="i", yaxs="i")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,max(rep$SPR_t)*3, max(rep$SPR_t)*3), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      lines(LBSPR_outs$SPR, col="red", lwd=2)
      lines(LBSPR_outs2$SPR, col="orange", lwd=2)
      
      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,3), xaxs="i", yaxs="i", xlab="SPR", ylab="F/F30")
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,3,3), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(rep, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.na(F30)==FALSE){
      	points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*cr_lh$M)/F30, col="red", pch=19, cex=2)
      	points(x=LBSPR_outs2$SPR[length(LBSPR_outs2$SPR)], y=(LBSPR_outs2$FM[length(LBSPR_outs2$FM)]*cr_lh$M)/F30, col="orange", pch=19, cex=2)

      }