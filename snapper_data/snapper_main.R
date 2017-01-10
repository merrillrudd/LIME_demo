rm(list=ls())
###################################
## Directories
###################################

## set your working directory
main_dir <- "C:\\Git_Projects\\LIME_demo\\snapper_data"

setwd(main_dir)
source("R_functions\\functions.R")

###################################
## Get data for LIME
###################################
## snapper data
lg <- read.csv(file.path(main_dir, "cr_snapper_filtered.csv"))

#### length frequency data
## life history list
## fixed values required: vbk, linf, lwa, lwb, M50
## fixed values that have defaults included: binwidth, CVlen, M, F1
## starting values that will be estimated parameters: S50, SigmaR
cr_lh <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=3, selex_input="age", M50=34, maturity_input="length", SigmaR=0.4, M=0.43, F1=0.34)

lf <- length_frequency(binwidth=1, linf=cr_lh$linf, lmat=cr_lh$Lmat, data=lg, plot=FALSE, weight=TRUE)

## some other details about the data
other <- manipulate_snap(lg, TRUE)

saveRDS(lf, "Snapper_LenFreq.rds")
saveRDS(other$obs_num, "Snapper_TotalObs.rds")
saveRDS(other$obs_day, "Snapper_Days.rds")

