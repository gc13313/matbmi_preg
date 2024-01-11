############################################################################
#                                                                          #
#              This script runs MR for each individual study               #
#                                                                          #
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 

setwd("mypath")


##
## Clear the work environment

rm(list = ls())


##
## Setting repository (UoB)

options(repos = c(CRAN ="http://www.stats.bris.ac.uk/R/"))


##
## Setting digits

options(digits = 10)


##
## Library directory

.libPaths()


##
## Updating packages 

#update.packages(ask = "FALSE")


##
## Install packages

#install.packages(c("data.table", "purrr", "devtools", "xlsx", "dplyr"))


##
## Load libraries

library(data.table)
library(dplyr)
library(purrr)
library(haven)
library(TwoSampleMR)
library(MRInstruments)


############################################################################
#  Run MR and sensitivity analyses (study-specific and pooled estimates)   #
############################################################################

##
## Import study specific outcome information

load("./FINAL_METANALYSIS_2021/data/snp-trait.RData")


##
## Harmonise SNP-exposure and SNP-outcome data

mr.all_dat <- map(all_outdat2, ~harmonise_data(giant_expdat, ., action = 2)) %>%
  bind_rows %>%
  mutate(
    snp32 = ifelse(SNP %in% snps32$snp, T, F),
    snp97 = ifelse(SNP %in% snps97$snp, T, F)
  )


##
## Extract study, outcome, model for each id

idout <- stringr::str_split(unique(mr.all_dat$id.outcome), "_x_", simplify = TRUE) %>% 
  data.table %>% 
  mutate(id.outcome = unique(mr.all_dat$id.outcome)) %>%
  rename(study = V1, outcome = V2, model = V3)

mr.all_dat <- merge(mr.all_dat, idout[, -"outcome"], by = "id.outcome") 


##
## MR methods

mr.all_res <- mr(mr.all_dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")) %>%
  merge(., idout[, -"outcome"], by = "id.outcome") 

write.csv(mr.all_res, file = "./FINAL_METANALYSIS_2021/data/mr_res.csv")


##
## Sensitivity analyses

fix_dat <- filter(mr.all_dat, study == "Pooled (FE)" & model == "unadjusted")

# Between SNP heterogeneity (ie Cochrane's Q)

het <- mr_heterogeneity(fix_dat, method_list="mr_ivw")

# Directional pleiotropy (ie MR-Egger intercept)

int <- mr_pleiotropy_test(fix_dat)


##
## save the workspace

rm(list=setdiff(ls(), c("mr.all_dat", "mr.all_res", "het", "int")))

save.image("./FINAL_METANALYSIS_2021/data/mr_res.RData")