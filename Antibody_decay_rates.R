################################################################################
# This code calculates antibody decay rates observed                           #
# between May 2020 and June 2021.                                              #
#                                                                              #
# Cite as:                                                                     #
# Lavezzo E. et al, Neutralising reactivity against SARS-CoV-2 Delta           #
# and Omicron variants by vaccination and infection history.                   #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Antibody_decay_rates.R'\n\n")


# set seed for bootstrapping
set.seed(1)

# source scripts
source("./functions_association_antibody_slopes.R")


# read data
library(xlsx)
ds <- read.xlsx("./DB.xlsx", "Vo'_June_2021", header = TRUE)


# formatting ------------------------------------------------------------------#

# select baseline groud truth
ds <- ds[which(ds$Groundtruth == 1), ]

# select unvaccinated subjects
Diasorin_noVax <- ds[which(ds$Vaccinated_at_least_7_days_before_sampling_time == "no"), ]

#select subject with no roche antibody titres doubled in november (no abbott doubled in november )

Roche_noD <- ds[which(ds$Roche_not_doubled == "TRUE"), ]

Diasorin_noD <- Diasorin_noVax[which(Diasorin_noVax$Diasorin_not_doubled == "TRUE"), ]
                

ds[ds=="NA"] = NA


ds$Abbott_quantitative_June_2021 <-
  as.numeric(as.character(ds$Abbott_quantitative_June_2021))

ds$Abbott_semiquantitative_May_2020 <-
  as.numeric(as.character(ds$Abbott_semiquantitative_May_2020))

Roche_noD$Roche_quantitative_May_2020 <-
  as.numeric(as.character(Roche_noD$Roche_quantitative_May_2020))

Roche_noD$Roche_quantitative_June_2021 <-
  as.numeric(as.character(Roche_noD$Roche_quantitative_June_2021))

Diasorin_noD$DiaSorin_S1_S2_quantitative_June_2021 <-
  as.numeric(as.character(Diasorin_noD$DiaSorin_S1_S2_quantitative_June_2021))

Diasorin_noD$DiaSorin_S1_S2_quantitative_May_2020 <-
  as.numeric(as.character(Diasorin_noD$DiaSorin_S1_S2_quantitative_May_2020))




# number of days between serosurveys
days_surveys <- as.numeric(as.Date("06/05/2021", "%m/%d/%Y") -
                             as.Date("05/01/2020", "%m/%d/%Y"))

################################################################################
# Decay rate                                                                   #
################################################################################
Abbott=ds

ds$slope_Abbott <- log(ds$Abbott_quantitative_June_2021/
                         ds$Abbott_semiquantitative_May_2020)/days_surveys

Roche_noD$slope_Roche <- log(Roche_noD$Roche_quantitative_June_2021/
                               Roche_noD$Roche_quantitative_May_2020)/days_surveys

Diasorin_noD$slope_Diasorin <- log(Diasorin_noD$DiaSorin_S1_S2_quantitative_June_2021/
                           Diasorin_noD$DiaSorin_S1_S2_quantitative_May_2020)/days_surveys

################################################################################
# Half life                                                                    #
################################################################################



ds$hl_Abbott <- log(0.5)/ds$slope_Abbott
Diasorin_noD$hl_Diasorin <- log(0.5)/Diasorin_noD$slope_Diasorin
Roche_noD$hl_Roche <- log(0.5)/Roche_noD$slope_Roche

# Abbott
compute_half_life(ds, "Abbott")
bootstrap_hl(ds, "Abbott")

# DiaSorin
compute_half_life(Diasorin_noD, "DiaSorin")
bootstrap_hl(Diasorin_noD, "DiaSorin")

# Roche
compute_half_life(Roche_noD, "Roche")
bootstrap_hl(Roche_noD, "Roche")
