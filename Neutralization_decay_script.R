# ################################################################################
# # This code calculates the halflife of neutralising antibodies observed        #
# # between May 2020 and June 2021 considering only subjects who were positive   #
# # in May 2020 and not increasing between May 2020 and November 2020 and        #
# # between November 2020 and June 2021.                                         #  
# #                                                                              #
# #                                                                              #
# # Cite as:                                                                     #
# # Lavezzo E et al, Neutralising reactivity against SARS-CoV-2 Delta            #
# # and Omicron variants by vaccination and infection history.                     #
# #                                                                              #
# # Description:                                                                 #
# # See README.md                                                                #
# ################################################################################
# 
cat("\n\n### Running 'scripts/Neutralization_decay_script.R'\n\n")
# 
# # source scripts
source("./functions_association_antibody_slopes.R")
# 

############################################################################
# function to transform categorical neutralization values into continuous values
#############################################################################

sample_FORW.fun = function(x)
{
  test = x/10
  val = 0
  
  if(test != 0)
    if(test < 1) val = sample(x:9, 1, prob=dgeom(x:9, 0.5)) else val = sample(x:(2*x-1), 1, prob=dgeom(x:(2*x-1)-x, 0.5/(2^log(test,2))))
  val
}
###################################################################################
# end of function
###################################################################################


# read data
library(xlsx)
ds <- read.xlsx("C:\\Users\\medcomp\\Desktop\\Omicron\\Scripts\\DB.xlsx", "Vo'_June_2021", header = TRUE)


# select baseline groud truth
ds <- ds[which(ds$Groundtruth == 1), ]

# select not doubled in nov
ds <- ds[which(ds$Neutralisation_not_doubling == "TRUE"), ]

# select unvaccinated subjects

ds <- ds[which(ds$Vaccinated_at_least_7_days_before_sampling_time == "no"), ]

# transform data into numeric data

for (i in 1:nrow(ds)) {
  
  element <- ds[i,"Neutralisation_May_2020"]
  temp<-unlist(strsplit(element, split = ":"))
  a<-temp[[2]]
  ds[i,"Neutralisation_May_2020"]<-a
}

for (i in 1:nrow(ds)) {
  
  element <- ds[i,"Neutralisation_June_2021"]
  temp<-unlist(strsplit(element, split = ":"))
  a<-temp[[2]]
  ds[i,"Neutralisation_June_2021"]<-a
}


ds$Neutralisation_May_2020 = as.numeric((as.character(ds$Neutralisation_May_2020)))
ds$Neutralisation_June_2021 = as.numeric((as.character(ds$Neutralisation_June_2021)))

#select subjects with neutraising activity in May (Neu>1:40)
ds <- ds[which(ds$Neutralisation_May_2020 > 40), ]

#define final dataset
Neu = data.frame(ds$Neutralisation_May_2020, ds$Neutralisation_June_2021)

# number of days between serosurveys
days_surveys <- as.numeric(as.Date("06/05/2021", "%m/%d/%Y") -
                             as.Date("05/01/2020", "%m/%d/%Y"))



# calculate halflife and CI
dt = Neu 
hl_db <- vector()
c=nrow(dt)
HL_N.sim = vector("numeric", length=4999)
Slopes=matrix(nrow = 4999, ncol = c)

for (idy in 1:c) 
  {
  for(idx in 1:4999)
  {   
    dt = data.frame(matrix(unlist(lapply(as.matrix(Neu), sample_FORW.fun)), ncol=ncol(Neu), nrow=nrow(Neu)))  
    HL_N.sim[idx] = median(days_surveys/log(dt[,1]/dt[,2],2))
    Slopes[idx, ]= log(dt[,2]/dt[,1])/days_surveys
  }
  
}
Final_slopes <- vector()
N=ncol(Slopes)
for (i in 1:N) {
  Final_slopes[i]=median(Slopes[,i])
}

Final_slopes
median(HL_N.sim)

rbind(quantile(HL_N.sim, c(0.025, 0.975)))

