###############
rm(list=ls(all=TRUE))

library(dplyr)
library(data.table)
library(lme4)
library(piecewiseSEM)

# custom function opposite of %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

## data 
wab <- data.frame(readRDS("./wab_dat_growth_format_20220321.rds"))  ## need to set dir

## basal areas (m^2)
wab$ba.08 <- pi*(wab$dbh.08/2000)^2
wab$ba.13 <- pi*(wab$dbh.13/2000)^2
wab$ba.18 <- pi*(wab$dbh.18/2000)^2 

## trees alive in 2008
wab_sub_A <-subset(wab,status.08=="A")

wab.dat <- subset(wab_sub_A, 
                  select = c("plotX", "plotY","ba.08","spcode.18","status.13","tag.18", "stemID.18"))

colnames(wab.dat) <- c("x", "y", "ba", "sp", "status", "tag", "stemtag")

datt <- data.table(wab.dat) ## faster with data.table 

#### focal tree data set 
wab_sub_f <- subset(wab,status.13=="A" & !is.na(dbh.13))


wab_sub_dat <- wab_sub_f

## select columns 
wab.focal <-subset(wab_sub_dat, 
                   select = c("plotX", "plotY","dbh.13","spcode.18","status.13","status.18","tag.18","quadrat_txt.18", "stemID.18"))

# less typing way to do it
wab.focal <- wab_sub_dat[, c("plotX", "plotY","dbh.13","spcode.18",
                             "status.13","status.18",
                             "tag.18","quadrat_txt.18", "stemID.18")]

## rename for idw basal area formula 
colnames(wab.focal) <- c("x", "y", "dbh", "sp", "status","status2", "tag","quad", "stemtag")

## set focal stems (exclude columns 00 & 14 as well as rows 00 & 41 to remove edge effects)

# somethings wrong with the gx gy code below
# wab.focal = wab.focal[which(substr(wab.focal$quad,1,2) %in% c("00","14") == F & substr(wab.focal$quad,3,4) %in% c("00","41") == F),]
f <- data.table(wab.focal)


big_dat <- readRDS("./Wab_long_neighbs_grouped_no_exp.RDS")  ## need to set dir

d2 <- big_dat

d2$dbh1 <- sqrt(4 * d2$ba / pi) * 1000

# table(d2$dist < (d2$dbh/2+d2$dbh1/2)/1000)

d2$dist <- ifelse(d2$dist < (d2$dbh/2+d2$dbh1/2)/1000,
                  (d2$dbh/2+d2$dbh1/2)/1000, d2$dist)

big_dat <- d2

#### now read in exponents 
params <- read.csv("./rev3_search1.csv")  ## fix dir 


###
s <- Sys.getenv("SLURM_ARRAY_TASK_ID")  ### to set for numbers of jobs
params_try = params[s, ]

Cdens_a = params_try$Cdens_a
Cdens_b = params_try$Cdens_b
Cdens_d = params_try$Cdens_d

Hdens_a = params_try$Hdens_a
Hdens_b = params_try$Hdens_b
Hdens_d = params_try$Hdens_d

LCdens_a = params_try$LCdens_a
LCdens_b = params_try$LCdens_b
LCdens_d = params_try$LCdens_d

LHdens_a = params_try$LHdens_a
LHdens_b = params_try$LHdens_b
LHdens_d = params_try$LHdens_d





neighborhood_params <- data.table(
  type = c("ca", "ha","cd", "hd"),
  alpha = c(Cdens_a, Hdens_a, LCdens_a, LHdens_a),
  beta = c(Cdens_b, Hdens_b, LCdens_b, LHdens_b)
)

# Join the neighborhood parameters with your data.table
dt <- merge(big_dat, neighborhood_params, by = "type", all.x = TRUE)

# Calculate dens using data.table separated by status, species, stemtag, type, focal, and neighborhood type
Dt2 <- dt[, dens := sum(ba^alpha / dist^beta), by = .(type, focal)]

result <- dcast(setDT(Dt2), focal ~ type, value.var = "dens", fun.aggregate =  mean)

# Rename the columns as per your requirement
setnames(result, c("focal", "ca", "ha", "cd", "hd"),  c("stemtag", "Cdens", "Hdens", "LCdens", "LHdens"))


## call it dat3
dat3 <- merge(f, result, by = "stemtag", all.x = T)

## remove missing 
dat3 <- dat3[dat3$status2 != "M",]

## survival = 1, mortality = 0
dat3$surv <- ifelse(dat3$status2 == "A", 1, 0)

## Additional data filtering and calculations
dat4 <- subset(dat3, dat3$sp !="AMESAN" & dat3$sp !="AME sp." & dat3$sp !="ILEMUC" & dat3$sp !="MALPUM" & dat3$sp !="PINSTR"
               & dat3$sp !="POPBAL" & dat3$sp !="LARLAR" & dat3$sp !="SORAME" &
                 dat3$sp !="TSUCAN" & dat3$sp !="QUERUB" & dat3$sp !="CORALT" & dat3$sp !="AMELAE")


## Remove extreme values trees mapped on other trees ##
dat4$tag <- substr(dat4$stemtag, 1, 6)
dat4 <- dat4[which(dat4$tag %in% c("119873", "111658", "113409", "094757", "075633", "075632", "027821") == F),]

## Remove smaller multi-stems from multi-stem individuals
dat4 = dat4[order(dat4$tag, -dat4$dbh),]
dlist = split(dat4, as.character(dat4$tag))
for(i in 1:length(dlist)) {
  test = dlist[[i]]
  if(any(test$status2 == "D") & any(test$status2 == "A")) {
    test = test[which(test$status2 == "A"),]
    test = test[which(test$dbh == max(test$dbh))[1],]
  } else {
    test = test[which(test$dbh == max(test$dbh))[1],]
  }
  dlist[[i]] = test
}
dat5 = do.call('rbind', dlist)
dat5 = dat5[which(substr(dat5$quad,1,2) %in% c("00","14") == F & substr(dat5$quad,3,4) %in% c("00","41") == F),]

dat4 = dat5

####
dat4$Cdens <- ifelse(is.na(dat4$Cdens), 0, dat4$Cdens)
dat4$Hdens <- ifelse(is.na(dat4$Hdens), 0, dat4$Hdens)
dat4$LCdens <- ifelse(is.na(dat4$LCdens), 0 , dat4$LCdens)
dat4$LHdens <- ifelse(is.na(dat4$LHdens), 0, dat4$LHdens)

## variables   
Cdens = dat4$Cdens
Hdens = dat4$Hdens
LCdens = dat4$LCdens
LHdens = dat4$LHdens


DenExp = Cdens_d
DenExp2 = Hdens_d
DenExp3 = LCdens_d
DenExp4 = LHdens_d


## scale across data 
conhetero = c(Cdens^(DenExp), Hdens^(DenExp2), LCdens^(DenExp3), LHdens^(DenExp4))
dat4$conspp = ((Cdens^(DenExp)) - mean(conhetero)) / sd(conhetero)
dat4$heterospp = ((Hdens^(DenExp2)) - mean(conhetero)) / sd(conhetero)
dat4$Lconspp = ((LCdens^(DenExp3)) - mean(conhetero)) / sd(conhetero)
dat4$Lheterospp = ((LHdens^(DenExp4)) - mean(conhetero)) / sd(conhetero)
dat4$dbh.2 <- scale(dat4$dbh)


## GLMM model with random slopes for species 
fm3 <- lme4::glmer(surv ~ dbh.2 + conspp + heterospp + Lconspp + Lheterospp + 
                     (conspp + heterospp + Lconspp + Lheterospp|sp) + (1|quad),
                   family = binomial, data = dat4,
                   control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))

# Save the output for this trial 
output = data.frame(s = s, 
                    Cdens_a = Cdens_a, Cdens_b = Cdens_b, Cdens_d = Cdens_d,
                    Hdens_a = Hdens_a, Hdens_b = Hdens_b, Hdens_d = Hdens_d,
                    LCdens_a = LCdens_a, LCdens_b = LCdens_b, LCdens_d = LCdens_d,
                    LHdens_a = LHdens_a, LHdens_b = LHdens_b, LHdens_d = LHdens_d,
                    logLik = logLik(fm3), AIC = AIC(fm3))

## create file name 
filename  <- paste("./revisions3/output_", s, ".csv", sep = "")

## save to directory it has the file name of the trial 
write.csv(output, filename)
