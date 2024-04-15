## code to setup neighborhoods and models 
## 
rm(list=ls(all=TRUE))
# install.package
library(dplyr)
library(data.table)

## data 
###############
rm(list=ls(all=TRUE))

library(dplyr)
library(data.table)
library(lme4)
library(piecewiseSEM)

# custom function opposite of %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

## read in neighborhood function 
source("./source/neighborhood_functional_forms_coefficeints_ljm.R")

## data 
## full data are availale upon request --we provide a subset
# wab <- data.frame(readRDS("./data/wab_dat_growth_format_20220321.rds"))

# subset 
wab <- readRDS("./Data/wabikon_subset.RDS")


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

### take a random sample to make MLEs quicker n = 1000
# set.seed(14)
# wab_sub_dat <- sample_n(wab_sub_f, 10000)
# wab_sub_dat <- subset(wab_sub_f, wab_sub_f$dbh.18 < 100) ## small trees 
wab_sub_dat <- wab_sub_f


summary(wab_sub_dat)

## select columns 
wab.focal <- subset(wab_sub_dat, 
                    select = c("plotX", "plotY","dbh.13","spcode.18","status.13","status.18","tag.18","quadrat_txt.18", "stemID.18"))

# less typing way to do it
wab.focal <- wab_sub_dat[, c("plotX", "plotY","dbh.13","spcode.18",
                             "status.13","status.18",
                             "tag.18","quadrat_txt.18", "stemID.18")]

## rename for idw basal area formula 
colnames(wab.focal) <- c("x", "y", "dbh", "sp", "status","status2", "tag","quad", "stemtag")

## set focal stems (exclude columns 00 & 14 as well as rows 00 & 41 to remove edge effects)
wab.focal = wab.focal[which(substr(wab.focal$quad,1,2) %in% c("00","14") == F & substr(wab.focal$quad,3,4) %in% c("00","41") == F),]
f <- data.table(wab.focal)

## apply function, set alpha and beta to 1
res.mat <- neighbours_2(dt=datt,f=f,r=20,alpha=0.6,beta=1.1) ## linear functional form

## call it dat3
dat3 <- cbind(f, res.mat)

## rename neighborhoods
colnames(dat3)[10:13] <- c('Cdens','Hdens','LCdens','LHdens')

## check mixed model 
dat3 <- dat3[dat3$status2!="M",]

## survival = 1, mortality = 0
dat3$surv <- ifelse(dat3$status2=="A",1,0)

saveRDS(dat3, "./Data/subset_wabikon_neigbhors.RDS")
