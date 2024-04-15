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
wab <- data.frame(readRDS("./data/wab_dat_growth_format_20220321.rds"))
##

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

# somethings wrong with the gx gy code below
# wab.focal = wab.focal[which(substr(wab.focal$quad,1,2) %in% c("00","14") == F & substr(wab.focal$quad,3,4) %in% c("00","41") == F),]
f <- data.table(wab.focal)

res.mat <- neighbours_2(dt = datt, f = f, r = 20, alpha = 0.6, beta = 1.1) ## linear functional form
# res.mat1 <- neighbours_2(dt = datt, f = f, r = 20, alpha = 1, beta = 1) ## linear functional form

## call it dat3
dat3 <- cbind(f, res.mat)
# dat3 <- cbind(f, res.mat1)   
## rename neighborhoods
colnames(dat3)[10:13] <- c('Cdens','Hdens','LCdens','LHdens')
    
## check mixed model 
dat3 <- dat3[dat3$status2 != "M",]
    
## survival = 1, mortality = 0
dat3$surv <- ifelse(dat3$status2 == "A", 1, 0)
    
    
dat4 <- subset(dat3, dat3$sp !="AMESAN" & dat3$sp !="AME sp." & dat3$sp !="ILEMUC" & dat3$sp !="MALPUM" & dat3$sp !="PINSTR"
                   & dat3$sp !="POPBAL" & dat3$sp !="LARLAR" & dat3$sp !="SORAME" &
                     dat3$sp !="TSUCAN" & dat3$sp !="QUERUB" & dat3$sp !="CORALT" & dat3$sp !="AMELAE")
    
## Remove extreme values (outliers identified by Luke Magee)
## distance = 0 for 
## trees mapped on other trees 
## Code with Alpha = 1.0 and Beta = 1.0: dat4 <- subset(dat4, dat4$Hdens != Inf & dat4$Cdens < 54 & dat4$Hdens < 2.5)  ## or huge outliers 
    dat4 <- dat4[which(dat4$tag %in% c("119873","111658", "113409", "094757", "075633", "075632", "027821") == F),]
    
    ## Remove smaller multi-stems from multi-stem individuals (keep largest stem
    dat4 = dat4[order(dat4$tag, -dat4$dbh),]
    dlist = split(dat4, as.character(dat4$tag))
    for(i in 1:length(dlist)) {	# loop to keep the largest alive stem for each individual
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
    
    dat4 = dat5
    
    

    
 
DenExp = 0.3
      
# load("./data/CTFSRPackage.rdata")

# dat4$plotID5 <- gxgy.to.quad(gx = dat4$x, gy = dat4$y, gridsize = 26,start = "zero", plotdim = c(300, 840))

## check moran's I 

# install.packages('spdep')
library(spdep)
library(sp)
library(spatialEco)
library(ggplot2)
# install.packages('adespatial')
library(adespatial)

# grid_sizes <- seq(5, 40, by = 5)
# 
# # Initialize an empty vector to store Moran's I results
# moran_results <- numeric(length(grid_sizes))
# moran_p <- numeric(length(grid_sizes))
# for (i in seq_along(grid_sizes)) {
  # Assign plotID based on the current grid size
# dat4$plotID <- gxgy.to.quad(gx = dat4$x, gy = dat4$y, gridsize = 30, start = "zero", plotdim = c(300, 840))
# 
# dat4$plotID5 <- gxgy.to.quad(gx = dat4$x, gy = dat4$y, gridsize = 5, start = "zero", plotdim = c(300, 840))

dat4 = dat4[which(substr(dat4$quad,1,2) %in% c("00","14") == F & substr(dat4$quad,3,4) %in% c("00","41") == F),]

# data_plt <- data.frame(
  Cdens = dat4$Cdens,
  Hdens = dat4$Hdens,
  LCdens = dat4$LCdens,
  LHdens = dat4$LHdens,
  # total = dat4$Cdens + dat4$Hdens,
  dbh = dat4$dbh
)



Cdens = dat4$Cdens
Hdens = dat4$Hdens
LCdens = dat4$LCdens
LHdens = dat4$LHdens

# Use GGally::ggpairs to create correlation plot
plt <- GGally::ggcorr(data_plt2, label = TRUE, label_size = 4, label_color = "black", lower = TRUE)

# Print the plot
print(plt)

# Rotate the text labels on y-axis
plt  + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1))

#####
conhetero = c(Cdens^(DenExp), Hdens^(DenExp), LCdens^(DenExp), LHdens^(DenExp)) # Collective mean across all living and legacy con- and heterospecific densities transformed with D = DenExp
dat4$conspp = ((Cdens^(DenExp)) - mean(conhetero)) / sd(conhetero) # Scaled living conspecific densities
dat4$heterospp = ((Hdens^(DenExp)) - mean(conhetero)) / sd(conhetero) # Scaled living heterospecific densities
dat4$Lconspp = ((LCdens^(DenExp)) - mean(conhetero)) / sd(conhetero) # Scaled legacy conspecific densities
dat4$Lheterospp = ((LHdens^(DenExp)) - mean(conhetero)) / sd(conhetero) # Scaled legacy heterospecific densities
dat4$dbh.2 <- scale(dat4$dbh)


data_plt2 <- data.frame(
  Cdens = dat4$conspp,
  Hdens = dat4$heterospp,
  LCdens = dat4$Lconspp,
  LHdens = dat4$Lheterospp,
  # total = dat4$Cdens + dat4$Hdens,
  dbh = dat4$dbh.2
)



colnames(dat4)
# length(unique(dat4$plotID5/plotID))

fm_free <- lme4::glmer(surv ~ dbh.2 + conspp + heterospp + Lconspp + Lheterospp + 
                  (conspp + heterospp + Lconspp + Lheterospp|sp) + (1|quad),
                   family = binomial, data = dat4, 
                  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))

fm_free1 <- lme4::glmer(surv ~ dbh.2 + conspp  + Lconspp + Lheterospp + 
                         (conspp + Lconspp + Lheterospp|sp) + (1|quad),
                       family = binomial, data = dat4, 
                       control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))

fm_check <- glm(surv ~ dbh + Cdens  + Hdens, family = "binomial", data = dat4)

summary(fm_check)

summary(fm_free)
library(sjPlot)

plot_model(fm_free, show.values = T, std.response = T )


# [1] 30476.25
library(car)
vif(fm_free)

d4 <- data.frame(dat4)
d4$res <- residuals(fm_free, type = "pearson")
# Calculate residuals
resi <- residuals(fm_free, type = "pearson")

# Create spatial weights matrix
W <- dnearneigh(d4[c("x", "y")], d1 = 1, d2 = 20)  ### 

W_listw <- nb2listw(W, style = "W")

# Perform Moran's I test
moran_test_result <- moran.mc(resi, W_listw, nsim = 199)

# Store Moran's I statistic in the vector


moran_results[i] <- moran_test_result$parameter
moran_p[i] <- moran_test_result$p.value
# }

####

# Assuming moran_results and moran_p are your vectors
library(ggplot2)

# Create a dataframe with the data
df <- data.frame(moran_results = moran_results,
                 moran_p = moran_p)

# Create a scatter plot with lines connecting successive points
ggplot(df, aes(x = grid_sizes , y =  moran_results)) +
  geom_point() +
  geom_line() +  # Add this line to connect successive points
  labs(x = "Grid size (m^2)",
       y = "Moran p-values",
       title = "")





# Assuming you have residuals in 'res' and coordinates in 'd4'
ggplot(data = d4, aes(x = x, y = y, fill = resi)) +
  geom_tile(width = 10, height = 10) +
  labs(title = "Tile Plot of Residuals vs. Coordinates",
       x = "X-coordinate",
       y = "Y-coordinate",
       fill = "Residuals") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal()



fm_nest <- lme4::glmer(surv ~ dbh.2 + conspp + heterospp + Lconspp + Lheterospp + 
                         (1|sp),
                       family = binomial, data = dat4, 
                       control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))

d4 <- data.frame(dat4)
d4$res <- residuals(fm_nest, type = "pearson")
# Calculate residuals
resi <- residuals(fm_nest, type = "pearson")

# Create spatial weights matrix
W <- dnearneigh(d4[c("x", "y")], d1 = .1, d2 = 20)  ### 

W_listw <- nb2listw(W, style = "W")

# Perform Moran's I test
moran_test_result <- moran.mc(resi, W_listw, nsim = 499)


