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

dat4 <- subset(dat3, dat3$sp !="AMESAN" & dat3$sp !="AME sp." & dat3$sp !="ILEMUC" & dat3$sp !="MALPUM" & dat3$sp !="PINSTR"
               & dat3$sp !="POPBAL" & dat3$sp !="LARLAR" & dat3$sp !="SORAME" &
                 dat3$sp !="TSUCAN" & dat3$sp !="QUERUB" & dat3$sp !="CORALT" & dat3$sp !="AMELAE")


## Remove extreme values (outliers identified by Luke Magee)
## one distance = 0 
## Code with Alpha = 1.0 and Beta = 1.0: dat4 <- subset(dat4, dat4$Hdens != Inf & dat4$Cdens < 54 & dat4$Hdens < 2.5)  ## or huge outliers 
dat4 <- dat4[which(dat4$tag %in% c("119873", "113409", "094757", "075633", "075632", "027821") == F),]

## Remove smaller multi-stems from multi-stem individuals (keep largest stem
dat4 = dat4[order(dat4$tag, -dat4$dbh),]
dlist = split(dat4, as.character(dat4$tag))
for(i in 1:length(dlist)) {				# loop to keep largest alive stem for each individual
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

mod = lme4::glmer(surv ~ dbh + conspp + heterospp + Lconspp + Lheterospp + 
                    + (1|plotnum), family = binomial,data = dat4,
                  control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e6)))

mod <- glmer(surv ~ scale(dbh) + scale(Cdens) + scale(Hdens) + scale(LCdens) + scale(LHdens) + (1|quad),
             family = binomial, data = data_subset,
             control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))


library(lme4)
Cdens = dat4$Cdens
Hdens = dat4$Hdens
LCdens = dat4$LCdens
LHdens = dat4$LHdens

DenExp = 0.30	

# Scale all neighborhood densities to same scale (helpful for model interpretation later on)
# conhetero = c(Cdens^(DenExp), Hdens^(DenExp), LCdens^(DenExp), LHdens^(DenExp))	# Collective mean across all living and legacy con- and heterospecific densities trasnformed with D = DenExp
# dat4$conspp = ((Cdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled living conspecific densities
# dat4$heterospp = ((Hdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled living heterospecific densities
# dat4$Lconspp = ((LCdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled legacy conspecific densities
# dat4$Lheterospp = ((LHdens^(DenExp)) - mean(conhetero)) / sd(conhetero)	# Scaled legacy heterospecific densities
# dat4$dbh.2 <- scale(dat4$dbh)

dat4$conspp = 	Cdens^DenExp
dat4$heterospp = 	Hdens^DenExp
dat4$Lconspp = 		LCdens^DenExp
dat4$Lheterospp = LHdens^DenExp



# Create an empty data frame to store the coefficient estimates
coefficient_estimates <- data.frame()

dat4 <- subset(dat4, dat4$sp !="TSUCAN" & dat4$sp !="QUERUB" & dat4$sp !="CORALT" & dat4$sp !="AMELAE")
str(dat4$sp)
# Loop over each unique species in the dat4 data set
unique_species <- unique(dat4$sp)



for (sp1 in unique_species) {
  # Subset the data for the current species
  data_subset <- subset(dat4, dat4$sp == sp1)
  
  # Fit the model for the current species
  mod <- glmer(surv ~ scale(dbh) + scale(conspp) + scale(heterospp) + scale(Lconspp) + scale(Lheterospp) +  (1|quad),
               family = binomial, data = data_subset,
               control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))
  
  # Extract the coefficient estimates for fixed effects
  coefficients <- fixef(mod)
  
  # Set the species name as a column in the coefficients data frame
  coefficients <- data.frame(species = sp1, intercept = coefficients[1], dbh = coefficients[2],
                             Cdens = coefficients[3],Hdens = coefficients[4], LCdens = coefficients[5],
                             LHdens = coefficients[6])
  
  # Add the coefficients to the data frame
  coefficient_estimates <- rbind(coefficient_estimates, coefficients)
}

table(dat4$sp, dat4$surv)
rownames(coefficient_estimates) <- NULL
library(dplyr)

est.2 <- coefficient_estimates %>%
  select(species, everything()) %>%
  mutate(across(2:7, round, 3))

length(which(est.2$Cdens > est.2$LCdens))
## standard errors 
standard_errors <- data.frame()

for (sp1 in unique_species) {
  # Subset the data for the current species
  data_subset <- subset(dat4, dat4$sp == sp1)
  
  # Fit the model for the current species
  mod <- glmer(surv ~ scale(dbh) + scale(conspp) + scale(heterospp) + scale(Lconspp) + scale(Lheterospp) +  (1|quad),
               family = binomial, data = data_subset,
               control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))
  
  # Extract the standard errors for fixed effects
  se <- summary(mod)$coefficients[, "Std. Error"]
  
  # Set the species name as a column in the standard errors data frame
  se <- data.frame(species = sp1, intercept = se[1], dbh = se[2],
                   Cdens = se[3], Hdens = se[4], LCdens = se[5],
                   LHdens = se[6])
  
  # Add the standard errors to the data frame
  standard_errors <- rbind(standard_errors, se)
}
rownames(standard_errors) <- NULL
se <-  standard_errors %>%
  select(species, everything()) %>%
  mutate(across(2:7, round, 3))

# Identify the columns to modify (excluding the first column 'species')
columns_to_modify <- names(se)[-1]  # Excludes the first column 'species'

# Add '_se' to the column names
new_column_names <- paste0(columns_to_modify, "_se")

# Rename the columns in the dataframe
names(se)[match(columns_to_modify, names(se))] <- new_column_names

# Display the updated result
print(se)


# Assuming 'est.2' and 'se' are your dataframes

# Merge the dataframes on the 'species' column
merged_data <- merge(est.2, se, by = "species")

# Create a new dataframe with values in the format you described
result <- data.frame(
  species = merged_data$species,
  intercept = paste0(merged_data$intercept, " (", merged_data$intercept_se, ")"),
  dbh = paste0(merged_data$dbh, " (", merged_data$dbh_se, ")"),
  Cdens = paste0(merged_data$Cdens, " (", merged_data$Cdens_se, ")"),
  Hdens = paste0(merged_data$Hdens, " (", merged_data$Hdens_se, ")"),
  LCdens = paste0(merged_data$LCdens, " (", merged_data$LCdens_se, ")"),
  LHdens = paste0(merged_data$LHdens, " (", merged_data$LHdens_se, ")")
)

# Display the result
print(result)
# write.csv(result, "./manuscript/glmer_sp.csv", row.names = F)


# Initialize an empty dataframe to hold the p-values
p_values <- data.frame()

for (sp1 in unique_species) {
  # Subset the data for the current species
  data_subset <- subset(dat4, dat4$sp == sp1)
  
  # Fit the model for the current species
  mod <- glmer(surv ~ scale(dbh) + scale(conspp) + scale(heterospp) + scale(Lconspp) + scale(Lheterospp) +  (1|quad),
               family = binomial, data = data_subset,
               control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e6)))
  
  # Extract the p-values for fixed effects
  p_val <- summary(mod)$coefficients[, "Pr(>|z|)"]
  
  # Set the species name as a column in the p-values data frame
  p_val <- data.frame(species = sp1, intercept = p_val[1], dbh = p_val[2],
                      Cdens = p_val[3], Hdens = p_val[4], LCdens = p_val[5],
                      LHdens = p_val[6])
  
  # Add the p-values to the data frame
  p_values <- rbind(p_values, p_val)
}








