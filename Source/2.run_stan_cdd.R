
# Load necessary R packages
library(vegan)
library(doBy)
library(arm)
library(reshape)
library(nlme)
library(spatstat)
library(lme4)
library(vegan)
library(boot)
library(lmtest)
library(abind)
library(optimx)
library(nloptr)
library(parallel)
library(MASS)
library(stats)
library(DescTools)
library(purrr)
library(rstan)
library(rstanarm)
# Code for multi-core
options(mc.cores = parallel::detectCores())
library(DHARMa)


## read in data (from data prep code)  --again we provide a subset here from the preparation script
dat3 <- readRDS("./Data/subset_wabikon_neigbhors.RDS")
# dat3 <- readRDS("./data/Data.legacy_effects_alpha0-6_beta1-1.rds")

## Remove species with not enough records 
dat4 <- subset(dat3, dat3$sp !="AMESAN" & dat3$sp !="AME sp." & dat3$sp !="ILEMUC" & dat3$sp !="MALPUM" & dat3$sp !="PINSTR"
               & dat3$sp !="POPBAL" & dat3$sp !="LARLAR" & dat3$sp !="SORAME" &
                 dat3$sp !="TSUCAN" & dat3$sp !="QUERUB" & dat3$sp !="CORALT" & dat3$sp !="AMELAE")


## some requirements 
dat4$sp <- as.factor(dat4$sp)
dat4$quad <- as.factor(dat4$quad)

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


### IN ADDITION MAKE SURE THESE RData FILES ARE IN R WORKING DIRECTORY:
# "ForestGEO_LogisticSurvival_Model_STAN.stan"


# Load other analysis functions

# Function to standardize conspecific and heterospecific densities on the same scale
std.conhet.density = function(conspp, heterospp, id, DenExp) {
  consppid = tapply(conspp, id, mean)							# Mean conspecific density across individuals
  heterosppid = tapply(heterospp, id, mean)						# Mean heterospecific density across individuals
  conhetero = c(consppid^(DenExp), heterosppid^(DenExp))				# Collective mean across all con- and heterospecific densities trasnformed with D = DenExp
  conspp.std = ((conspp^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled conspecific densities
  heterospp.std = ((heterospp^(DenExp)) - mean(conhetero)) / sd(conhetero)	# Scaled heterospecific densities
  return(list(conspp = conspp.std, heterospp = heterospp.std))
}


# Function to perform Bayesian bootstrap (Rubin 1981, Gustafson 2007)
# see https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution for proof of gamma distribution equivalency to Dirichlet distribution
# sample.weights is an optional agrument used to weight Bayesian bootstrap average predictive comparisons by plot area
bayesian.bootstrap = function(test_sample, nIts, sample.weights = NA) {
  weights <- matrix(rexp(length(test_sample) * nIts, 1), ncol = length(test_sample), byrow = TRUE)	# Matrix to hold weights sampled from Dirichlet distribution with equal weights
  if(is.numeric(sample.weights)) {weights <- t(t(weights) * sample.weights)}	# Weight by sample.weights (e.g. forest-plot area in ha)
  weights <- weights / rowSums(weights)							# Make sure weights sum to 1 (Dirichlet distribution)
  result = rowSums(t(t(weights) * test_sample))						# Multiply and add test_sample to weights for each iteration
  return(result)
}




#######################################################################################
### BAYESIAN MODEL OF TREE SURVIVAL	

DenExp = 0.30			# Density exponent (non-linearity parameter in neighborhood function; selected from GLMER analysis loop)
d = dat5

# Prepare data for analysis 
surv = d$surv
hazard = 1 - surv
uniID = 1:nrow(d)
Cdens = d$Cdens
Hdens = d$Hdens
LCdens = d$LCdens
LHdens = d$LHdens

# Scale all neighborhood densities to same scale (helpful for model interpretation later on)
conhetero = c(Cdens^(DenExp), Hdens^(DenExp), LCdens^(DenExp), LHdens^(DenExp))	# Collective mean across all living and legacy con- and heterospecific densities trasnformed with D = DenExp
conspp = ((Cdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled living conspecific densities
heterospp = ((Hdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled living heterospecific densities
Lconspp = ((LCdens^(DenExp)) - mean(conhetero)) / sd(conhetero)		# Scaled legacy conspecific densities
Lheterospp = ((LHdens^(DenExp)) - mean(conhetero)) / sd(conhetero)	# Scaled legacy heterospecific densities

allspp = c(scale((Cdens + Hdens)^(DenExp)))
dbh.raw = d$dbh/10
dbh = (dbh.raw - mean(dbh.raw)) / sd(dbh.raw)					# Standardized DBH
sp = as.numeric(factor(as.character(d$sp)))					# Species numbers (need numbers for STAN)
sp.info = unique(data.frame(sp.name = d$sp, sp = sp))				# Species info (relating species number to species name)
n = nrow(d)											# Total samples
plotnum = as.numeric(factor(as.character(factor(as.character(d$quad)))))# Quadrats indicator for quadrat random effect
n.plot = length(unique(plotnum))
plot.info = unique(data.frame(plot = d$quad, plotnum = plotnum))

# Predictor matrices for model (first column = intercept, subsequent columns = covariates)
x = cbind(rep(1,times = n), dbh, conspp, heterospp, Lconspp, Lheterospp)# Fixed effects model matrix
K = ncol(x)											# Number of fixed effects parameters
SX = cbind(rep(1,times = n), conspp, heterospp, Lconspp, Lheterospp)	# Species predictor matrix (1st column = random intercept for each sp; subsequent columns = random slopes)
SK = ncol(SX)										# Num of species predictors
SJ = length(unique(sp))									# Num of species


# Bundle data for STAN
dat <- list(y = surv, x = x, K = K, quadnum = plotnum, nquad = n.plot, sp = sp, SJ = SJ, SK = SK, SX = SX, N = n)

# GLMER model with ML for ranef variance estimates (to obtain initial values for STAN HMC model below)
set.seed(314)
m1 = lme4::glmer(surv ~ dbh + conspp + heterospp + Lconspp + Lheterospp + 
                   (conspp + heterospp + Lconspp + Lheterospp|sp) + (1|plotnum), family = binomial,
                 control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e6)))
summary(m1)
logLik(m1)



### Fit STAN model 
# HMC Parameters
nchains = 5					# Number of chains to run
postburn = 4000				# Number of post-burn-in posterior samples
burnin = 500				# Number of burn-in samples
its = postburn + burnin			# Total number of samples to be drawn from posterior
thin = 10					# Rate at which to thin the posterior samples		
hmc_seed = 55406				# HMC seed for reproducibility

# Initial values from LMER model  -- not required but can be helpful 
# vc <- VarCorr(m1)
# sigma_SB_lmer = as.matrix(Matrix::bdiag(vc$sp))
# inits <- replicate(nchains, list(
#   Sz = t(data.matrix(ranef(m1)$'sp')),
#   SL_Omega = sigma_SB_lmer,
#   Stau_unif = runif(SK),
#   beta = fixef(m1),
#   sigma_QUAD = attr(summary(m1)$varcor$plotnum,"stddev")[[1]],
#   QUAD = ranef(m1)$'plotnum'[[1]], 
#   SB = data.matrix(ranef(m1)$'sp'), 
#   sigma_SB = sigma_SB_lmer, 
#   Stau = runif(SK)
# ), simplify = F)


# Run model
begin = Sys.time()
fit <- stan(file = "./source/ForestGEO_LogisticSurvival_Model_STAN.stan", 		# Stan model
            data = dat,   										# named list of data
            # init = inits,										# initial values for parameters
            chains = nchains,             							# number of Markov chains
            warmup = burnin,          								# number of warmup iterations per chain
            iter = its,            								# total number of iterations per chain
            cores = nchains,              							# number of cores (could use one per chain)
            thin = thin,										# Thin rate
            seed = hmc_seed,									# HMC seed for reproducibility
            control = list(adapt_delta = 0.99)						# HMC sampling parameter
)
end = Sys.time()
(duration = end - begin)

#  save(fit, file = "./Data/stan_output_test.RData")
# 
# load("./Data/stan_output_test.RData")

median(fit@sim$samples)
summary(fit)
getwd()

# HMC traceplots
traceplot(fit, pars = c("beta"), inc_warmup = TRUE, nrow = 3)
traceplot(fit, pars = c("sigma_QUAD"), inc_warmup = TRUE, nrow = 3)
traceplot(fit, pars = c("sigma_SB"), inc_warmup = TRUE, nrow = 3)
traceplot(fit, pars = c("SL_Omega"), inc_warmup = TRUE, nrow = 3)
traceplot(fit, pars = c("Sz"), inc_warmup = TRUE, nrow = 3)
traceplot(fit, pars = c("Stau"), inc_warmup = TRUE, nrow = 3)



# test for auto-correlation in chains
draws <- extract(fit, permuted = FALSE)
acf(c(draws[,,"beta[1]"]), lag.max = 100)
acf(c(draws[,,"beta[2]"]), lag.max = 100)
acf(c(draws[,,"beta[3]"]), lag.max = 100)
acf(c(draws[,,"beta[4]"]), lag.max = 100)
acf(c(draws[,,"beta[5]"]), lag.max = 100)
acf(c(draws[,,"beta[6]"]), lag.max = 100)
acf(c(draws[,,"sigma_QUAD"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[1,1]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[2,2]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[3,3]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[4,4]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[5,5]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[1,2]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[1,3]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[1,4]"]), lag.max = 100)
acf(c(draws[,,"sigma_SB[1,5]"]), lag.max = 100)






# Extract estimates
draws <- extract(fit, permuted = FALSE)
beta.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,4) == "beta")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,4) == "beta"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,4) == "beta")])),3,c)
SB.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,2) == "SB")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,2) == "SB"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,2) == "SB")])),3,c)
SB.m = array(SB.m, dim = c(nrow(SB.m),SJ,SK), dimnames = list(paste("it",c(1:nrow(SB.m)),sep=""),paste("sp",c(1:SJ),sep=""),paste("SB",c(1:SK),sep="")))
quad.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,4) == "QUAD")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,4) == "QUAD"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,4) == "QUAD")])),3,c)
sigmaSB.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,8) == "sigma_SB")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,8) == "sigma_SB"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,8) == "sigma_SB")])),3,c)
sigmaQUAD.m = c(draws[,,"sigma_QUAD"])

ll.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,7) == "log_lik")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,7) == "log_lik"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,7) == "log_lik")])),3,c)
fittedvalues = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,6) == "fitted")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,6) == "fitted"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,6) == "fitted")])),3,c)

Sz.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,2) == "Sz")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,2) == "Sz"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,2) == "Sz")])),3,c)
Stau.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,4) == "Stau")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,4) == "Stau"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,4) == "Stau")])),3,c)
SL_Omega.m = apply(array(draws[,,which(substr(attributes(draws)$dimnames$parameters,1,8) == "SL_Omega")], dim = c((dim(draws)[1]), dim(draws)[2], length(which(substr(attributes(draws)$dimnames$parameters,1,8) == "SL_Omega"))),dimnames = list(c(),c(),attributes(draws)$dimnames$parameters[which(substr(attributes(draws)$dimnames$parameters,1,8) == "SL_Omega")])),3,c)
(SL_Cholesky = diag(apply(Stau.m,2,mean)[6:10])%*%matrix(apply(SL_Omega.m,2,mean),nrow=5,ncol=5))
matrix(apply(sigmaSB.m,2,mean),nrow=5,ncol=5)
cov2cor(matrix(apply(sigmaSB.m,2,mean),nrow=5,ncol=5))

getME(m1,"theta")	# Cholesky factors from the lmer model

# Estimates
t(apply(beta.m,2,quantile,c(0.025,0.975,0.05,0.95,0.5))); apply(beta.m,2,mean)
t(apply(sigmaSB.m,2,quantile,c(0.025,0.975,0.05,0.95,0.5))); apply(sigmaSB.m,2,mean)
quantile(sigmaQUAD.m^2,c(0.025,0.975,0.05,0.95,0.5)); mean(sigmaQUAD.m^2)




# write.csv(beta.m, "./manuscript/betas.csv", row.names = F)






######################################
### GOODNESS OF FIT TESTS

# Residual plots using 'DHARMa' R package

# Prepare data and model
log_lik = ll.m
fitted.m = inv.logit(fittedvalues)
dev.m = rowSums(log_lik)*-2
Z = length(sigmaQUAD.m)

set.seed(314)
rep2.m = list()
for(i in 1:nrow(fittedvalues)){rep2.m[[i]] = as.numeric(rbinom(n = length(inv.logit(fittedvalues[i,])), size = 1, prob = inv.logit(fittedvalues[i,])))}	# simulate survival or death for each observation at each posterio sample for 'DHARMa' package
rep2.m = do.call('rbind', rep2.m)


# DHARMa scaled residuals
sim2 = createDHARMa(simulatedResponse = t(rep2.m), observedResponse = as.vector(surv), fittedPredictedResponse = apply(inv.logit(fittedvalues), 2, median), integerResponse = T)
plot(sim2)
# testOutliers(sim2, type = "bootstrap")
plotResiduals(sim2, xlab = "Scaled predicted survival probability")
plotResiduals(sim2, form = conspp, xlab = "Scaled living conspecific density")
plotResiduals(sim2, form = heterospp, xlab = "Scaled living heterospecific density")
plotResiduals(sim2, form = Lconspp, xlab = "Scaled legacy conspecific density")
plotResiduals(sim2, form = Lheterospp, xlab = "Scaled legacy heterospecific density")
plotResiduals(sim2, form = dbh, xlab = "Scaled minimum DBH")
plotResiduals(sim2, form = factor(sp), xlab = "Species ID")
# plotResiduals(sim2, form = factor(plotnum), xlab = "Quadrat ID")



#pdf("HJA_5-15cmDBH_Survival_Bayesian_residual_plots_main_20210914.pdf", height = 5, width = 9, useDingbats=FALSE)
plot(sim2)
#dev.off()

#pdf("HJA_5-15cmDBH_Survival_Bayesian_residual_plots_predictors_20210914.pdf", height = 9, width = 8, useDingbats=FALSE)
set.seed(314)
par(mfrow=c(3,2))
par(mar=c(4,4,2,2))
plotResiduals(sim2, form = conspp, xlab = "Scaled living conspecific density")
plotResiduals(sim2, form = heterospp, xlab = "Scaled living heterospecific density")
plotResiduals(sim2, form = Lconspp, xlab = "Scaled legacy conspecific density")
plotResiduals(sim2, form = Lheterospp, xlab = "Scaled legacy heterospecific density")
plotResiduals(sim2, form = dbh, xlab = "Scaled minimum DBH")
plotResiduals(sim2, form = factor(sp), xlab = "Species ID")
#testDispersion(sim2)
#dev.off()











# Leave-one-out fit metrics (LOO)
# *Run Goodness-of-fit R code above first, the LOO code below depends on object 'llid' which is the log-likelihood matrix for each individual & each posterior sample

# Leave-one-out fit metrics
library(loo)
begin = Sys.time()
r_eff <- relative_eff(exp(ll.m), chain_id = rep(1:5,each=400), cores = 5)
loo_1 <- loo(ll.m, r_eff = r_eff, cores = 5)
end = Sys.time()
(duration = end - begin)
print(loo_1)








# logLik and dev at mean posterior parameters (for DIC calculation)
beta.mean = apply(beta.m, 2, mean)
SB.mean = apply(SB.m,c(2,3),mean)
quad.mean = apply(quad.m, 2, mean)
SBbysp = matrix(NA,nrow=n,ncol=SK)					# Prepare model matrices for predictions
for(i in 1:n) {SBbysp[i,] = SB.mean[sp[i],]}
QUADbycensus = c()
for(i in 1:n) {QUADbycensus[i] = quad.mean[plotnum[i]]}
fitted.meanpar = inv.logit((x %*% beta.mean) + rowSums(SX*SBbysp) + QUADbycensus)	# Fitted instantaneous hazard from model at mean posterior values
logLik.meanpar = (log(fitted.meanpar) * surv) + (log(1-fitted.meanpar) * (1-surv))	# Log-likelihood for each individual at each interval at mean posterior values
logLik.meanpar.sum = sum(logLik.meanpar)								# Summed log-likelihood at mean posterior values
dev.meanpar = -2 * logLik.meanpar.sum								# Model deviance at mean posterior values
pD = mean(dev.m) - dev.meanpar									# pD for model
DIC = dev.meanpar + (2*pD)										# DIC for model
dev.meanpar; pD; DIC









## For Posterior Predictive Check (based on all samples from the posterior parameter distributions)
## Requires object 'predS.phi.m' from residual plot code above
rep.m = list()
sim_log_lik = list()
prop.predict = c()
prop.predict.death = c()
prop.predict.live = c()

set.seed(314)
for(i in 1:nrow(fitted.m)){
  rep.m[[i]] = as.numeric(rbernoulli(n = length(fitted.m[i,]), p = (fitted.m[i,])))			# Simulated survival for each observation
  sim_log_lik[[i]] = (log(fitted.m[i,]) * rep.m[[i]]) + (log(1 - fitted.m[i,]) * (1 - rep.m[[i]]))		# Simulated log-likelihood for each obs (logistic regression likelihood)
}
rep.m = do.call('rbind', rep.m)
sim_log_lik = do.call('rbind', sim_log_lik)
for(i in 1:nrow(rep.m)) {
  prop.predict[i] = mean(surv == rep.m[i,])				# Proportion of individuals where survival/mortlaity is correctly simulated based on model at each posterior sample
}
for(i in 1:nrow(rep.m)) {
  prop.predict.death[i] = 1-mean(rep.m[i,][which(surv == 0)])		# Proportion of dead individuals where mortlaity is correctly simulated based on model at each posterior sample
}
for(i in 1:nrow(rep.m)) {
  prop.predict.live[i] = mean(rep.m[i,][which(surv == 1)])	# Proportion of surviving individuals where survival is correctly simulated based on model at each posterior sample
}

summary(prop.predict)						# posterior classification accuracy
summary(prop.predict.live)					# posterior classification accuracy of survival
summary(prop.predict.death)					# posterior classification accuracy of death
summary((prop.predict.death + prop.predict.live)/2)	# posterior balanced classification accuracy


# LogLik posterior predictive test
loglik.obs = rowSums(log_lik)
loglik.sim = rowSums(sim_log_lik)
summary(loglik.obs)
summary(loglik.sim)
#pdf("HJA_5to15cmDBH_ExpSurvival_logLik_PPC_BAYES_20220129.pdf", height = 4, width = 5, useDingbats=FALSE)
hist(loglik.sim, breaks=30, xlab = "Log-likelihood of simulated data", las = 1, main = "5-15 cm DBH survival log-likelihood PPC")
abline(v = mean(loglik.obs), lwd = 2, col = "red")
sum(loglik.sim > mean(loglik.obs)) / nrow(rep.m)
#dev.off()


# Mortality posterior predictive test
mean(1 - surv)
mean.hazard.sim = apply(1 - rep.m,1,mean)
summary(mean.hazard.sim)
#pdf("HJA_5to15cmDBH_ExpSurvival_Mortality_PPC_BAYES_20220129.pdf", height = 4, width = 5, useDingbats=FALSE)
hist(mean.hazard.sim, breaks=20, xlab = "Proportion of individuals dead in simulated data", las = 1, main = "5-15 cm DBH survival mortality PPC")
abline(v = mean(1-surv), lwd = 2, col = "red")
sum(mean.hazard.sim > mean(1-surv)) / nrow(rep.m)
#dev.off()


# Posterior predictive tests for mortality levels in each species, each plot, and each species-by-plot combo
sp.test = list()
for(i in 1:nrow(rep.m)) {
  sp.test[[i]] = tapply(1 - rep.m[i,], sp, mean)
}
sp.test = do.call('rbind', sp.test)
obs.sp.test = tapply(1-surv, sp, mean)
sp.test.p = c()
for(j in 1:ncol(sp.test)) {
  sp.test.p[j] = min(sum(sp.test[,j] >= obs.sp.test[j]), sum(sp.test[,j] <= obs.sp.test[j])) / length(sp.test[,j])
  if(obs.sp.test[j] == 0 | obs.sp.test[j] == 1) {sp.test.p[j] = (sum(sp.test[,j] == obs.sp.test[j]) / length(sp.test[,j]))/2}
}
summary(sp.test.p); sum(sp.test.p < 0.025 | sp.test.p > 0.975)
tapply(1-surv, sp, mean)
apply(sp.test,2,mean)
round(tapply(1-surv, sp, mean) - apply(sp.test,2,mean), 4)











###################################################################################
# CDD-HDD average predictive comparisons and visualize model

DenExp = 0.3			# D parameter (nonlinearity parameter) selected based on data
StdAdultDen = 1.225848		# Mean total living tree density across all tree neighborhoods 
StdLAdultDen = 0.08160748	# Mean total legacy tree density across all tree neighborhoods 

# Values for calculation of average marginal effects (adding one standardized conspecific tree at a standardized distance to compare CDD and legacy-CDD across species)
StdTreeSize = 3.0			# Standardized conspecific tree size (minimum of maximum tree diameters (in cm DBH) across species)
StdTreeSizeBA = pi*((StdTreeSize/200)^2)	# Standardized conspecific tree size in BA (StdTreeSize converted to basal area in m^2)
StdTreeDist = 1.0			# Standardized conspecific tree distance from focal tree (in meters)
Alpha = 0.6				# Alpha parameter from neighborhood density eq. (size-dependence of density effect, higher values give greater weight to larger trees)
Beta = 1.1				# Beta parameter from neighborhood density eq. (distance-dependence of density effect, higher values give greater weight to closer trees)			
oneStdtree = (StdTreeSizeBA^Alpha)/(StdTreeDist^Beta)	# One standardized tree at standardized distance calculated with the neighborhood density eq.

ZeroTree = (0 - mean(conhetero)) / sd(conhetero)					# Standardized density index with zero trees (only used to calculate separate conspecific and heterospecific effects, not for CDD-HDD)
OneConTree = ((oneStdtree^DenExp) - mean(conhetero)) / sd(conhetero)		# Standardized density index at mean tree density (for comparing conspecific stand to heterospecific stand)
AllHetTree = ((StdAdultDen^DenExp) - mean(conhetero)) / sd(conhetero)		# Standardized density index at half of mean tree density (for determining survival in a stand composed of half conspecifics and half heterospecifics)
AllLHetTree = ((StdLAdultDen^DenExp) - mean(conhetero)) / sd(conhetero)		# Standardized density index at half of mean tree density (for determining survival in a stand composed of half conspecifics and half heterospecifics)
AllHetTreeMinusOneConTree = (((StdAdultDen - oneStdtree)^DenExp) - mean(conhetero)) / sd(conhetero)		# Standardized density index at half of mean tree density (for determining survival in a stand composed of half conspecifics and half heterospecifics)
AllLHetTreeMinusOneConTree = (((StdLAdultDen - oneStdtree)^DenExp) - mean(conhetero)) / sd(conhetero)		# Standardized density index at half of mean tree density (for determining survival in a stand composed of half conspecifics and half heterospecifics)

# Species table to hold estimates for each species
sp.info2 = sp.info[order(sp.info$sp),]

# Data for avg. perdictive comparisons (Hi = "high" living conspecific density; defined as mean density of trees, but some living trees are conspecific, standardized to oneStdtree set above)
xHi = cbind(rep(1,times = n), dbh, OneConTree, AllHetTreeMinusOneConTree, Lconspp, Lheterospp)	
SXHi = cbind(rep(1,times = n), OneConTree, AllHetTreeMinusOneConTree, Lconspp, Lheterospp)		# Species predictor matrix

# Data for avg. perdictive comparisons (Lo = "low" living conspecific density; defined as mean density of trees, all living trees are heterospecifics)
xLo = cbind(rep(1,times = n), dbh, ZeroTree, AllHetTree, Lconspp, Lheterospp)	
SXLo = cbind(rep(1,times = n), ZeroTree, AllHetTree, Lconspp, Lheterospp)		# Species predictor matrix

# Data for avg. perdictive comparisons for legacy CDD (Hi = "high" legacy conspecific density; defined as mean density of trees, but some legacy trees are conspecific, standardized to oneStdtree set above)
xLHi = cbind(rep(1,times = n), dbh, conspp, heterospp, OneConTree, AllLHetTreeMinusOneConTree)	
SXLHi = cbind(rep(1,times = n), conspp, heterospp, OneConTree, AllLHetTreeMinusOneConTree)		# Species predictor matrix

# Data for avg. perdictive comparisons for legacy CDD (Lo = "low" legacy conspecific density; defined as mean density of trees, all legacy trees are heterospecifics)
xLLo = cbind(rep(1,times = n), dbh, conspp, heterospp, ZeroTree, AllLHetTree)	
SXLLo = cbind(rep(1,times = n), conspp, heterospp, ZeroTree, AllLHetTree)		# Species predictor matrix


# Calculate average predictive comparisons for each species (Gelman & Pardoe 2007, Gustafson 2007)

realSurvLivingConHi.sp = list()
realSurvLivingConLo.sp = list()
realSurvLegacyConHi.sp = list()
realSurvLegacyConLo.sp = list()
realLivingCDD.sp = list()
realLegacyCDD.sp = list()

begin.time = Sys.time()
for(j in 1:SJ) {
  Z = length(sigmaQUAD.m)
  realSurvLivingConHi.m = matrix(NA,nrow=Z,ncol=n)
  realSurvLivingConLo.m = matrix(NA,nrow=Z,ncol=n)
  realSurvLegacyConHi.m = matrix(NA,nrow=Z,ncol=n)
  realSurvLegacyConLo.m = matrix(NA,nrow=Z,ncol=n)
  realLivingCDD.m = matrix(NA,nrow=Z,ncol=n)
  realLegacyCDD.m = matrix(NA,nrow=Z,ncol=n)
  for(z in 1:Z) {							# Prepare model matrices for predictions
    SBbysp = matrix(NA,nrow=n,ncol=SK)
    for(i in 1:n) {SBbysp[i,] = SB.m[z,sp.info2$sp[j],]}
    QUADbycensus = c()
    for(i in 1:n) {QUADbycensus[i] = quad.m[z,plotnum[i]]}
    predSurvHi = inv.logit((xHi %*% beta.m[z,]) + rowSums(SXHi*SBbysp) + QUADbycensus)	# Predicted survival in high living conspecific density situation
    predSurvLo = inv.logit((xLo %*% beta.m[z,]) + rowSums(SXLo*SBbysp) + QUADbycensus)	# Predicted survival in low living conspecific density situation
    LpredSurvHi = inv.logit((xLHi %*% beta.m[z,]) + rowSums(SXLHi*SBbysp) + QUADbycensus)	# Predicted survival in high legacy conspecific density situation
    LpredSurvLo = inv.logit((xLLo %*% beta.m[z,]) + rowSums(SXLLo*SBbysp) + QUADbycensus)	# Predicted survival in low legacy conspecific density situation
    realSurvLivingConHi.m[z,] = predSurvHi^(1/5)				# Annual survival at 'high' living conspecific density (survival in presence of a standardized size/dist living conspecific at mean tree density -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
    realSurvLivingConLo.m[z,] = predSurvLo^(1/5)				# Annual survival at zero living conspecific density (survival in presence of no standardized size/dist living conspecific at mean tree density -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
    realSurvLegacyConHi.m[z,] = LpredSurvHi^(1/5)				# Annual survival at 'high' legacy conspecific density (survival in presence of a standardized size/dist legacy conspecific at mean tree density -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
    realSurvLegacyConLo.m[z,] = LpredSurvLo^(1/5)				# Annual survival at zero legacy conspecific density (survival in presence of no standardized size/dist legacy conspecific at mean tree density -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
    realLivingCDD.m[z,] = (predSurvHi^(1/5)) - (predSurvLo^(1/5))		# Calculate living CDD-HDD by comparing predicted survival in a stand with a standardized size/dist living conspecific to predicted survival in a stand with only living heterospecifics (each at mean tree density) -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
    realLegacyCDD.m[z,] = (LpredSurvHi^(1/5)) - (LpredSurvLo^(1/5))	# Calculate legacy CDD-HDD by comparing predicted survival in a stand with a standardized size/dist legacy conspecific to predicted survival in a stand with only legacy heterospecifics (each at mean tree density) -- taken to 1/5 power to calculate in annual survival (census interval was 5 yrs)
  }
  realSurvLivingConHi.sp[[j]] = apply(realSurvLivingConHi.m, 1, mean)	# Calculate average predictive comparison for living CDD-HDD (by averaging over all individuals for each posterior sample)
  realSurvLivingConLo.sp[[j]] = apply(realSurvLivingConLo.m, 1, mean)	# Calculate average predictive comparison for living CDD-HDD (by averaging over all individuals for each posterior sample)
  realSurvLegacyConHi.sp[[j]] = apply(realSurvLegacyConHi.m, 1, mean)	# Calculate average predictive comparison for living CDD-HDD (by averaging over all individuals for each posterior sample)
  realSurvLegacyConLo.sp[[j]] = apply(realSurvLegacyConLo.m, 1, mean)	# Calculate average predictive comparison for living CDD-HDD (by averaging over all individuals for each posterior sample)
  realLivingCDD.sp[[j]] = apply(realLivingCDD.m, 1, mean)	# Calculate average predictive comparison for living CDD-HDD (by averaging over all individuals for each posterior sample)
  realLegacyCDD.sp[[j]] = apply(realLegacyCDD.m, 1, mean)	# Calculate average predictive comparison for legacy CDD-HDD (by averaging over all individuals for each posterior sample)
}
end.time = Sys.time()
(duration = end.time - begin.time)

realSurvLivingConHi.sp = do.call('cbind', realSurvLivingConHi.sp)
realSurvLivingConLo.sp = do.call('cbind', realSurvLivingConLo.sp)
realSurvLegacyConHi.sp = do.call('cbind', realSurvLegacyConHi.sp)
realSurvLegacyConLo.sp = do.call('cbind', realSurvLegacyConLo.sp)
realLivingCDD.sp = do.call('cbind', realLivingCDD.sp)
realLegacyCDD.sp = do.call('cbind', realLegacyCDD.sp)

#save(realSurvLivingConHi.sp, realSurvLivingConLo.sp, realSurvLegacyConHi.sp, realSurvLegacyConLo.sp, realLivingCDD.sp, realLegacyCDD.sp, sp.info2, file = "ForestGEO_Wabikon_Legacy_CDD_LogisticSurvival_cdd_matrices_test.RData")
