################################################################################
############################# SETTING THE SCENE ################################
################################################################################

# ------------------------------ EXPLANATION  ----------------------------------
# This is the R script used for the data analysis presented in the manuscript
# Rehren, J., Pennino, M.G., Coll, M., Jiddawi, N., Muhando, C. 
# Supporting spatial management of data-poor small-scale fisheries with a 
# Bayesian approach. Frontiers in Marine Science. (submitted october 2020)

# ------------------------------------------------------------------------------

# --------------------------- LOAD PACKAGES AND FUNCTIONS ----------------------
# 1.) CLEAR SCREEN AND SET WORKING DIRECTORY ----
ls()
rm(list=ls())

# 2.) LOAD PACKAGES ----

library(dplyr)
library(tidyverse)
library(INLA)

# ------------------------------------------------------------------------------

################################################################################
################################# DATA ANALYSIS ################################
################################################################################

# 1.) LOAD THE DATA FRAME DATA ----

df <- read.table("Data/Data frame_Frontiers2020_Rehren.csv", 
                  sep=";", header=T)
head(df)

# remove space after species names
df$latinName <- trimws(df$latinName)

# 2.) LOG-TRANSFORMATION OF WPUE ----

df$logwpue <- log(df$wpue)

# 3.) SUBSET THE DATA BY SPECIES ----

dfL <- split(df, df$latinName, drop=T)


# 4.) REMOVE NA's ----

datL <- list()
for (i in names(dfL)){
  datL[[i]] <- dfL[[i]][!is.na(dfL[[i]]$wpue), ]
  datL[[i]] <- datL[[i]][!is.na(datL[[i]]$seagrass), ] # remove missing data
  datL[[i]] <- datL[[i]][!is.na(datL[[i]]$depth), ] # remove missing data
}

# 5.) OUTLIERS ----
# The following data points were excluded from analysis given their
# influence / leverage, which was determined using the conditional predictive
# ordinate. 

datL$`Scarus ghobban`       <- datL$`Scarus ghobban`[-48,]
datL$`Lethrinus lentjan`    <- datL$`Lethrinus lentjan`[-327,]
datL$`Lethrinus lentjan`    <- datL$`Lethrinus lentjan`[-355,]
datL$`Lethrinus mahsena`    <- datL$`Lethrinus mahsena`[-67,]
datL$`Lutjanus fulviflamma` <- datL$`Lutjanus fulviflamma`[-255,]
ind <- which(datL$`Siganus sutor`$wpue>9)
datL$`Siganus sutor`        <- datL$`Siganus sutor`[-ind,]

# 6.) MESH ----

locs <- lapply(datL, function (x) cbind(x$X.utm, x$Y.utm))

boundaries <- lapply(locs, function (x) inla.nonconvex.hull(x))

maxEdges <- lapply(datL, function (x) {
  diff(range(x$Y.utm))/15
}
)

meshes <- list()
for (i in names(boundaries)){
  meshes[[i]] <- inla.mesh.2d(boundary = boundaries[[i]], 
                              max.edge=c(1, 5) * maxEdges[[i]], 
                              cutoff = maxEdges[[i]] / 5)
}

# 7.) A MATRIX ----

A <- mapply(inla.spde.make.A, mesh=meshes, loc = locs)

# 8.) SPDE ----

# Sigma values for each species used in the inla.make.spde function
names(datL)
sigma <- list(2, 1.7, 2, 1.7, 1.7, 2)

spdes <- list()
for (i in seq_along(meshes)){
  spdes[[i]] <- inla.spde2.pcmatern(mesh = meshes[[i]], 
                                    prior.range = c(1000, 0.05), 
                                    prior.sigma = c(sigma[[i]], 0.05))
}
names(spdes) <- names(datL) 

# 9.) W INDICES ----

wIndices <- lapply(spdes, function (x) inla.spde.make.index(name = 'w', 
                                                            n.spde = x$n.spde)) 

# 10.) COVARIATE MATRIX ----

# We transform the categorical covariates into factors

for (i in names(datL)){
  datL[[i]]$freef             <- as.factor(datL[[i]]$reef)
  datL[[i]]$fseason           <- as.factor(datL[[i]]$season)
  datL[[i]]$flunarCycle       <- as.factor(datL[[i]]$lunarCycle)
  datL[[i]]$fmetier           <- as.factor(datL[[i]]$metier)
}

xm <- list()
for (i in names(datL)){
  xm[[i]] <- model.matrix(~ depth.int.std + seagrass.int.std + 
                            temperature.std + freef + 
                            flunarCycle + fseason + fmetier, data = datL[[i]])
}


dims <- lapply(xm, function(x) dim(x))
n <- lapply(datL, function(x) nrow(x))

x <- list()
for (i in names(xm)){
  x[[i]] <- data.frame(depth.int.std         = xm[[i]][,2],
                       seagrass.int.std      = xm[[i]][,3],
                       temperature.std       = xm[[i]][,4],
                       freef                 = xm[[i]][,5],
                       flunarCycle           = xm[[i]][,6:8],
                       fseason               = xm[[i]][,9],
                       fmetier               = xm[[i]][,10:dims[[i]][2]])
}

# 11.) STACK ----

stacks <- list()

for (i in names(datL)){
  stacks[[i]] <- inla.stack(
    tag = "est",
    data = list(y = datL[[i]]$logwpue),  
    A = list(1, 1, A[[i]], 1),                  
    effects = list(   
      beta0 = rep(1, n[[i]]),
      X = x[[i]],
      w = wIndices[[i]],
      metier = datL[[i]]$fmetier)
  )
}

# 12.) MODEL FORMULAR ----

listForm <- list()

listForm[[1]] <- c("y ~ -1 + beta0 + depth.int.std + seagrass.int.std + 
                   temperature.std + flunarCycle.flunarCycle4 + fseason",
                   "f(metier, model='iid')","f(w, model=spdes$'Leptoscarus vaigiensis')")

listForm[[2]] <- c("y ~ -1 + beta0 + depth.int.std + flunarCycle.flunarCycle4", 
                   "f(metier, model='iid')", "f(w, model = spdes$'Lethrinus lentjan')")

listForm[[3]] <- c("y ~ -1 + beta0 + seagrass.int.std + temperature.std +  freef", 
                   "f(metier, model='iid')",
                   "f(w, model = spdes$'Lethrinus mahsena')")

listForm[[4]] <- c("y ~ -1 + beta0 + depth.int.std +  seagrass.int.std", 
                   "f(metier, model='iid')", 
                   "f(w, model = spdes$'Lutjanus fulviflamma')")

listForm[[5]] <- c("y ~ -1 + beta0 + temperature.std +  
                   flunarCycle.flunarCycle4", "f(metier, model='iid')",  
                   "f(w, model = spdes$'Scarus ghobban')")

listForm[[6]] <- c("y ~ -1 + beta0 + flunarCycle.flunarCycle2", 
                   "f(metier, model='iid')", "f(w, model = spdes$'Siganus sutor')")

names(listForm) <- names(datL)

form <- list()
for (i in names(listForm)){
  form[[i]] <- paste(listForm[[i]], collapse = " + ")
}

# 13.) INLA ----

modList <- list()
for (i in names(datL)){
    modList[[i]] <- inla(formula = eval(parse(text=form[[i]])),
                              family = "gaussian", 
                              data = inla.stack.data(stacks[[i]]),
                              control.compute = list(dic = TRUE, waic=T, cpo=T),
                              control.predictor = list(A = inla.stack.A(stacks[[i]]))) 
  }

