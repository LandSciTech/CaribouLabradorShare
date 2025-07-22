#Example combining Labrador observed and simulated data.

#devtools::load_all(path = "../caribouMetrics/")
#devtools::install_github("LandSciTech/bboudata")
#devtools::install_github("LandSciTech/bboutoolsMultiPop")
#devtools::install_github("LandSciTech/caribouMetrics",ref="BbouIntegration")

###############################
#Initial setup
library(readxl)
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())
library(RColorBrewer)
library(bboutools)
library(mcmcr)
#library(caribouMetrics) #Note this must be installed from the BbouIntegration branch
devtools::load_all(path = "../caribouMetrics/")


figDir ="./figs/"
dir.create(figDir,recursive=T)
dia_shp <- 23
err_col <- "grey50"
baseDir <- "."

#Note - set niters to 100 to run quickly when testing. Set to 1000 for complete results.
niters <- 100

N0 <- 10000
scns=list()
scns$lQuantile=0.99
correlateRates = T #Force correlation among demographic rates to examine extreme cases

bbouResultFile = "../CaribouLabradorShare/data/bbouResultsLabrador.Rds"
scns$obsAnthroSlope = 0 #set NA to use bboutools nimble model, set number to use jags model with informative priors
scns$projAnthroSlope = 0 #set NA to use bboutools nimble model, set number to use jags model with informative priors
eParsIn = list()
eParsIn$collarOnTime=4
eParsIn$collarOffTime=4
eParsIn$collarNumYears=4
labFontSize = 10; breakInterval=5

#devtools::load_all(path = "../caribouMetrics/")
simBig <- getSimsInitial(bbouResultFile,cPars = list(correlateRates=correlateRates),forceUpdate=T)

###############
#Example monitoring scenario - repeat everything that was done again
#Kara - you could replace this with some other monitoring scenario.
cowCounts <- subset(simBig$recruit_data,!is.na(Cows)&(Year>=2019),select=c("PopulationName","Year", "Cows","CowsBulls", "UnknownAdults","Yearlings"))
simStartYr <- max(cowCounts$Year)+2
cowCounts$Year <- simStartYr + cowCounts$Year - min(cowCounts$Year)

pops<- unique(cowCounts$PopulationName)
firstPop <- T
for(p in pops){
  #p = pops[1]
  startBit <- subset(simBig$surv_data,(!is.na(Mortalities))&(PopulationName==p))
  startBit <- startBit[order(startBit$Year,startBit$Month),]
  startBit$nextStartTotal <- c(startBit$StartTotal[2:nrow(startBit)],NA)
  startBit$numStarts <- startBit$nextStartTotal-startBit$StartTotal+startBit$Mortalities+startBit$Malfunctions
  startBit$numStarts[is.na(startBit$numStarts)] <- 0
  startBit$numStarts[1] <- startBit$StartTotal[1]
  startBit <- startBit  %>% group_by(Annual,PopulationName) %>% summarize(numStarts=sum(numStarts))
  if(firstPop){
    freqStartsByYr <- startBit
    firstPop <- F
  }else{
    freqStartsByYr <- rbind(freqStartsByYr,startBit)
  }
}
freqStartsByYr$Annual <- as.numeric(as.factor(freqStartsByYr$Annual))
freqStartsByYr$Year <- simStartYr + freqStartsByYr$Annual - min(freqStartsByYr$Annual)

#Needed inputs are cowCounts and freqStartsByYear tables.
eParsIn$cowCounts <- cowCounts;eParsIn$freqStartsByYear <- freqStartsByYr

scns$obsYears<- max(freqStartsByYr$Year)-min(simBig$recruit_data$Year)
scns$startYear <- min(simBig$summary$Year)
scns$projYears <- max(simBig$summary$Year)-scns$obsYears-scns$startYear

#########################################################
#Example simulation
#devtools::load_all(path = "../caribouMetrics/")
posteriorResult = caribouMetrics:::runScnSet(scns,eParsIn,simBig,printProgress=F,niters=niters,nthin=10)

recPosterior =  plotRes(posteriorResult, "Recruitment", lowBound=-0.1,highBound = 1.2,
                        breakInterval=breakInterval,
                        labFontSize=labFontSize)
plot(recPosterior)
ggsave(paste0(baseDir,"/figs/eRec.png"),
       width = 9.6*0.779, height = 4, units = "in",
       dpi = 1200)

survPosterior =  plotRes(posteriorResult, "Adult female survival", lowBound=0,
                         breakInterval=breakInterval,
                         labFontSize=labFontSize)
plot(survPosterior)
ggsave(paste0(baseDir,"/figs/eSurv.png"),
       width = 9.6*0.779, height = 4, units = "in",
       dpi = 1200)


lambdaPosterior =  plotRes(posteriorResult, "Population growth rate", lowBound=0,
                           breakInterval=breakInterval,
                           labFontSize=labFontSize)+
  ylim(c(0, 1.5))
plot(lambdaPosterior)
ggsave(paste0(baseDir,"/figs/eLam.png"),
       width = 9.6*0.779, height = 4, units = "in",
       dpi = 1200)

lambdaPosterior =  plotRes(posteriorResult, "Expected growth rate", lowBound=0,
                           breakInterval=breakInterval,
                           labFontSize=labFontSize)+
  ylim(c(0, 1.5))
plot(lambdaPosterior)
ggsave(paste0(baseDir,"/figs/eLamBar.png"),
       width = 9.6*0.779, height = 4, units = "in",
       dpi = 1200)

