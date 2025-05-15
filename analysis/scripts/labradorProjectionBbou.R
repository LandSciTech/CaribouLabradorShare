########
#Before running, install these packages from github
#devtools::install_github("LandSciTech/bboudata")
#devtools::install_github("LandSciTech/bboutoolsMultiPop")
#devtools::install_github("LandSciTech/caribouMetrics",ref="BbouIntegration")
#devtools::install_github("LandSciTech/CaribouDemographyBasicApp")
library(bboutools)
library(caribouMetrics)
library(CaribouDemographyBasicApp)

#These R packages can be installed from cran if you don't have them.
library(readxl)
library(tidyverse)
library(mcmcr)

figDir ="./figs/"
dir.create(figDir,recursive=T)

bbouResultFile = "./data/bbouResultsLabrador.Rds"
results = readRDS(bbouResultFile)

names(results)

str(results$surv_fit)

results$surv_fit$data

subset(results$surv_fit$data,!is.na(MortalitiesCertain))


surv_fit <- results$surv_fit
surv_pred <- bb_predict_survival(surv_fit)
bb_plot_year_survival(surv_pred) +
  expand_limits(y = c(0, 1))

recruit_fit <- results$recruit_fit
bb_plot_year_recruitment(recruit_fit)

#####################
#project using bboutools

predict_lambda <- bb_predict_growth(survival = surv_fit, recruitment = recruit_fit)

png(paste0(figDir,"/bbouLambda.png"),
    height = 6, width = 10.56, units = "in",res=600)
bb_plot_year_growth(predict_lambda)
dev.off()

predict_change <- bb_predict_population_change(survival = surv_fit, recruitment = recruit_fit)
bb_plot_year_population_change(predict_change)

######################
#project using caribouMetrics - mcmc samples
surv_pred = bb_predict_survival (surv_fit,year=T,month=F,conf_level=F)
rec_pred =bb_predict_calf_cow_ratio(recruit_fit,year=T,conf_level=F)

#devtools::load_all(path = "../caribouMetrics/")
outmcmc = caribouPopSimMCMC(results$parTab,rec_pred,surv_pred,initYear = 2024)

outmcmc$Anthro = NA
outmcmc$fire_excl_anthro=NA

########################
#project using caribouMetrics - from summary parameters
parTab = results$parTab
write.table(parTab,"./data/Lab_pop_data.csv",row.names=F,sep=",")

pops=unique(parTab$pop_name)
numSteps = 20;numPops=length(unique(outmcmc$id))
for(i in 1:length(pops)){
  #i=2
  p = pops[i]
  cpars = subset(parTab,pop_name==p)
  scn_nm = p;addl_params=list()
  ddBit = doSim(numSteps, numPops, unique(parTab$N0), cpars$R_bar, cpars$S_bar, cpars$R_sd, cpars$S_sd,
                R_iv_cv=cpars$R_iv_cv, S_iv_cv=cpars$S_iv_cv,
                R_iv_sd = cpars$R_iv_sd,S_iv_sd=cpars$S_iv_sd,
             scn_nm=scn_nm, addl_params=NULL,type="logistic")
  if(i==1){
    dd = ddBit
  }else{
    dd=rbind(dd,ddBit)
  }
}

#################
#now compare various projections
ddSamp = subset(dd,type=="samp")
ddSamp$PopulationName=ddSamp$scn

summarySet = subset(ddSamp,select=c(N,id,time,PopulationName,R_t,S_t,lambda))
summarySet$method="summaries"

mcmcSet = subset(outmcmc,select=c(N,id,time,PopulationName,R_t,S_t,lambda))
mcmcSet$method="mcmc"

dAll = rbind(summarySet,mcmcSet)
base=ggplot(subset(dAll,id<=35),aes(x=time,y=N,group=id))+geom_line(alpha=0.3)+facet_grid(PopulationName~method)
print(base)

base=ggplot(subset(dAll,id<=100),aes(x=time,y=R_t,group=id))+geom_line(alpha=0.3)+facet_grid(PopulationName~method)
print(base)

base=ggplot(subset(dAll,id<=100),aes(x=time,y=S_t,group=id))+geom_line(alpha=0.3)+facet_grid(PopulationName~method)
print(base)

dAllGrp = dAll |> group_by(PopulationName,id,method) |> summarize (
  R_bar=mean(R_t),S_bar=mean(S_t),lambda=exp(mean(log(lambda))),
  R_Annual = sd(R_t),S_Annual = sd(S_t))

base=ggplot(dAllGrp,aes(x=method,y=R_bar))+geom_violin()+geom_boxplot()+facet_wrap(~PopulationName)
print(base) #confirming similarity of distributions

base=ggplot(dAllGrp,aes(x=method,y=S_bar))+geom_violin()+geom_boxplot()+facet_wrap(~PopulationName)
print(base) #confirming similarity of distributions

base=ggplot(dAllGrp,aes(x=method,y=lambda))+geom_violin()+geom_boxplot()+facet_wrap(~PopulationName)
print(base) #confirming similarity of distributions

base=ggplot(dAllGrp,aes(x=method,y=R_Annual))+geom_violin()+geom_boxplot()+facet_wrap(~PopulationName)
print(base) ##confirming similarity of distributions

base=ggplot(dAllGrp,aes(x=method,y=S_Annual))+geom_violin()+geom_boxplot()+facet_wrap(~PopulationName)
print(base) ##confirming similarity of distributions

#In summary - I believe summary method adequately matches the mcmc method.
