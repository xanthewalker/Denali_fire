rm(list=ls())

library(nlme)
library(ggplot2)
library(dplyr)
library(gtools)
library(glmmTMB)
library(MuMIn)

##Burn Depth, Residual SOL, Proportion SOL, ALT 
#####

setwd("G:\\Other computers\\My Laptop\\Documents\\NAU\\Data\\Denali\\Combustion")


site<-read.csv("Site_Denali_ALL_02122021.csv") 
names(site)
head(site)


comb=read.csv("gridpoint_meter_C_combustion_04232020.csv")
names(comb)              
head(comb)
siteC=merge(site, comb, by="gridpoint", all=TRUE)
head(siteC)


sol=read.csv("soil.thaw.depth.csv")
head(sol)

siteD=merge(siteC, sol, by=c("gridpoint", "m"), all=TRUE)

siteD=mutate(siteD, moisture=ifelse(moisture=="subxeric", "mesic", 
                                    ifelse(moisture=="subxeric-mesic", "mesic", 
                                           ifelse(moisture=="xeric-subxeric", "mesic", 
                                                  ifelse(moisture=="mesic", "mesic", 
                                                         ifelse(moisture=="mesic-subhygric", "mesic-subhygric",
                                                                ifelse(moisture=="subhygric", "subhygric", NA)))))))

head(siteD)
names(siteD)

siteX=dplyr::select(siteD, grid.x, gridpoint, m, sol.y,  burn.depth.y,  density_spruce, burn.year1, new.td)
siteX=na.omit(siteX)


siteX=mutate(siteX, burn.depth.y=ifelse(burn.depth.y<0, 0, burn.depth.y))
min(siteX$sol.y)

siteX=mutate(siteX, prefireSOL=sol.y+burn.depth.y)
siteX=mutate(siteX, propSOL=burn.depth.y/prefireSOL)

siteX<-mutate(siteX, group1=quantcut(density_spruce, q=3))
levels(siteX$group1)
siteX<-mutate(siteX, Density=plyr::revalue(group1, c("0"="low", "(0,0.0249]"="mid", "(0.0249,0.229]"="high")))
head(siteX)

siteX=subset(siteX, grid.x=="WC"|grid.x=="ET")
siteX[] <- lapply(siteX, function(x) if(is.factor(x)) factor(x) else x)

head(siteX)
head(siteX)
siteX=subset(siteX, burn.year1=="burn")



##burn depth model
y1<-lme(data=siteX, burn.depth.y~prefireSOL*density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y1)

y<-lme(data=siteX, burn.depth.y~prefireSOL+density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y)
anova(y, y1)
plot(y)
AIC(y)
r.squaredGLMM(y)

y4<-lme(data=siteX, burn.depth.y~density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y4)
AIC(y4)

y6<-lme(data=siteX, burn.depth.y~prefireSOL, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
AIC(y6)
anova(y, y6)

y7<-lme(data=siteX, burn.depth.y~1, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y6)
anova(y1, y2)

library(MuMIn)
model.sel(y1, y,  y4, y6, y7)
summary(y)
r.squaredGLMM(y)

##residual SOL depth
y<-lme(data=siteX, sol.y~prefireSOL*density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y)


y1<-lme(data=siteX, sol.y~prefireSOL+density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y1)
plot(y1)
anova(y, y1)


y2<-lme(data=siteX, sol.y~prefireSOL, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
anova(y1, y2)
plot(y)
anova(y, y1)

plot(allEffects(y1))
summary(y1)
r.squaredGLMM(y1)

##proportion SOL combusted

y<-lme(data=siteX, propSOL~prefireSOL*density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y)
plot(y)


y1<-lme(data=siteX, propSOL~prefireSOL+density_spruce, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y1)
anova(y, y1)


y2<-lme(data=siteX, propSOL~prefireSOL, random=~1|grid.x/gridpoint, na.action=na.omit, method="ML")
summary(y1)

anova(y2, y1)

plot(allEffects(y1))
summary(y1)
r.squaredGLMM(y1)

##thaw depth
names(siteX)

td=read.csv("thaw_depth.csv")
names(td)      
str(td)
head(td)
td=mutate(td, new.td=thaw.d)
forsoltd=dplyr::select(td, grid, gridpoint, new.td, sol)

head(td)
TD.all=read.csv("Denali_ThawDepths_ALL_05252021.csv")
head(TD.all)
TD=filter(TD.all, grid %in% c("ET", "WC", "HC"))
head(TD)
TD$thaw.d=as.numeric(TD$thaw.d)

sol=read.csv("soil.thaw.depth.csv")
head(sol)
levels(sol$gridpoint)
sol=subset(sol, burn=="burn")
sol1=filter(sol, grid %in% c("ET", "WC"))
sol2=dplyr::select(sol1, grid, gridpoint, new.td, sol)
cuz=rbind(forsoltd, sol2)
head(cuz)
cuz=na.omit(cuz)
cuz1=mutate(cuz, new.td=ifelse(new.td==sol, NA, new.td))
head(cuz1)
cuz1=na.omit(cuz1)

site<-read.csv("Site_Denali_ALL_02122021.csv") 
names(site)
head(site)
levels(site$burn.year1)

site1=dplyr::select(site, gridpoint, density_spruce)


head(siteM)

td1=merge(cuz1, site1, by=c("gridpoint"))
names(td1)

td1=dplyr::select(td1, grid, gridpoint, sol, new.td, density_spruce)
head(td1)

td2=na.omit(td1)


td2<-mutate(td2, group1=quantcut(density_spruce, q=3))
td2<-mutate(td2, Density=plyr::revalue(group1, c("0"="low", "(0,0.0249]"="mid", "(0.0249,0.229]"="high")))


td.mod1 <-glmer(new.td~sol*density_spruce+ (1|grid/gridpoint), 
                family=Gamma(link="log"), nAGQ=0, data = td2, na.action=na.omit)

td.mod <-glmer(new.td~sol+density_spruce+ (1|grid/gridpoint), 
               family=Gamma(link="log"), nAGQ=0, data = td2, na.action=na.omit)


anova(td.mod, td.mod1)

td.mod3 <-glmer(new.td~sol  + (1|grid/gridpoint), 
                family=Gamma(link="log"), nAGQ=0, data = td2, na.action=na.omit)

anova(td.mod, td.mod3)

td.mod4 <-glmer(new.td~density_spruce  + (1|grid/gridpoint), 
                family=Gamma(link="log"), nAGQ=0, data = td2, na.action=na.omit)

anova(td.mod, td.mod4)


model.sel(td.mod,  td.mod3, td.mod4)
plot(allEffects(td.mod))
summary(td.mod)
r.squaredGLMM(td.mod)



#####
#####  Total C Combustion, Prop total C combusted
#####
resid=read.csv("measured_residual_C_BD_04232020.csv")
head(resid)

resid=mutate(resid, total.gC.m2=total.gC.m2/1000)

comb=read.csv("gridpoint_meter_C_combustion_04232020.csv")
names(comb)              
head(comb)
levels(comb$gridpoint)

comb=mutate(comb, carbon_comb=carbon_comb/1000)




try=merge(resid, comb, by=c("gridpoint", "m"))
head(try)
names(try)


try1=mutate(try, prefireSOLC=carbon_comb+total.gC.m2)
head(try1)
try1=mutate(try1, prop.SOL.C.con=carbon_comb/prefireSOLC)            
head(try)
names(try)
try2=na.omit(try1)


#
names(siteD)
siteD=merge(siteX, try2, by=c("gridpoint", "m"))
head(siteD)

names(siteD)
soil=dplyr::select(siteD, grid.x, gridpoint, burn.year1,  carbon_comb, total.gC.m2, density_spruce)
head(soil)


tree=dplyr::select(site, grid, gridpoint, density_spruce, tree.carbon.cons.m2, tree.bio.carbon.m2, burn.year1)
tree <- tree %>%
  mutate(tree.carbon.cons.m2 = if_else(is.na(tree.carbon.cons.m2), 0, tree.carbon.cons.m2), 
         tree.bio.carbon.m2=if_else(is.na(tree.bio.carbon.m2), 0, tree.bio.carbon.m2))
head(tree)

tree=na.omit(tree)
tree.burn=subset(tree, burn.year1=="burn")
tree.burn=subset(tree.burn, grid=="ET"|grid== "WC")

soil.burn=subset(soil, burn.year1=="burn")
soil.burn=na.omit(soil.burn)
soil.burn=mutate(soil.burn, prefireSOLC=carbon_comb+total.gC.m2)
soil.burn=mutate(soil.burn, prop.SOL.C.con=carbon_comb/prefireSOLC)  


soil.burn1=soil.burn%>%
  dplyr::group_by(grid.x, gridpoint)%>%
  dplyr::summarize(
    SOL.c.cons= mean(carbon_comb), 
    prefire.SOL.c=mean(prefireSOLC))

head(soil.burn1)
names(tree.burn)

soil.burn1=mutate(soil.burn1, grid=grid.x)

siteZ=merge(tree.burn, soil.burn1, by=c("grid", "gridpoint"))
head(siteZ)

siteZ=mutate(siteZ, SOL.c.cons=ifelse(SOL.c.cons>prefire.SOL.c, prefire.SOL.c, SOL.c.cons))

siteZ=mutate(siteZ, total.C.cons=SOL.c.cons+tree.carbon.cons.m2)
siteZ=mutate(siteZ, total.C.prefire=prefire.SOL.c+tree.bio.carbon.m2)
siteZ=mutate(siteZ, prop.C.con=total.C.cons/total.C.prefire)
siteZ=mutate(siteZ, prop.SOL.C.con=SOL.c.cons/total.C.cons)
head(siteZ)

siteZ<-mutate(siteZ, group1=quantcut(density_spruce, q=3))
levels(siteZ$group1)
siteZ<-mutate(siteZ, Density=plyr::revalue(group1, c("0"="low", "(0,0.0249]"="mid", "(0.0249,0.229]"="high")))
head(siteZ)


##total C combusted
y<-lme(data=siteZ, total.C.cons~total.C.prefire*density_spruce, random=~1|grid, na.action=na.omit, method="ML")
summary(y)

y1<-lme(data=siteZ, total.C.cons~total.C.prefire+density_spruce, random=~1|grid, na.action=na.omit, method="ML")
summary(y1)

anova(y, y1)


y6<-lme(data=siteZ, total.C.cons~total.C.prefire, random=~1|grid, na.action=na.omit, method="ML")
AIC(y6)

y2<-lme(data=siteZ, total.C.cons~density_spruce, random=~1|grid, na.action=na.omit, method="ML")
AIC(y2)

y3<-lme(data=siteZ, total.C.cons~1, random=~1|grid, na.action=na.omit, method="ML")
AIC(y3)

model.sel(y, y6, y1, y2, y3)
plot(y1)
plot(allEffects(y1))
summary(y1)
r.squaredGLMM(y1)

##proportion of total C combusted

y4 <-glmer(prop.C.con~total.C.prefire*density_spruce  + (1|grid), 
           family=Gamma(link="log"), nAGQ=0, data = siteZ, na.action=na.omit)
plot(y4)
AIC(y4)
summary(y3)

y3 <-glmer(prop.C.con~total.C.prefire+density_spruce  + (1|grid), 
           family=Gamma(link="log"), nAGQ=0, data = siteZ, na.action=na.omit)
plot(y3)
AIC(y3)
summary(y3)


y2 <-glmer(prop.C.con~density_spruce  + (1|grid), 
           family=Gamma(link="log"), nAGQ=0, data = siteZ, na.action=na.omit)
plot(y3)
AIC(y3)


y1 <-glmer(prop.C.con~total.C.prefire  + (1|grid), 
           family=Gamma(link="log"), nAGQ=0, data = siteZ, na.action=na.omit)
plot(y3)
AIC(y3)

y <-glmer(prop.C.con~1 + (1|grid), 
          family=Gamma(link="log"), nAGQ=0, data = siteZ, na.action=na.omit)
plot(y3)
AIC(y3)


model.sel(y4, y3, y2, y1, y)

library(DHARMa)
plot(simulateResiduals(y3))
r.squaredGLMM(y3)

#####
#####  Turnover and natural regeneration
#####
setwd("G:\\Other computers\\My Laptop\\Documents\\NAU\\Data\\Denali\\R Code")

to=read.csv("DNP_turnover_04212021.csv")
names(to)
head(to)
site=read.csv("Site_Denali_ALL_04232020.csv")
names(site)
head(site)
site.to=merge(site, to, by="gridpoint")
names(site.to)

levels(site.to$moisture)
library(dplyr)

site.to$moisture<-factor(site.to$moisture, levels=c("xeric", "subxeric", "mesic-subxeric", "mesic", "mesic-subhygric", "subhygric"))
site.to=mutate(site.to, moisture=ifelse(moisture=="subxeric", "dry", 
                                        ifelse(moisture=="subxeric-mesic", "dry", 
                                               ifelse(moisture=="xeric-subxeric", "dry", 
                                                      ifelse(moisture=="mesic", "dry", 
                                                             ifelse(moisture=="mesic-subhygric", "moist",
                                                                    ifelse(moisture=="subhygric", "wet", NA)))))))

levels(site.to$moisture)
site.to[] <- lapply(site.to, function(x) if(is.factor(x)) factor(x) else x)
library(ggplot2)
names(t9)

t9=dplyr::select(site.to, total, mean.org.tot,density_spruce,burn.depth, grid)
t9=na.omit(t9)
library(nlme)
M1<-lme(total~ mean.org.tot+burn.depth+density_spruce,random=~ 1|grid, method  ="ML", data=t9, na.action=na.exclude)
plot(M1)
summary(M1)
library(MuMIn)
r.squaredGLMM(M1)


M2<-lme(total~ mean.org.tot+density_spruce,random=~ 1|grid, method  ="ML", data=t9, na.action=na.exclude)
plot(M2)
anova(M1, M2)

M3<-lme(total~ mean.org.tot,random=~ 1|grid, method  ="ML", data=t9, na.action=na.exclude)
plot(M3)
anova(M3, M2)

M4<-lme(total~ 1,random=~ 1|grid, method  ="ML", data=t9, na.action=na.exclude)
plot(M3)
anova(M3, M4)

summary(M1)
r.squaredGLMM(M1)
plot(allEffects(M1))
library(effects)


plot(allEffects(M1))
library(effects)
#####
#######Natural regeneration
#####

site<-read.csv("Site_Denali_ALL_04232020.csv")
names(site)


data7<-mutate(site, burn.year1=ifelse(burn.year=="2015", "burn", 
                                      ifelse(burn.year=="2013", "burn", "unburned")))


data7=subset(data7, grid=="WC"|grid=="ET")
data8=subset(data7, burn.year1=="burn")
data7<-mutate(data7, Density=plyr::revalue(group1, c("0"="low", "(0,0.0249]"="mid", "(0.0249,0.229]"="high")))
head(data7)
data7[] <- lapply(data7, function(x) if(is.factor(x)) factor(x) else x)

bd=read.csv("grid_direction_burndepth.csv")
head(bd)
bd1=bd%>%
  dplyr::group_by(gridpoint)%>%
  dplyr::summarise(burn.depth_good=mean(burn.depth))
count=bd%>%
  group_by(gridpoint)%>%
  dplyr::summarize(count=n())

head(bd1)
names(bd1)
names(data7)
data8=dplyr::select(data7, gridpoint, grid, burn.year:bulk.dens, Density, density_ALL, density_spruce, density_picmar, density_picgla, mean.org.tot, sol, picmar.seedling:spruce.seedling, picmar.sdlng.dens, picgla.sdlng.dens, spruce.seedlng.dens, burn.depth:burn.year1, tree.bio.m2:tree.carbon.cons.m2, picgla_tree.spec.bio.carbon.m2, picmar_tree.spec.bio.carbon.m2, picmar_tree.spec.carbon.cons.m2, picgla_tree.spec.carbon.cons.m2, bd.tot, bd.fine)
data9=merge(data8, bd1, by="gridpoint", all=TRUE)
head(data9)



data9=mutate(data9, burn.D=ifelse(burn.year1=="unburned", 0, burn.depth_good))
head(data9)

data9=mutate(data9, burn.D=ifelse(burn.D<0, 0, burn.D))

names(data9)
names(data7)

names(data9)
nat=dplyr::select(data9, density_spruce, mean.org.tot, bd.tot, burn.D, grid, Density, spruce.seedling, burn.year1)
nat=na.omit(nat)
head(nat)
nat=subset(nat, burn.year1=="burn")
head(nat)
levels(nat$grid)

data9=na.omit(nat)
names(data9)

###FIRST MODEL CONIFER SEEDLINGS 
##nibinom2 w/zero inflation
fit_binomS2=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                    data=data9,  
                    ziformula = ~1,  # allows structural zeros to depend on all variables
                    family=nbinom2) ## uses long link
summary(fit_binomS2)
drop1(fit_binomS2)

#### so results are on natural log scale##
####http://environmentalcomputing.net/interpreting-coefficients-in-glms/ 

#### nbinom 2 # no zero infl.
fit_binomS2b=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                     data=data9,   
                     ziformula = ~0,  # 
                     family=nbinom2) ## uses log link
summary(fit_binomS2b)
plot(simulateResiduals(fit_binomS2b))

summary(fit_binomS2b)
drop1(fit_binomS2b)

###### poisson w/zero inflations
fit_poiS1=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                  data=data9, 
                  ziformula = ~1, 
                  family=poisson)## 
summary(fit_poiS1)
library(DHARMa)
plot(simulateResiduals(fit_poiS1))
summary(fit_poiS1)
drop1(fit_poiS1)


###### poisson # no zero infl
fit_poiS1b=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                   data=data9,  
                   ziformula = ~0, 
                   family=poisson)## 
summary(fit_poiS1b)
plot(simulateResiduals(fit_poiS1b))
summary(fit_poiS1b)
drop1(fit_poiS1b)

#### nbinom1 w/ zero inflation
fit_bS1=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                data=data9, 
                ziformula = ~1, 
                family=nbinom1)## 
summary(fit_bS1)
plot(simulateResiduals(fit_bS1))
summary(fit_bS1)
drop1(fit_bS1)

#### nbinom1 - no zero inf.
fit_bS1b=glmmTMB(spruce.seedling ~mean.org.tot+density_spruce+burn.D+(1|grid),
                 data=data9,  
                 ziformula = ~0, 
                 family=nbinom1)
summary(fit_bS1b)
drop1(fit_bS1b)
plot(simulateResiduals(fit_bS1b))


model.sel(fit_binomS2,fit_binomS2b,fit_bS1,fit_bS1b,fit_poiS1,fit_poiS1b) 
plot(allEffects(fit_binomS2))
plot(allEffects(fit_binomS2b))
plot(allEffects(fit_bS1))
plot(allEffects(fit_bS1b))
plot(allEffects(fit_poiS1))
plot(allEffects(fit_poiS1b))

##fit_binomS2bmost parsimonious 


a=glmmTMB(spruce.seedling ~mean.org.tot+burn.D+density_spruce+(1|grid),
          data=data9,   
          ziformula = ~0,  # 
          family=nbinom2) 

summary(a)
drop1(a)
b=glmmTMB(spruce.seedling ~density_spruce+burn.D+(1|grid),
          data=data9,   
          ziformula = ~0,  # 
          family=nbinom2) 



c=glmmTMB(spruce.seedling ~mean.org.tot+burn.D+(1|grid),
          data=data9,   
          ziformula = ~0,  # 
          family=nbinom2) 


d=glmmTMB(spruce.seedling ~1+(1|grid),
          data=data9,  
          ziformula = ~0, 
          family=nbinom2)

e=glmmTMB(spruce.seedling ~density_spruce+(1|grid),
          data=data9,   
          ziformula = ~0,  # 
          family=nbinom2) 


model.sel(a, b, c, d, e) 



summary(a)
r.squaredGLMM(a)


#####
######Seed rain 
#####

site<-read.csv("Site_Denali_ALL_04232020.csv")
names(site)

seed<-read.csv("DNP_SeedTrap.csv")
names(seed)

seed$burn<-factor(seed$burn, levels=c("unburned", "burn"))
seed$date<-factor(seed$date, levels=c("Fall 2016", "Spring 2017", "Fall 2017", "Spring 2018", "Fall 2018"))

sum(seed$spruce == 0)  #Number of zeros
100 * sum(seed$spruce == 0) / nrow(seed)  
#53% are zeros :()

seed1=seed%>%
  filter(burn=="unburned")%>%
  summarize(mean=mean(density))

seed1=seed%>%
  filter(burn=="unburned")%>%
  group_by(gridpoint)%>%
  summarize(sum=sum(density))

seed2=mutate(seed1, sum3=sum/3)

seed1=seed%>%
  filter(burn=="burn")%>%
  select(density)%>%
  summarize_all(funs(mean, sd, min, max))

#
###names(seed)

seed.site<-merge(site, seed, by="gridpoint", all=TRUE)
head(seed.site)

seed.site=dplyr::select(seed.site, gridpoint, grid.x, spruce, density_ALL,density_spruce, date, burn, picmar_tree.spec.bio.m2, picgla_tree.spec.bio.m2, picmar_tree.spec.carbon.cons.m2, picgla_tree.spec.carbon.cons.m2)


seed.site <- seed.site %>%
  mutate(picmar_tree.spec.bio.m2 = if_else(is.na(picmar_tree.spec.bio.m2), 0, picmar_tree.spec.bio.m2), 
         picgla_tree.spec.bio.m2 = if_else(is.na(picgla_tree.spec.bio.m2), 0, picgla_tree.spec.bio.m2), 
         picmar_tree.spec.carbon.cons.m2 = if_else(is.na(picmar_tree.spec.carbon.cons.m2), 0, picmar_tree.spec.carbon.cons.m2),
         picgla_tree.spec.carbon.cons.m2 = if_else(is.na(picgla_tree.spec.carbon.cons.m2), 0, picgla_tree.spec.carbon.cons.m2))

head(seed.site)
seed.site=na.omit(seed.site)
library(glmmTMB)
##FIRST MODEL CONIFER SEEDLINGS 
#nibinom2 w/zero inflation
fit_binomS2=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                    data=seed.site,  
                    ziformula = ~1,  # allows structural zeros to depend on all variables
                    family=nbinom2) ## uses long link
##so results are on natural log scale##


#### nbinom 2 # no zero infl.
fit_binomS2b=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2) ## uses log link

#### poisson w/zero inflations
fit_poiS1=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                  data=seed.site, 
                  ziformula = ~1, 
                  family=poisson)## 

#### poisson # no zero infl
fit_poiS1b=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                   data=seed.site, 
                   ziformula = ~0, 
                   family=poisson)## 

###nbinom1 w/ zero inflation
fit_bS1=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                data=seed.site, 
                ziformula = ~1, 
                family=nbinom1)## 
#### nbinom1 - no zero inf.
fit_bS1b=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                 data=seed.site, 
                 ziformula = ~0, 
                 family=nbinom1)


model.sel(fit_binomS2,fit_binomS2b,fit_bS1,fit_bS1b,fit_poiS1,fit_poiS1b) 
##fit_binomS2b most parsimonious 


fit_binomS2b=glmmTMB(spruce~density_spruce*burn+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2) 

summary(fit_binomS2b)
AIC(fit_binomS2b)
r.squaredGLMM(fit_binomS2b)

fit_binomS2c=glmmTMB(spruce~density_spruce+burn+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2) 

summary(fit_binomS2c)
r.squaredGLMM(fit_binomS2c)
AIC(fit_binomS2c)
fit_binomS2d=glmmTMB(spruce~burn+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2) 


fit_binomS2e=glmmTMB(spruce~density_spruce+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2)

fit_binomS2f=glmmTMB(spruce~1+(1|grid.x)+(1|date),
                     data=seed.site,  
                     ziformula = ~0,  # 
                     family=nbinom2)

AIC(fit_binomS2c)


model.sel(fit_binomS2b,fit_binomS2c, fit_binomS2d, fit_binomS2e, fit_binomS2f) 

###
summary(fit_binomS2c)
#####
##Experimental regeneration
#######
rm(list=ls())

library(nlme)
library(ggplot2)
library(dplyr)
library(MuMIn)
library(glmmTMB)
setwd("C:\\Users\\xjw5\\Dropbox\\NAU\\Data\\Denali\\R code")

site<-read.csv("Denali_sitedata_04082018.csv")
names(site)

exp<-read.csv("Exp_Regen_Data.csv")
names(exp)
head(exp)

exp=mutate(exp, seed.ratio=250/X2018_regen)
head(exp)

exp
exp[,22][is.infinite(exp[,22])] = NA
mean(exp$seed.ratio, na.rm=TRUE)
min(exp$seed.ratio, na.rm=TRUE)
max(exp$seed.ratio, na.rm=TRUE)


exp=mutate(exp, seed.ratio1=250/exp.regen)
head(exp)

exp
exp[,23][is.infinite(exp[,23])] = NA
mean(exp$seed.ratio1, na.rm=TRUE)
min(exp$seed.ratio1, na.rm=TRUE)
max(exp$seed.ratio1, na.rm=TRUE)


summary<-exp%>%
  group_by(seed, scar)%>%
  summarize(sum.exp=sum(exp.regen), 
            sum.exp18=sum(X2018_regen))

summary1=mutate(summary, mean=(sum.exp+sum.exp18)/2)

names(exp)

summary<-exp%>%
  summarize(sum.exp=sum(X2018_regen))


exp$burn.year<-factor(exp$burn.year, levels=c("control", "burned"))

p<-ggplot(exp, aes(x =seed, y = exp.regen, fill=scar)) + geom_boxplot() +labs(x="", y="Experimental regeneration (#)") 
p<- p+ theme_bw() +scale_fill_manual(values=c("goldenrod1", "dodgerblue3"))+facet_wrap(~burn.year)
p +theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"))

418+235
567+305
653/30000
653/872
418/567
235/305

p<-ggplot(exp, aes(x =seed, y = X2018_regen, fill=scar)) + geom_boxplot() +labs(x="", y="Experimental regeneration (#)") 
p<- p+ theme_bw() +scale_fill_manual(values=c("goldenrod1", "dodgerblue3"))+facet_wrap(~burn.year)
p +theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"))


p<-ggplot(exp, aes(y =exp.regen, x = X2018_regen)) + geom_point(aes(fill=scar), pch=21, size=2) +labs(x="2018 (#)", y="2017 (#)") 
p<- p+ theme_bw()+facet_wrap(~burn.year*seed)  +scale_fill_manual(values=c("goldenrod1", "dodgerblue3"))
p +theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold")) +geom_abline()


head(exp)
exp1=dplyr::select(exp, burn.year, grid, gridpoint, plot, Block, Quad, seed, scar, bd, exp.regen, X2018_regen)

library(reshape)
exp.melt=melt(exp1, id=c("burn.year", "gridpoint", "grid", "plot","Block","Quad","seed","scar","bd"))
head(exp.melt)
exp.melt=dplyr::rename(exp.melt, year=variable)
exp.melt=dplyr::rename(exp.melt, regen=value)


y<-lme(data=exp.melt, regen~seed*burn.year*scar, weights=varIdent(form=~seed), random=~1|grid/gridpoint/year, na.action=na.omit, method="ML")
summary(y)
plot(y)
y1 <- update(y, .~. -seed:burn.year:scar)
plot(y)
anova(y, y1)

plot(x=exp.melt$seed, y=resid(y, type="pearson"), ylab="Pearson Residuals", main="Incubation Time", xlab="Density")
abline(h = 0, lty = 2)
plot(x=exp.melt$burn.year, y=resid(y), ylab="Pearson Residuals", main="Incubation Time", xlab="Density")
abline(h = 0, lty = 2)
plot(x=exp.melt$scar, y=resid(y), ylab="Pearson Residuals", main="Incubation Time", xlab="Density")
abline(h = 0, lty = 2)

head(exp.melt)
p<-ggplot(exp.melt, aes(x =seed, y = regen, fill=scar)) + geom_boxplot() +labs(x="", y="Experimental regeneration (#)") 
p<- p+ theme_bw() +scale_fill_manual(values=c("goldenrod1", "dodgerblue3"))+facet_grid(burn.year~year)
p +theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"))


## nibinom2 w/zero inflation
fit_binomS2=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                    data=exp.melt, 
                    ziformula = ~1,  # allows structural zeros to depend on all variables
                    family=nbinom2) ## uses long link
#### so results are on natural log scale##


#### nbinom 2 # no zero infl.
fit_binomS2b=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                     data=exp.melt,
                     ziformula = ~0,  # 
                     family=nbinom2) ## uses log link

### poisson w/zero inflations
fit_poiS1=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                  data=exp.melt,
                  ziformula = ~1, 
                  family=poisson)## 

### poisson # no zero infl
fit_poiS1b=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                   data=exp.melt, 
                   ziformula = ~0, 
                   family=poisson)## 
summary(fit_poiS1b)
#### nbinom1 w/ zero inflation
fit_bS1=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                data=exp.melt,
                ziformula = ~1, 
                family=nbinom1)## 
#### nbinom1 - no zero inf.
fit_bS1b=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                 data=exp.melt, 
                 ziformula = ~0, 
                 family=nbinom1)



model.sel(fit_binomS2,fit_binomS2b,fit_bS1,fit_bS1b,fit_poiS1,fit_poiS1b) 



#### nbinom1 - no zero inf.
fit_bS1b=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                 data=exp.melt, 
                 ziformula = ~1, 
                 family=nbinom1)
summary(fit_bS1)
plot(simulateResiduals(fit_bS1b))


fit_bS1b=glmmTMB(regen~seed*burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
                 data=exp.melt, 
                 ziformula = ~0, 
                 family=nbinom1)
summary(fit_bS1b)
AIC(fit_bS1b)

x=glmmTMB(regen~seed*burn.year+scar+(1|grid/gridpoint/Block) +(1|year),
          data=exp.melt, 
          ziformula = ~0, 
          family=nbinom1)
summary(fit_bS1b)

y=glmmTMB(regen~seed+burn.year*scar+(1|grid/gridpoint/Block) +(1|year),
          data=exp.melt, 
          ziformula = ~0, 
          family=nbinom1)

z=glmmTMB(regen~seed*scar+burn.year+(1|grid/gridpoint/Block) +(1|year),
          data=exp.melt, 
          ziformula = ~0, 
          family=nbinom1)



w=glmmTMB(regen~seed+burn.year+scar+(1|grid/gridpoint/Block) +(1|year),
          data=exp.melt, 
          ziformula = ~0, 
          family=nbinom1)
summary(fit_bS1b)

model.sel(x,y,z,w,fit_bS1b )


x=glmmTMB(regen~seed*burn.year+scar+(1|grid/gridpoint/Block) +(1|year),
          data=exp.melt, 
          ziformula = ~0, 
          family=nbinom1)
summary(x)
AIC(x)
r.squaredGLMM(x)

x1=glmmTMB(regen~seed*burn.year+(1|grid/gridpoint/Block) +(1|year),
           data=exp.melt, 
           ziformula = ~0, 
           family=nbinom1)
AIC(x1)

anova(x, x1)

x2=glmmTMB(regen~seed+burn.year+(1|grid/gridpoint/Block) +(1|year),
           data=exp.melt, 
           ziformula = ~0, 
           family=nbinom1)
AIC(x2)
anova(x, x2)

library(DHARMa)
plot(simulateResiduals(x))


library(effects)
plot(allEffects(x))
summary(x)
r.squaredGLMM(x)
