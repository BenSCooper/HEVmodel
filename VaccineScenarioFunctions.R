  # single patch Hep E model with vaccination, allowing for people to move into and out of patch at a specified rate
# SEIR model structure, with probabilities for imports to be in each state (S, E , I, R)
# Also allows for vaccination - where we have one or two doses (with different effectiveness)
# The two doses are typically delayed by one 1 month and there is a delay from vaccination to immunity.
# In contrast to " model 1 - evaluate vaccine.r" where  vaccination is typically carried out at certain time periods over a 
# specified time period (which under some circumstances caused numerical problems due ot widely varying rates) in this version 
# We assume immusation happens instantaneously at a specified time (so this time should account for the delay from vaccination
# till immunity) 
# For the moment we assume the each dose of vaccine has some probablity of generating immunity. 
# Later we could also consider model structures where each dose generates partial immunity (see recent papers from Gabriela Gomes to see the difference)

library(deSolve)



get.prob.death.if.symptomatic<-function(assumption=1, pregnant=FALSE,age.grp=NA){
  # function returns probablity of death from HEV given symptomatic infection  
  # assumption is the assumption number. Currently implmented values are
  # assumption =1 : probabilities derived from posterior distributions reported by Rein RB et al Hepatology 2012 paper   
  # assumption =2 : probabilities derived  from results reported in Teshale EH et al, EID 2012 DOI: 10.3201/eid1601.090764 
  # age.grp is only used if assumption ==2 
  
  if(assumption==1){ # use numbers from Rein et al 2012
    if(!pregnant){
      prob<-rbeta(1, 340.7925, 17595.6551) # draw a random number from beta distribution with same mean and 95% CrI as reported in Rein et al: 0.019 (0.017, 0.021)
    } else { #  i.e pregnant
      prob<-rbeta(1, 143.6822, 581.9855)  # 95% CrI as reported in Rein et al: 0.198 (0.17, 0.23)
    }
  } else if (assumption ==2){  # From Teshale (allow probabilities to depend on age and pregnant status)
    stop("Assumption 2 not yet implemented")
  } else {
    stop("Invalid assumption option")
  }
  return(prob)
}

get.prob.symptomatic.if.infected<-function(assumption=1, pregnant=FALSE,age.grp=NA){
  # function returns probablity of symptoms of HEV infection given  infection  
  # assumption is the assumption number. Currently implmented values are
  # assumption =1 : probabilities derived from posterior distributions reported by Rein RB et al Hepatology 2012 paper   
  # assumption =2 : probabilities derived  from results reported in Teshale EH et al, EID 2012 DOI: 10.3201/eid1601.090764 
  # age.grp is only used if assumption ==2 
  
  if(assumption==1){ # use numbers from Rein et al 2012 (same for pregnant and non-pregant)
    prob<-rbeta(1, 133.7627, 541.8066) # draw a random number from beta distribution with same mean and 95% CrI as reported in Rein et al: 0.198 (0.17, 0.23)
  } else if (assumption ==2){ # From Teshale (allow probabilities to depend on age and pregnant status)
    stop("Assumption 2 not yet implemented")
  } else {
    stop("Invalid assumption option")
  }
  return(prob)
}  
  






seasonal.transmission.parameter<-function(beta.max, beta.min, phi, t){
  # returns seasonal transmission parameter at time t, assumed to have a max value of beta.max, a min value of beta.min and phase of phi
  # patterns of seasonality can be specified with arbitrary function here   
  # t is assumed to be specified in units of one year   
  #return(beta.max) # use this line for no seasonal variation 
  return((beta.min + beta.max)/2  + (beta.max-beta.min)*sin(2*pi*t-phi)/2) # sine wave varying between beta.min and beta.max with phase of phi
}



immunisation.event.forv3<-function(t, y,parms){
  # update state variable for HepEwithVAxmod model
  #  Here vacc.eff1 is prob of immunity from vaccine after first dose. 
  #  and vacc.eff2 is prob of immunity after having received two dose
  #  Assume that vaccination has no effected on someone on the E, I or R classes so don't need to split these into vaccinated and not
  with(as.list(c(y,parms)),{
   vacc.efficacy.given.prior.vacc<-max(0,(vacc.eff2-vacc.eff1)/(1-vacc.eff1) )# i.e prob of getting immunity given susceptible and previously vaccinated
   VS <- VS*(1- proportion.to.vaccinate* vacc.efficacy.given.prior.vacc) + SV * (1- vacc.eff1)* proportion.to.vaccinate         # susceptibles in vaccinatable population , who have been vaccinated but are still susceptible
   RV<-RV+ vacc.eff1*SV*proportion.to.vaccinate + vacc.efficacy.given.prior.vacc*proportion.to.vaccinate*VS

   SV<-SV * (1- proportion.to.vaccinate) # susceptibles in vaccinatable population , not yet vaccinated 
   V2<- V2 + (V1-V2)*proportion.to.vaccinate # number having recieved at least two doses  
   V1<- V1+ (SV +EV +IV +RV+ -V1)* proportion.to.vaccinate # number having receIved at least one dose = no. vaccinatable who have not been vaccinated * proportion of those to vaccinate + number who had previously received at least one dose. 
    return(c(SNV,SV, ENV,EV,INV, IV,RNV, RV, VS, V1, V2,CumIncNV, CumIncV ))
 })
}


HepEwithVAxmod.v3 <- function(t, x, parms) { 
  # This is the standard SEIR model. 
  #v3  allows for a different efficacy after second dose of vaccine - but in this case 
  # vaccination is implemented as an event in the above function immunisation.event (avoids numerical problems when vaccination is a rate)
  # collects output data on the number of infections in the vaccinated and unvaccinated populations  
  # Also to allow for the fact that a proportion of the population are not elligible for vaccination we divde susceptibles into 
  # 2 groups SV (susceptible vaccine candidates) and SNV susceptible non vaccine candidates
  # we also assume that the time it takes to generate immunity is fixed (rather than exponentially distributed) so 
  # delay from the vaccination time to immunity is dealt with when specify the time of the immunisation (should account for delay from vaccination)
  
  
  # function returning a list holding the model equations and the rates of import and export of people into the opulation at diferent time points. 
  # t is time, 
  #x is the vector of state variables (SNV,SV VS, VE,  EV, ENV, IV, INV, RV, RNV) # where SV is susceptible and unvaccinated  but vaccinatable, SNV is  Susceptible unvaccinated and unvaccinatable, VS is unusccesfully vaccinated (and still scuscptible) 
  # and VE is successfully vaccinated and waiting for immunity, EV and ENV are respectively latently infected and vaccinatable and latently infectied and non-vaccinatable. Similarly for RV and RNV (both infec)
  # we also track i) the number of people vaccinated with the first dose of vaccine V1
  # the number of poeple vaccinated with the second dose, V2, and cumulative number of HepE in fections, CumInc
  
  
  #  (* Note that if we have an imperfect vaccine  a real vaccination program wouldn't vaccinate people more than once in a each vaccination episode 
  # (where an episode is a  short period lasting a few days when a camp is vaccinated) - so to allow for this fact if vaccinated unsuccesfully 
  # vaccinated people move to a compartment VS (vaccinated and still susceptible). If susccessfully vaccinated they move to VE where they may still be infected. 
  # but after a while progress to R from VE as immunity developes. 
  # We assme that patients in the VE compartment can still be infected, and that vaccine has no effect on anyone who is in the E, I or R state when vaccinated. 
  # We want to model the use of up to two doses of vaccine (typically about a month apart). When doing this we assume that for both the first and second vaccination
  # episodes everyone in the target age group is equally likely to be vaccinated (and therefore that someone missing the first vaccine is no more or less likely to miss the
  #next vaccine). Those who are vaccinated effectiveley but not yet immune (VE) can still be infected.
  # Finally, in this version of the model we assume that each dose of vaccine to someone susceptible (S or VS) either gives full protection (with some fixed probabaility that is 
  # the same for the first or second dose) or no protection.
  
  
  # parms is the vector of parameters
  # For vaccination, we assume that vaccine only has any effect if used in those who are in S state of VS (vaccinated but still suceptible due to vaccine failure)
  
  with(as.list(c(parms, x)), {
    # beta.max is max seasonal tranmission parameter
    # beta.min is min seasonal transmission paramter,
    # phi is phase angle which determines what time of year peak(s) occurs
    # N is current population size
    # gamma is rate of prgression to exposed to infectious
    # rho is rate of progression from infectious to recovered.
    # vacc.eff1 is vaccine efficacy of first dose
    # vacc.eff2 is vaccine efficacy of second dose
    #  proportion.vaccinatable is the proportion of the population in the vaccine candidate population (e.g. proportion over 15 and not pregnant)
    #  proportion.to.vaccinate is the proportion of the vaccine candidate population to be vaccinated at each vaccination episode
    # p.imp.S  is proportion of imported cases are are susceptible
    # p.imp.E  is proportion of imported cases are are latently infected
    # p.imp.I  is proportion of imported cases are are infectious
    # p.imp.R  is proportion of imported cases are are immune
    
    # (note that p.imp.S, p.imp.E, p.imp.I, p.imp.R should sum to one)
    
    N<-SNV + SV +ENV + EV+INV +IV +RNV +RV +VS    
    import<-movements.into.patch(t)    # imports to patch at a given time 
    export<-movements.out.of.patch(t)    # exports  out of patch at a given time
        
    beta<-seasonal.transmission.parameter(beta.max,beta.min,phi,t)

    #  print(c("vac1and2",vac1,vac2))
    # vacc.eff2given1 is probability of immunisation from vaccine 2 given lack of immunisaiton from 1. 
    # derived from the probablity of successfull immunisation after 1 or 2 doses vacc.eff1 and  vacc.eff2
    # using (1-vacc.eff2)=(1-vacc.eff1)(1-vacc.eff2given1)
    
    # Derivatives 
    
    dSNV <- -beta*SNV*(INV +IV)/N    + import*p.imp.S* (1-proportion.vaccinatable) - export*SNV/N
    dSV <- -beta*SV*(INV +IV)/N    + import*p.imp.S* proportion.vaccinatable - export*SV/N
    dENV <-  beta*SNV*(INV +IV)/N - gamma*ENV  + import*p.imp.E - export*ENV/N
    dEV <-  beta*(VS +SV)*(INV +IV)/N - gamma*EV   - export*EV/N  #vaccinatable group (no necessarily vaccinated)
    dINV <- gamma*(ENV) - rho*INV  + import*p.imp.I* (1-proportion.vaccinatable) - export*INV/N
    dIV <- gamma*(EV) - rho*IV  + import*p.imp.I* proportion.vaccinatable - export*IV/N
    dRNV <- rho*INV +  import*p.imp.R*(1-proportion.vaccinatable) - export*RNV/N
    dRV <- rho*IV +  import*p.imp.R* proportion.vaccinatable - export*RV/N
    dVS<- - beta*VS*(INV +IV)/N  - export*VS/N
    dV1<-  - export*V1/N
    dV2<-  - export*V2/N
    dCumIncNV<-gamma*ENV  #cumulative incidence of infection in not vaccinatable group
    dCumIncV<-gamma*EV    #cumulative incidence of infection in vaccinatable group
    
    # Return values 
    list(c(dSNV,dSV,dENV,dEV,dINV, dIV,dRNV, dRV, dVS, dV1, dV2,dCumIncNV, dCumIncV ))
    
  })  # end with
} # end HepEwithVAxmod.v3



HepEwithVAxmod.SEEIR <- function(t, x, parms) { 
  # This is the  SEEIR model. 
  # allows for a different efficacy after second dose of vaccine - but in this case 
  # vaccination is implemented as an event in the above function immunisation.event (avoids numerical problems when vaccination is a rate)
  # collects output data on the number of infections in the vaccinated and unvaccinated populations  
  # Also to allow for the fact that a proportion of the population are not elligible for vaccination we divde susceptibles into 
  # 2 groups SV (susceptible vaccine candidates) and SNV susceptible non vaccine candidates
  # we also assume that the time it takes to generate immunity is fixed (rather than exponentially distributed) so 
  # delay from the vaccination time to immunity is dealt with when we specify the time of the immunisation (should account for delay from vaccination)
  
  
  # function returning a list holding the model equations and the rates of import and export of people into the opulation at diferent time points. 
  # t is time, 
  #x is the vector of state variables (SNV,SV VS, VE,  EV, ENV, IV, INV, RV, RNV) # where SV is susceptible and unvaccinated  but vaccinatable, SNV is  Susceptible unvaccinated and unvaccinatable, VS is unusccesfully vaccinated (and still scuscptible) 
  # and VE is successfully vaccinated and waiting for immunity, E1V and E1NV are respectively latently infected and vaccinatable and latently infectied and non-vaccinatable (similarly for E2V and E2NV).
  # Similarly for RV and RNV (both infec)
  # we also track i) the number of people vaccinated with the first dose of vaccine V1
  # the number of poeple vaccinated with the second dose, V2, and cumulative number of HepE in fections, CumInc
  
  
  #  (* Note that if we have an imperfect vaccine  a real vaccination program wouldn't vaccinate people more than once in a each vaccination episode 
  # (where an episode is a  short period lasting a few days when a camp is vaccinated) - so to allow for this fact if vaccinated unsuccesfully 
  # vaccinated people move to a compartment VS (vaccinated and still susceptible). If susccessfully vaccinated they move to VE where they may still be infected. 
  # but after a while progress to R from VE as immunity developes. 
  # We assme that patients in the VE compartment can still be infected, and that vaccine has no effect on anyone who is in the E, I or R state when vaccinated. 
  # We want to model the use of up to two doses of vaccine (typically about a month apart). When doing this we assume that for both the first and second vaccination
  # episodes everyone in the target age group is equally likely to be vaccinated (and therefore that someone missing the first vaccine is no more or less likely to miss the
  #next vaccine). Those who are vaccinated effectiveley but not yet immune (VE) can still be infected.
  # Finally, in this version of the model we assume that each dose of vaccine to someone susceptible (S or VS) either gives full protection (with some fixed probabaility that is 
  # the same for the first or second dose) or no protection.
  
  
  # parms is the vector of parameters
  # For vaccination, we assume that vaccine only has any effect if used in those who are in S state of VS (vaccinated but still suceptible due to vaccine failure)
  
  with(as.list(c(parms, x)), {
    # beta.max is max seasonal tranmission parameter
    # beta.min is min seasonal transmission paramter,
    # phi is phase angle which determines what time of year peak(s) occurs
    # N is current population size
    # gamma is rate of prgression to exposed to infectious
    # rho is rate of progression from infectious to recovered.
    # vacc.eff1 is vaccine efficacy of first dose
    # vacc.eff2 is vaccine efficacy of second dose
    #  proportion.vaccinatable is the proportion of the population in the vaccine candidate population (e.g. proportion over 15 and not pregnant)
    #  proportion.to.vaccinate is the proportion of the vaccine candidate population to be vaccinated at each vaccination episode
    # p.imp.S  is proportion of imported cases are are susceptible
    # p.imp.E  is proportion of imported cases are are latently infected
    # p.imp.I  is proportion of imported cases are are infectious
    # p.imp.R  is proportion of imported cases are are immune
    
    # (note that p.imp.S, p.imp.E, p.imp.I, p.imp.R should sum to one)
    
    N<-SNV + SV +E1NV + E1V + E2NV + E2V + INV +IV +RNV +RV +VS    
    import<-movements.into.patch(t)    # imports to patch at a given time 
    export<-movements.out.of.patch(t)    # exports  out of patch at a given time
    
    beta<-seasonal.transmission.parameter(beta.max,beta.min,phi,t)
    
    #  print(c("vac1and2",vac1,vac2))
    # vacc.eff2given1 is probability of immunisation from vaccine 2 given lack of immunisaiton from 1. 
    # derived from the probablity of successfull immunisation after 1 or 2 doses vacc.eff1 and  vacc.eff2
    # using (1-vacc.eff2)=(1-vacc.eff1)(1-vacc.eff2given1)
    
    # Derivatives 
    
    dSNV <- -beta*SNV*(INV +IV)/N    + import*p.imp.S* (1-proportion.vaccinatable) - export*SNV/N
    dSV <- -beta*SV*(INV +IV)/N    + import*p.imp.S* proportion.vaccinatable - export*SV/N
    dE1NV <-  beta*SNV*(INV +IV)/N - gamma*E1NV  + import*p.imp.E - export*E1NV/N
    dE1V <-  beta*(VS +SV)*(INV +IV)/N - gamma*E1V   - export*E1V/N  #vaccinatable group (no necessarily vaccinated)
    dE2NV <-  gamma*E1NV  - gamma*E2NV - export*E2NV/N
    dE2V <-  gamma*E1V - gamma*E2V   - export*E2V/N  #vaccinatable group (no necessarily vaccinated)    
    dINV <- gamma*(E2NV) - rho*INV  + import*p.imp.I* (1-proportion.vaccinatable) - export*INV/N
    dIV <- gamma*(E2V) - rho*IV  + import*p.imp.I* proportion.vaccinatable - export*IV/N
    dRNV <- rho*INV +  import*p.imp.R*(1-proportion.vaccinatable) - export*RNV/N
    dRV <- rho*IV +  import*p.imp.R* proportion.vaccinatable - export*RV/N
    dVS<- - beta*VS*(INV +IV)/N  - export*VS/N
    dV1<-  - export*V1/N
    dV2<-  - export*V2/N
    dCumIncNV<-gamma*E2NV  #cumulative incidence of infection in not vaccinatable group
    dCumIncV<-gamma*E2V    #cumulative incidence of infection in vaccinatable group
    
    # Return values 
    list(c(dSNV,dSV,dE1NV,dE1V, dE2NV,dE2V,dINV, dIV,dRNV, dRV, dVS, dV1, dV2,dCumIncNV, dCumIncV ))
    
  })  # end with
} # end HepEwithVAxmod.SEEIR



# time vector ..specifies times at which to evaluate model (units are one year)
times <- seq(0, 10, by = 0.01)

# ****************************************************
#  Specify movements in and out of the patch here

# for now just use arbitrary function for importation rate at specific times and interpolate this to get imports to patch at different times
# note that all rates are expressed in time units of 1 year


import.rates <- as.data.frame(list(times = times, import = rep(0, length(times)))) # set import element to values >0 to allow imports

movements.into.patch<-approxfun(import.rates, rule=2) # import rates specified by interpolation. Here rule 2 specifies that rates outside interval are taken from the closest extreme in the interval

export.rates <- as.data.frame(list(times = times, export = rep(0, length(times)))) # set export element to values >0 to allow exports from patch

movements.out.of.patch<-approxfun(export.rates, rule=2) # import rates specified by interpolation. Here rule 2 specifies that rates outside interval are taken from the closest extreme in the interval
# ****************************************************



# Specify initial values for the model compartments
N<-100000
S0<-N-1   #initial number susceptible etc
E0<-1
EV0<-1
ENV0<-0
IV0<-0
INV0<-0
RV0<-0
RNV0<-0
VS0<-0
VE0<-0
V10<-0
V20<-0
CumInc0<-0


# Specify parameter values 

R0max<-2 #maximum seasonal R0
R0min<-2 #minimum seasonal R0
latent.period<-0.04  # in units of 1 year
infectious.period<-0.1  # in units of 1 year
mean.time.from.vacc.to.immunity<-0.04   # in units of 1 year


parameters<-c(
  beta.min=R0min/infectious.period,
  beta.max=R0max/infectious.period,
  phi=0,  # phase angle of seasonal forcing in radians (specifies when beak occurs)
  gamma=1/latent.period,
  rho=1/infectious.period,
  p.imp.S=1,  # proportion imported cases who are susceptible etc
  p.imp.E=0,
  p.imp.I=0,
  p.imp.R=0,
  vacc.eff1=0.8,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
  vacc.eff2=0.8,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
  proportion.to.vaccinate=.9
  )

# set starting values for state variables
proportion.vaccinatable<-0.8
I0<-0; VS0 <-0 ; V1<-0; V20<-0
proportion.initially.immune<-0

xstart<-c(
  SNV=S0*(1-proportion.vaccinatable),
  SV=S0*proportion.vaccinatable,
  ENV=1,
  EV=0,
  INV=I0*(1-proportion.vaccinatable),
  IV=I0* proportion.vaccinatable,
  RNV= proportion.initially.immune*(1-proportion.vaccinatable),
  RV= proportion.initially.immune*proportion.vaccinatable,
  VS=VS0,
  V1=V10,
  V2=V20,
  CumIncNV=0,
  CumIncV=0
  )

# # Now call ode to solve the equations and store in out \
# params<-as.list(parameters)
# 
# # for some resason parameters don't seem to be passed to event function (check documentation to see if they should)
# vacc.eff1<-params$vacc.eff1 
# vacc.eff2<-params$vacc.eff2
# 
# out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,params)
# out2<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,params,events=list(func=immunisation.event.forv3, time=c(0.5,1.5)))
# 
# plot(out) #display output of a single run

# #Now with sliders
# 
# require(manipulate)
# 
# # define function required for plotting with sliders
# plot.modelv3<-function(Slider.R0max,  Slider.R0min, Slider.vacc.eff, Slider.vacc.coverage,Slider.vax.time1, Slider.vax.time2,delay.to.immunity){
#   parms<-parameters # take default parameters for anything not specified by  alisder
#   infectious.period<-1/as.list(parms)[[match("rho",names(parms))]]
#   beta.max<-Slider.R0max/infectious.period
#   beta.min<-Slider.R0min/infectious.period
#   pos<-match("beta.max",names(parms)) # position in parameter vector of beta.max
#   parms[pos]<-beta.max
#   pos<-match("beta.min",names(parms)) # position in parameter vecgtor of beta.min
#   parms[pos]<-beta.min
#   pos<-match("vacc.eff",names(parms)) 
#   parms[pos]<-Slider.vacc.eff
#   pos<-match("proportion.to.vaccinate",names(parms)) 
#   parms[pos]<-Slider.vacc.coverage  
#    
#   vacc.t1<-Slider.vax.time1
#   vacc.t2<-Slider.vax.time2
#   out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,parms,events=list(func=immunisation.event.forv3,time=c(vacc.t1+delay.to.immunity,vacc.t2+delay.to.immunity)))
#   time<-out[,1]
#   I<-out[,6] # currently infected
#   Cum<-out[,11] + out[,12] #cumulative cases
#   N<-max(Cum)
#   plot(time, I, type='l', xlab="Year",ylab="Prevalence",col="red",main="",ylim=c(0,N))
#   lines(time, Cum,col="plum")
#   legend("topright",c("I","Cumulative cases"),lty=1,col=c("red","plum") )   
# }
# 
# 
# 
# dev.off()
# manipulate(
#   plot.modelv3(R0max,  R0min, vacc.eff, vacc.coverage,Slider.vax.time1,Slider.vax.time2,delay.to.immunity),  
#   R0min=slider(0,5,2,step=.01),
#   R0max=slider(1,5,2,step=.01),
#   vacc.eff=slider(0,0.99,0.8,step=.01),
#   vacc.coverage=slider(0,0.99,0.0,step=.01),
#   Slider.vax.time1=slider(0,5,2,step=.01),
#   Slider.vax.time2=slider(0,5,2.01,step=.01),
#   delay.to.immunity=slider(0,4/52,0,step=.01)
#   )  



# Calculate infections, cases and 


# Estimate 2 dose vaccine efficacy

P<-100
n.samples<-length(posterior.samples[,1])
runs.to.plot<-  c(1:n.samples)[order(runif(n.samples))][1:P]



# now add a line for what Zhu call the first two doses subset (though not totally sure how this is defined)
n_v<- 54986  # number vaccinated  in "first doses subset"
n_p<- 54973  # number getting placebo  in "first doses subset"
c_v<- 0  # HEV cases in those vaccinated in "first doses subset"
c_p<- 5  # HEV cases in those getting placebo  in "first doses subset"

ARV<-c_v/n_v # attack rate vaccinated
ARU<-c_p/n_p # attack rate unvaccinated

VE<-100*(ARU-ARV)/ARU
VE1.post.sample<-rep(0,n.samples)  #  default assumption of no efficacy of a single dose of vaccine

p_v.post.sample<-rbeta(n.samples, 1+ c_v, 1+ n_v-c_v)
p_p.post.sample<-rbeta(n.samples,1+ c_p, 1+ n_p-c_p)

VE2subset.post.sample<-100*(p_p.post.sample -p_v.post.sample)/p_p.post.sample

N <-100 # number of samples

# scenario 1 - no vaccine, no watsan 
# Specify initial values for the model compartments
N<-10000
S0<-N-1   #initial number susceptible etc
E0<-1
EV0<-1
ENV0<-0
I0<-0
R0<-0
VS0<-0
VE0<-0
V10<-0
V20<-0
CumInc0<-0

proportion.pregnant<-.03 
relative.risk.of.symptoms.given.infected.if.pregnant<-2   # can we find data for this??? This just a guess based on Khuroo 1980
xstart<-c(
  S=S0,
  E=E0,
  I=I0,
  R=R0,
  VS=VS0,
  VE=VE0,
  V1=V10,
  V2=V20,
  CumInc=CumInc0
  )



results.scenario1<-data.frame(infections=0, cases=0, deaths.preg=0, deaths.nonpreg=0 )


mean.time.from.vacc.to.immunity <-2/52  # assume on average it takes two weeks to develop immunity

# 
# for(j in 1:length(runs.to.plot)){
#   i<-runs.to.plot[j]
#   params<-modv1.360to420k.posterior.sample[i,]  # sample parameters from posterios
#   VE<-VE2subset.post.sample[i]
#   prob.death.pregnant<-rbeta(1, 1+54, 132-54)
#   
#  prob.case.if.preg<-rbeta(10000,5+1, 1+359-5)   #  these numbers based on Jamam data 
#  prob.case.if.not.preg<-rbeta(10000,92+1, 1+23934-92)
#  relative.risk.of.symptoms.given.infected.if.pregnant <-prob.case.if.preg/prob.case.if.not.preg
# 
#   prob.death.not.pregnant<-0.02  # really we need a distribution for this. 
#   parameters<-c(
#     beta.min=params$R01* params$rho ,
#     beta.max=params$R01* params$rho,
#     phi=0,  # phase angle of seasonal forcing in radians (specifies when beak occurs)
#     gamma=params$gamma,
#     rho=params$rho,
#     p.imp.S=1,  # proportion imported cases who are susceptible etc
#     p.imp.E=0,
#     p.imp.I=0,
#     p.imp.R=0,
#     vacc.eff=VE/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
#     vac.immun.rate=1/mean.time.from.vacc.to.immunity,
#     proportion.to.vaccinate=0.6,
#     start.time1=Inf, # first vaccine dosing start time 
#     stop.time1=Inf,   # first vaccine dosing stop time
#     start.time2=Inf, # second vaccine dosing start time
#     stop.time2=Inf   # second vaccine dosing stop time
#     )
#   out<-ode(y=xstart,times=times, func=HepEwithVAxmod,parameters)
#   total.infected<-max(out[,10])
#   total.symptomatic<-total.infected*params$proportion.symptomatic
#   # assuming same risk of infection if pregant as others but a different risk of symptoms so... 
#   prob.symptoms.if.not.pregnant<- params$proportion.symptomatic/(1+proportion.pregnant*(relative.risk.of.symptoms.given.infected.if.pregnant-1))
#   prob.symptoms.if.pregnant<- relative.risk.of.symptoms.given.infected.if.pregnant*prob.symptoms.if.not.pregnant
#   Number.cases.pregnant<-total.infected*proportion.pregnant*prob.symptoms.if.pregnant
#   Number.cases.not.pregnant<-total.infected*(1-proportion.pregnant)*prob.symptoms.if.not.pregnant
#   Number.deaths.pregnant<-prob.death.pregnant*Number.cases.pregnant 
#   Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.cases.not.pregnant 
#   
#   results.scenario1[j,]<-c(infections=total.infected, cases=total.symptomatic, deaths.preg=Number.deaths.pregnant, deaths.nonpreg=Number.deaths.not.pregnant)
#   
#   
# }


run.vaccine.scenario<-function(num.its=1,  times=seq(0,3,1/52), cases.b4.vacc=NA, start.V1=Inf, stop.V1=Inf, start.V2=Inf,stop.V2=Inf,vacc.eff1=NA,vacc.eff2=NA,immun.delay1=0, immun.delay2=0, proportion.vaccinatable=0.8,proportion.of.these.to.vax=1,pregnant.vaccinated=FALSE, proportion.pregnant=0.04){
  
  # WARNING..solution to ODE sysmtem seems to be with unstable  when vaccine rates are very high....i.e. if vaccination takes place in a very short time period so either use this version wiht extra care or use run.vaccine.scenario.v2 instead
  # where within that time period (before start and stop of vaccination) vaccination occurs.
  # The quick fix in this code is to i) insist start and top times differ by at least 0.02 years (about a week)
  # Two possible better fixes when time are 
  #i) use ODE commands to move people between compartments instantaneously  (this is now implemented in run.vaccine.scenario.v2)
  #ii) adjust settings of solver (or use another deSolve solver) to cope with stiff problems. See manual. Currently using 
  # lsodes which is designed for stiff systems...see  a little better than lsoda
  # For now that just leave it so that delivering a dose takes 0.02 of a year , but need to use this function wth 
  # care, particularly when high coverage levels are required. 
  #
  #
  # cases.b4.vacc is equal to the number of cases before vaccine starts. After these cases
  # assume dose 1 is given after a delay start.v2 and second dose after a delay start.V2. If cases.b4.vacc is NA there is no reactive vaccination
  #   num.its is the number of iterations
  #	start.V1 etc give start and stop times for delivering first and second dase (time units of one year).
  #   vacc.eff1 and vacc.eff2 are percentage who are successfully immunized after 1 or 2 doses (between 0 and 100)
  #  immun.delay1 and immun.delay2 and delays from vaccination to immunity being established
  #	proportion.vaccinatable is proportion of the population who ar eligible to receive the vaccine
  #	proportion.of.these.to.vax is proportion of the population eligible to actually be vaccinated
  #  note that if pregnant.vaccinated is true then proportion.pregnant must be <= proportion.vaccinatable 
  # if pregnant.vaccinated is false then proportion.pregnant must be <= (1-proportion.vaccinatable )
  
  #   note that these numbers don't account for binomial variation in number of deaths etc ...so we are looking at expected values 
  # assumes posterior.samples,  VE1.post.sample and VE2subset.post.sample  is globally defined (latter 2 are vaccine efficacies after 1 or 2 doses)
  n.samples<-length(posterior.samples[,1])
  runs.to.plot<-  c(1:n.samples)[order(runif(n.samples))][1:num.its]
  results<-data.frame(infections=0, cases=0, cases.pregnant =0, cases.not.pregnant=0 ,deaths.preg=0, deaths.nonpreg=0 )
  #if((stop.V1-start.V1)<0.019999 ) stop("Start and stop vaccination times too close together and very high vaccination rates may make estimates unstable. Consider revising code if you need to consider this scenario.")
  #if((stop.V2-start.V2)<0.019999 ) stop("Start and stop vaccination times too close together and very high vaccination rates may make estimates unstable. Consider revising code if you need to consider this scenario.")
  for(j in 1:length(runs.to.plot)){
    i<-runs.to.plot[j]
    if(is.na(vacc.eff1))   VE1<-VE1.post.sample[i] else VE1<-vacc.eff1
    if(is.na(vacc.eff2))   VE2<-VE2subset.post.sample[i] else VE2<-vacc.eff2
    
    prob.symptomatic.infection.pregnant <-get.prob.symptomatic.if.infected(assumption = 1, pregnant = TRUE)
    prob.symptomatic.infection.not.pregnant <-prob.symptomatic.infection.pregnant  # baseline assumption is no difference (lacking data to definitively say)
    
    prob.death.not.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=FALSE) # given symptomatic infection
    prob.death.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=TRUE) # given symptomatic infection
    
    
    params<-posterior.samples[i,]  # 
    parameters<-c(
      beta.min=params$R02* params$rho ,
      beta.max=params$R02* params$rho,
      phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
      gamma=params$gamma,
      rho=params$rho,
      p.imp.S=1,  # proportion imported cases who are susceptible etc
      p.imp.E=0,
      p.imp.I=0,
      p.imp.R=0,
      vacc.eff1=VE1/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
      vacc.eff2=VE2/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
      
      vac.immun.delay1= immun.delay1,
      vac.immun.delay2= immun.delay2,
      proportion.vaccinatable= proportion.vaccinatable,
      proportion.to.vaccinate= proportion.of.these.to.vax,
      start.time1= start.V1, # first vaccine dosing start time 
      stop.time1= stop.V1,   # first vaccine dosing stop time
      start.time2= start.V2, # second vaccine dosing start time
      stop.time2= stop.V2   # second vaccine dosing stop time
      )
    xstart<-c(
      SNV=S0*(1-proportion.vaccinatable),
      SV=S0*proportion.vaccinatable,
      ENV=1,
      EV=0,
      INV=INV0,
      IV=IV0,
      RV=RV0,
      RNV=RNV0,
      VS=VS0,
      V1=V10,
      V2=V20,
      CumInc.NV=0,
      CumInc.V=0
      )
    
    
    if(!is.na(cases.b4.vacc)) {
      #   then  run model and determine when first day when there  are  cases.b4.vacc
      params.novax<-parameters
      params.novax[10] <-0 #params.novax$vacc.eff1<-0
      params.novax[11] <-0 #p#params.novax$vacc.eff2<-0
      #  print(params.novax)
      no.vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v2, rtol = 1e-7, atol = 1e-7,params.novax,method="bdf")
      no.vax.cum.inc<-(no.vax.out[,11]+no.vax.out[,12])*params$proportion.symptomatic
      #print(no.vax.cum.inc)
      n.time.points<-length(no.vax.cum.inc)
      first.time.gt.Xcases<-no.vax.out[min(n.time.points,1+findInterval(c(-Inf, cases.b4.vacc), no.vax.cum.inc)[2]),1]
      print(c("xxx",first.time.gt.Xcases))
      parameters<-as.list(parameters)
      parameters$start.time1<-parameters$start.time1 + first.time.gt.Xcases
      parameters$stop.time1<-parameters$stop.time1 + first.time.gt.Xcases
      parameters$start.time2<-parameters$start.time2 + first.time.gt.Xcases
      parameters$stop.time2<-parameters$stop.time2 + first.time.gt.Xcases
    }
    out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v2,parameters,rtol = 1e-7, atol = 1e-7,method="lsoda")
    # currently use lsodes as fast vax rate can lead to numerical instability and bizarre rates
    # sems to be problem if one dose of vaccine given to most people in less than about 3 days. 
    # ideally should recod to avoid this but to have an option for instantaneous vaccination instead 
    total.infected<-max(out[,11]) + max(out[,12])
    total.infected.vaccinatable<-max(out[,12])
    total.infected.notvaccinatable<-max(out[,11])
    
    total.cases<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    
    # Need to split this output into vaccinatable and not vaccinatable populations so we get total infections and symptomatic in each
    #  using pregnant.vaccinated  proportion.pregnant
    # prob.symptoms.if.not.pregnant<- params$proportion.symptomatic/(1+proportion.pregnant*(relative.risk.of.symptoms.given.infected.if.pregnant-1))
    
    if(pregnant.vaccinated) {
      Number.infected.pregnant<-total.infected.vaccinatable*proportion.pregnant/ proportion.vaccinatable
      prob.pregnant.infected<-total.infected.vaccinatable/(S0*proportion.vaccinatable)
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable-proportion.pregnant
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)
      
    } else {
      Number.infected.pregnant<-total.infected.notvaccinatable*proportion.pregnant/ (1-proportion.vaccinatable)	
      prob.pregnant.infected<-total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
      prob.not.pregnant.infected<-total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)-proportion.pregnant
    }
    Number.infected.not.pregnant<-total.infected-Number.infected.pregnant
    # note that for the output we need to distinguish cases (those which are reported) and symptomatic infections (which may or may not be reported)
   
    relative.risk.of.beingcase.given.infected.if.pregnant <- 2  # assume that those who are pregnant are twice as likely to be reported as cases (irresepective of relative risk of symptoms)  
    
    #note in line below  params$proportion.symptomatic is a bit of misnomer - it really represents the proportion of infections which are reported as cases
    
    prob.case.if.not.pregnant<-params$proportion.symptomatic*(1+ Number.infected.not.pregnant/Number.infected.pregnant) /(relative.risk.of.beingcase.given.infected.if.pregnant + Number.infected.not.pregnant/Number.infected.pregnant)
    prob.case.if.pregnant<- relative.risk.of.beingcase.given.infected.if.pregnant *prob.case.if.not.pregnant
    
    Number.symptomatic.pregnant<-prob.symptomatic.infection.pregnant* Number.infected.pregnant
    Number.symptomatic.not.pregnant<-prob.symptomatic.infection.not.pregnant* Number.infected.not.pregnant
    
    total.symptomatic<-Number.symptomatic.pregnant + Number.symptomatic.not.pregnant
    
    Number.cases.pregnant = prob.case.if.pregnant * Number.infected.pregnant
    Number.cases.not.pregnant<-prob.case.if.not.pregnant * Number.infected.not.pregnant
    Number.deaths.pregnant<-prob.death.pregnant*Number.symptomatic.pregnant 
    Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.symptomatic.not.pregnant 
    
    # print(j)
    results[j,]<-c(infections=total.infected, cases=total.cases, cases.pregnant= Number.cases.pregnant , cases.not.pregnant= Number.cases.not.pregnant, symptomatic=total.symptomatic,  deaths.preg=Number.deaths.pregnant, deaths.nonpreg=Number.deaths.not.pregnant)
    
  } # end for 
  results
}


run.vaccine.scenario.v2<-function(num.its=1,  times=seq(0,3,1/52), cases.b4.vacc=NA,  immunisation.times =Inf,vacc.eff1=NA,vacc.eff2=NA,immun.delay1=1/24, proportion.vaccinatable=0.8,proportion.of.these.to.vax=1,pregnant.vaccinated=FALSE, proportion.pregnant=0.04, campnumber=1){
   # Before running this make sure that the dataframe posterior.samples holds samples from the distribution of parameters values that we want to use inthe simulations.
  # In this version, immunisation takes place instantaneously at immunisation.times
    # (so to allow for a delay from vaccination to immunity should adjust these times)
  # cases.b4.vacc is equal to the number of cases before vaccine starts. 
  # If cases.b4.vacc is NA there is no reactive vaccination
  #   num.its is the number of iterations
  #   vacc.eff1 and vacc.eff2 are percentage who are successfully immunized after 1 or 2 doses (between 0 and 100)
    #	proportion.vaccinatable is proportion of the population who ar eligible to receive the vaccine
  #	proportion.of.these.to.vax is proportion of the population eligible to actually be vaccinated
  #  note that if pregnant.vaccinated is true then proportion.pregnant must be <= proportion.vaccinatable 
  # if pregnant.vaccinated is false then proportion.pregnant must be <= (1-proportion.vaccinatable )
  
  #   note that these numbers don't account for binomial variation in number of deaths etc ...so we are looking at expected values 
  # assumes posterior.samples,  VE1.post.sample and VE2subset.post.sample  is globally defined (latter 2 are vaccine efficacies after 1 or 2 doses)
  n.samples<-length(posterior.samples[,1])
  runs.to.plot<-  c(1:n.samples)[order(runif(n.samples))][1:num.its]
  results<-data.frame(infections=0, cases=0, cases.pregnant =0, cases.not.pregnant=0 ,deaths.preg=0, deaths.nonpreg=0 )
  #if((stop.V1-start.V1)<0.019999 ) stop("Start and stop vaccination times too close together and very high vaccination rates may make estimates unstable. Consider revising code if you need to consider this scenario.")
  #if((stop.V2-start.V2)<0.019999 ) stop("Start and stop vaccination times too close together and very high vaccination rates may make estimates unstable. Consider revising code if you need to consider this scenario.")
  for(j in 1:length(runs.to.plot)){
    i<-runs.to.plot[j]
    if(is.na(vacc.eff1))   VE1<-VE2.post.sample[i] else VE1<-vacc.eff1
    if(is.na(vacc.eff2))   VE2<-VE3subset.post.sample[i] else VE2<-vacc.eff2
  #    print(c("aaaa"))
   # print(VE1)
   # print(VE2)
 #   print(c("xxxx"))
    
    prob.symptomatic.infection.pregnant <-get.prob.symptomatic.if.infected(assumption = 1, pregnant = TRUE)
    prob.symptomatic.infection.not.pregnant <-prob.symptomatic.infection.pregnant  # baseline assumption is no difference (lacking data to definitively say)
    
    prob.death.not.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=FALSE) # given symptomatic infection
    prob.death.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=TRUE) # given symptomatic infection
    
    params<-posterior.samples[i,]  # 
    if(campnumber==1) campR0<-params$R01
    if(campnumber==2) campR0 <-params$R02
    if(campnumber==3) campR0 <-params$R03
 
    parameters<-c(
      beta.min= campR0* params$rho , #  currently no seasonality
      beta.max= campR0* params$rho,
      phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
      gamma=params$gamma,
      rho=params$rho,
      p.imp.S=1,  # proportion imported cases who are susceptible etc
      p.imp.E=0,
      p.imp.I=0,
      p.imp.R=0,
      vacc.eff1=VE1/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
      vacc.eff2=VE2/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose

      proportion.vaccinatable= proportion.vaccinatable,
      proportion.to.vaccinate= proportion.of.these.to.vax
      )
    xstart<-c(
      SNV=S0*(1-proportion.vaccinatable),
      SV=S0*proportion.vaccinatable,
      ENV=1,
      EV=0,
      INV=I0 * (1- proportion.vaccinatable ),
      IV=I0 * proportion.vaccinatable ,
      RNV=R0 * (1- proportion.vaccinatable ),
      RV=R0 * proportion.vaccinatable ,
      VS=VS0,
      V1=V10,
      V2=V20,
      CumIncNV=0,
      CumIncV=0
      )
    
    
    if(!is.na(cases.b4.vacc)) {
      #   then  run model and determine when first day when there  are  cases.b4.vacc - and delay immunisation times by this much
      params.novax<-parameters
      params.novax[10] <-0 #params.novax$vacc.eff1<-0
      params.novax[11] <-0 #p#params.novax$vacc.eff2<-0
      #  print(params.novax)
      no.vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3, params.novax)
      no.vax.cum.inc<-(no.vax.out[,13]+no.vax.out[,14])*params$proportion.symptomatic
      #print(no.vax.cum.inc)
      n.time.points<-length(no.vax.cum.inc)
      first.time.gt.Xcases<-no.vax.out[min(n.time.points,1+findInterval(c(-Inf, cases.b4.vacc), no.vax.cum.inc)[2]),1]
      print(c("xxx",first.time.gt.Xcases))
      new.immunisation.times <-times[findInterval(immunisation.times+ first.time.gt.Xcases +immun.delay1,times)]
      
      out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,parameters,events=list(func=immunisation.event.forv3,time= new.immunisation.times))

    } else {
      #print(parameters)
      new.immunisation.times <-times[findInterval(immunisation.times+ immun.delay1,times)]
      out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,parameters,events=list(func=immunisation.event.forv3,time= new.immunisation.times))
    }
    # above performs instantaneous immunisation at immunisation.tiems
    total.infected<-max(out[,13]) + max(out[,14])
    total.infected.vaccinatable<-max(out[,14])
    total.infected.notvaccinatable<-max(out[,13])
    
    total.symptomatic<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    total.cases<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    
    # Need to split this output into vaccinatable and not vaccinatable populations so we get total infections and symptomatic in each
    #  using pregnant.vaccinated  proportion.pregnant
    # prob.symptoms.if.not.pregnant<- params$proportion.symptomatic/(1+proportion.pregnant*(relative.risk.of.symptoms.given.infected.if.pregnant-1))
    
    if(pregnant.vaccinated) {
      Number.infected.pregnant<-total.infected.vaccinatable*proportion.pregnant/ proportion.vaccinatable
      prob.pregnant.infected<-total.infected.vaccinatable/(S0*proportion.vaccinatable)
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable-proportion.pregnant
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)
      
    } else {
      Number.infected.pregnant<-total.infected.notvaccinatable*proportion.pregnant/ (1-proportion.vaccinatable)
      prob.pregnant.infected<-total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)-proportion.pregnant
    }
    Number.infected.not.pregnant<-total.infected-Number.infected.pregnant
    
    # note that for the output we need to distinguish cases (those which are reported) and symptomatic infections (which may or may not be reported)
    
    relative.risk.of.beingcase.given.infected.if.pregnant <- 2  # assume that those who are pregnant are twice as likely to be reported as cases (irresepective of relative risk of symptoms)  
    
    #note in line below  params$proportion.symptomatic is a bit of misnomer - it really represents the proportion of infections which are reported as cases
    
    prob.case.if.not.pregnant<-params$proportion.symptomatic*(1+ Number.infected.not.pregnant/Number.infected.pregnant) /(relative.risk.of.beingcase.given.infected.if.pregnant + Number.infected.not.pregnant/Number.infected.pregnant)
    prob.case.if.pregnant<- relative.risk.of.beingcase.given.infected.if.pregnant *prob.case.if.not.pregnant
    
    Number.symptomatic.pregnant<-prob.symptomatic.infection.pregnant* Number.infected.pregnant
    Number.symptomatic.not.pregnant<-prob.symptomatic.infection.not.pregnant* Number.infected.not.pregnant
    
    total.symptomatic<-Number.symptomatic.pregnant + Number.symptomatic.not.pregnant
    
    Number.cases.pregnant = prob.case.if.pregnant * Number.infected.pregnant
    Number.cases.not.pregnant<-prob.case.if.not.pregnant * Number.infected.not.pregnant
    Number.deaths.pregnant<-prob.death.pregnant*Number.symptomatic.pregnant 
    Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.symptomatic.not.pregnant 
    
    # print(j)
    results[j,]<-c(infections=total.infected, cases=total.cases, cases.pregnant= Number.cases.pregnant , cases.not.pregnant= Number.cases.not.pregnant, symptomatic=total.symptomatic,  deaths.preg=Number.deaths.pregnant, deaths.nonpreg=Number.deaths.not.pregnant)
    
  } # end for 
  results
}



run.vaccine.scenario.v3<-function(num.its=1,  times=seq(0,3,1/52), cases.b4.vacc=NA,  immunisation.times =Inf,vacc.eff1=NA,vacc.eff2=NA,immun.delay1=2/52, proportion.vaccinatable=0.8,proportion.of.these.to.vax=1,pregnant.vaccinated=FALSE, proportion.pregnant=0.04, campnumber=1){
	# like v2 except it calculates the change in outcomes compared to no vaccine
  # Before running this make sure that the dataframe posterior.samples holds samples from the distribution of parameters values that we want to use inthe simulations.
  # In this version, immunisation takes place instantaneously at immunisation.times
  # but there is a delay given by immun.delay1 (half a month by default) before immunity .
  # cases.b4.vacc is equal to the number of cases before vaccine starts. 
  # If cases.b4.vacc is NA there is no reactive vaccination
  #   num.its is the number of iterations
  #   vacc.eff1 and vacc.eff2 are percentage who are successfully immunized after 1 or 2 doses (between 0 and 100)
  #	proportion.vaccinatable is proportion of the population who ar eligible to receive the vaccine
  #	proportion.of.these.to.vax is proportion of the population eligible to actually be vaccinated
  #  note that if pregnant.vaccinated is true then proportion.pregnant must be <= proportion.vaccinatable 
  # if pregnant.vaccinated is false then proportion.pregnant must be <= (1-proportion.vaccinatable )
  
  #   note that these numbers don't account for binomial variation in number of deaths etc ...so we are looking at expected values 
  # assumes posterior.samples,  VE1.post.sample and VE2subset.post.sample  is globally defined (latter 2 are vaccine efficacies after 1 or 2 doses)
  n.samples<-length(posterior.samples[,1])
  runs.to.plot<-  c(1:n.samples)[order(runif(n.samples))][1:num.its]
  results<-data.frame(infections=0, cases=0, cases.pregnant =0, cases.not.pregnant=0 ,deaths=0, deaths.preg=0, deaths.nonpreg=0 )
 # store percentage reduction in above associated with vaccine
 
  for(j in 1:length(runs.to.plot)){
    i<-runs.to.plot[j]
     # by default assume first vaccine effect is after first two doses so no effect from 1 dose. Timing of VE1 (which is really second vaccine dose) should therefore account for time needed to give first dose, delay to second dose, and delay to immunity 
    if(is.na(vacc.eff1))   VE1<-VE2.post.sample[i] else VE1<-vacc.eff1 
    if(is.na(vacc.eff2))   VE2<-VE3.post.sample[i] else VE2<-vacc.eff2

    
    prob.symptomatic.infection.pregnant <-get.prob.symptomatic.if.infected(assumption = 1, pregnant = TRUE)
    prob.symptomatic.infection.not.pregnant <-prob.symptomatic.infection.pregnant  # baseline assumption is no difference (lacking data to definitively say)
    
    prob.death.not.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=FALSE) # given symptomatic infection
    prob.death.pregnant<-get.prob.death.if.symptomatic(assumption=1, pregnant=TRUE) # given symptomatic infection
    
    
    params<-posterior.samples[i,]   

    if(campnumber==1) campR0 <-params$R01
    if(campnumber==2) campR0 <-params$R02
    if(campnumber==3) campR0 <-params$R03
 
    parameters<-c(
      beta.min= campR0* params$rho , #  currently no seasonality
      beta.max= campR0* params$rho,
      phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
      gamma=params$gamma,
      rho=params$rho,
      p.imp.S=1,  # proportion imported cases who are susceptible etc
      p.imp.E=0,
      p.imp.I=0,
      p.imp.R=0,
      vacc.eff1=VE1/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose
      vacc.eff2=VE2/100,   # vaccine efficacy: in this implementation we define this as the proportion of people successfully immunized by one dose

      proportion.vaccinatable= proportion.vaccinatable,
      proportion.to.vaccinate= proportion.of.these.to.vax
      )
  #  print("parameter vector for this run")
  #  print(parameters)
    

    xstart<-c(
      SNV=S0*(1-proportion.vaccinatable),  # susceptibles not eligible to be vaccinated given current policy
      SV=S0*proportion.vaccinatable,       # susceptibles  eligible to be vaccinated
      ENV=1,
      EV=0,
      INV=I0 * (1- proportion.vaccinatable ),
      IV=I0 * proportion.vaccinatable ,
      RNV=R0 * (1- proportion.vaccinatable ),
      RV=R0 * proportion.vaccinatable ,
      VS=VS0,
      V1=V10,
      V2=V20,
      CumIncNV=0,
      CumIncV=0
      )
    
      params.novax<-parameters
      params.novax[10] <-0 #params.novax$vacc.eff1<-0
      params.novax[11] <-0 #p#params.novax$vacc.eff2<-0
      #  print(params.novax)
      no.vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3, params.novax)

     # Now run the model with vaccination
    if(!is.na(cases.b4.vacc)) {
      #   then  run model and determine when first day when there  are  cases.b4.vacc - and delay immunisation times by this much
      no.vax.cum.inc<-(no.vax.out[,13]+no.vax.out[,14])*params$proportion.symptomatic
      n.time.points<-length(no.vax.cum.inc)
      first.time.gt.Xcases<-no.vax.out[min(n.time.points,1+findInterval(c(-Inf, cases.b4.vacc), no.vax.cum.inc)[2]),1]

      new.immunisation.times <-times[findInterval(immunisation.times+ first.time.gt.Xcases +immun.delay1,times)]
      out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,parameters,events=list(func=immunisation.event.forv3,time= new.immunisation.times))

    } else {
      #print(parameters)
      immunisation.times <-times[findInterval(immunisation.times + immun.delay1,times)] #this ensures immunisation time is at included time point
      out<-ode(y=xstart,times=times, func=HepEwithVAxmod.v3,parameters,events=list(func=immunisation.event.forv3,time= immunisation.times))
    }
   # print("Summary of model without vaccine")
   # print(summary(no.vax.out))
   # print("Summary of model with vaccine")
   # print(summary(out))

    # above performs instantaneous immunisation at immunisation.tiems
    total.infected<-max(out[,13]) + max(out[,14])
    total.infected.vaccinatable<-max(out[,14])
    total.infected.notvaccinatable<-max(out[,13])
    
    total.symptomatic<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    total.cases<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    
    # Need to split this output into vaccinatable and not vaccinatable populations so we get total infections and symptomatic in each
    #  using pregnant.vaccinated  proportion.pregnant
    # prob.symptoms.if.not.pregnant<- params$proportion.symptomatic/(1+proportion.pregnant*(relative.risk.of.symptoms.given.infected.if.pregnant-1))
    
    if(pregnant.vaccinated) {
      Number.infected.pregnant<-total.infected.vaccinatable*proportion.pregnant/ proportion.vaccinatable
      prob.pregnant.infected<-total.infected.vaccinatable/(S0*proportion.vaccinatable)
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable-proportion.pregnant
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)
      
    } else {
      Number.infected.pregnant<-total.infected.notvaccinatable*proportion.pregnant/ (1-proportion.vaccinatable)
      prob.pregnant.infected<-total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)-proportion.pregnant
    }
    Number.infected.not.pregnant<-total.infected-Number.infected.pregnant
    
    prob.symptoms.if.not.pregnant<-total.infected*params$proportion.symptomatic/(Number.infected.not.pregnant+ Number.infected.pregnant* relative.risk.of.symptoms.given.infected.if.pregnant)
    #prob.symptoms.if.not.pregnant<-params$proportion.symptomatic*(1+ Number.infected.not.pregnant/Number.infected.pregnant)/(relative.risk.of.symptoms.given.infected.if.pregnant+Number.infected.not.pregnant/Number.infected.pregnant)
    # prob.not.pregnant.infected<-proportion.vaccinatable.and.not.pregnant*total.infected.vaccinatable/(S0*proportion.vaccinatable) +  #proportion.not.vaccinatable.and.not.pregnant* total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
    
    #prob.symptoms.if.not.pregnant<-params$proportion.symptomatic/((1-proportion.pregnant)* prob.not.pregnant.infected +proportion.pregnant*  #relative.risk.of.symptoms.given.infected.if.pregnant* prob.pregnant.infected) 
    prob.symptoms.if.pregnant<- relative.risk.of.symptoms.given.infected.if.pregnant*prob.symptoms.if.not.pregnant
    
    Number.cases.pregnant = prob.symptoms.if.pregnant * Number.infected.pregnant
    Number.cases.not.pregnant<-prob.symptoms.if.not.pregnant * Number.infected.not.pregnant
    Number.deaths.pregnant<-prob.death.pregnant*Number.cases.pregnant 
    Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.cases.not.pregnant  
    Number.deaths<-Number.deaths.pregnant+ Number.deaths.not.pregnant
# store results for vax policy here
    results.vax<-c(infections=total.infected, cases=total.symptomatic, cases.pregnant= Number.cases.pregnant , cases.not.pregnant= Number.cases.not.pregnant, deaths=Number.deaths ,deaths.preg=Number.deaths.pregnant, deaths.nonpreg=Number.deaths.not.pregnant)


# now calculate corresponding numbers using same parameters except when there is no vaccine

    total.infected<-max(no.vax.out[,13]) + max(no.vax.out[,14])
    total.infected.vaccinatable<-max(no.vax.out[,14])
    total.infected.notvaccinatable<-max(no.vax.out[,13])
    
    total.symptomatic<-total.infected*params$proportion.symptomatic #here proportion.symptomatic is in whole pop (including pregant)
    
    # Need to split this output into vaccinatable and not vaccinatable populations so we get total infections and symptomatic in each
    #  using pregnant.vaccinated  proportion.pregnant
    # prob.symptoms.if.not.pregnant<- params$proportion.symptomatic/(1+proportion.pregnant*(relative.risk.of.symptoms.given.infected.if.pregnant-1))
    
    if(pregnant.vaccinated) {
      Number.infected.pregnant<-total.infected.vaccinatable*proportion.pregnant/ proportion.vaccinatable
      prob.pregnant.infected<-total.infected.vaccinatable/(S0*proportion.vaccinatable)
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable-proportion.pregnant
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)
      
    } else {
      Number.infected.pregnant<-total.infected.notvaccinatable*proportion.pregnant/ (1-proportion.vaccinatable)
      prob.pregnant.infected<-total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
      proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable
      proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)-proportion.pregnant
    }
    Number.infected.not.pregnant<-total.infected-Number.infected.pregnant
    
    prob.symptoms.if.not.pregnant<-total.infected*params$proportion.symptomatic/(Number.infected.not.pregnant+ Number.infected.pregnant* relative.risk.of.symptoms.given.infected.if.pregnant)
    #prob.symptoms.if.not.pregnant<-params$proportion.symptomatic*(1+ Number.infected.not.pregnant/Number.infected.pregnant)/(relative.risk.of.symptoms.given.infected.if.pregnant+Number.infected.not.pregnant/Number.infected.pregnant)
    # prob.not.pregnant.infected<-proportion.vaccinatable.and.not.pregnant*total.infected.vaccinatable/(S0*proportion.vaccinatable) +  #proportion.not.vaccinatable.and.not.pregnant* total.infected.notvaccinatable/(S0*(1-proportion.vaccinatable))
    
    #prob.symptoms.if.not.pregnant<-params$proportion.symptomatic/((1-proportion.pregnant)* prob.not.pregnant.infected +proportion.pregnant*  #relative.risk.of.symptoms.given.infected.if.pregnant* prob.pregnant.infected) 
    prob.symptoms.if.pregnant<- relative.risk.of.symptoms.given.infected.if.pregnant*prob.symptoms.if.not.pregnant
    
    Number.cases.pregnant = prob.symptoms.if.pregnant * Number.infected.pregnant
    Number.cases.not.pregnant<-prob.symptoms.if.not.pregnant * Number.infected.not.pregnant
    Number.deaths.pregnant<-prob.death.pregnant*Number.cases.pregnant 
    Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.cases.not.pregnant 
    Number.deaths<-Number.deaths.pregnant+ Number.deaths.not.pregnant

   # print(Number.infected.pregnant) 
  #    print("aaa")
   # print(prob.symptoms.if.pregnant)
  #      print("bbb")
  #  print( Number.infected.not.pregnant)
  #      print("ccc")
  #  print( prob.symptoms.if.not.pregnant)
 #   print("xxx")
    # print(j)
    # output is % reduction in cases deaths etc due to vax
    results.novax<-c(infections=total.infected, cases=total.symptomatic, cases.pregnant= Number.cases.pregnant , cases.not.pregnant= Number.cases.not.pregnant, deaths= Number.deaths, deaths.preg=Number.deaths.pregnant, deaths.nonpreg=Number.deaths.not.pregnant)

    results[j,]<-100*(results.novax-results.vax)/results.novax
    
  } # end for 
  results
}




