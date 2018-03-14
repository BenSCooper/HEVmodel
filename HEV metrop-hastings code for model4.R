# HEV metrop-hastings code for model4.R 
#( SEIIR model with moderately informative priors where transmission rate from first and second I compartment is the same)
# Ben Cooper 2016

# Fit the model to all three camps together using Metropolis-Hastings Algorithm
# estimating common natural history parameters and the proportion symptomaticc
# Model 4 does this with an SEIIR model with equal transmissibility in the two I compartments
#(by default with same parameters for leaving I1 and I2 so infectious period has erlang distribution)
# Model 5 does this with  an SEIIR model with two levels of transmissibility (early and latet infectious)
#camp specific R0s, camp-specific time of first case but estimating 
# gamma and proportion sympomatic which are the same for all three camps.
# Assume data follow Negative binomial distribution around predicted means 
# and use compiled C code to run the ODE model (much faster)
# ************************************************************************************************


library(deSolve)

# first compiled and load the C code for the SEIR model

system("R CMD SHLIB SEEIR.c")
dyn.load("SEEIR.so")

system("R CMD SHLIB SEIIR.c")
dyn.load("SEIIR.so")

system("R CMD SHLIB SEIRforced.c")
dyn.load("SEIRforced.so")

Uganda<-read.csv("UgandaHepECases.csv")
Uganda[is.na(Uganda)]<-0

logit<-function(y){log(y/(1-y))}
inv.logit<-function(x){(exp(x)/(1+exp(x)))}

HepEmod.nowatsan <- function(t, x, parms) { 
  #modified version of model to allow transmission to vary according to Watsan interventions
  # function returning a list holding the model equations and the rates of import and export of people into the opulation at diferent time points. 
  # Assumes a fixed population size and no vaccination
  # t is time  
  #x is the vector of state variables (S, E, I, R) 
  # parms is the vector of parameters
  
  with(as.list(c(parms, x)), {
    # beta.max is max  transmission parameter (when no water)
    # beta.min is min seasonal transmission paramter (when unlimited water)
    # r.water is the rate term in the model relating watsan intervention effects on water to beta
    # r.latrines is the rate term in the model relating watsan intervention effects on latrines to beta
    
    # phi is phase angle which determines what time of year peak(s) occurs
    # N is current population size
    # gamma is rate of progression to exposed to infectious
    # rho is rate of progression from infectious to recovered.
    
    N<-S+E+I+R 
    
    # beta<-seasonal.transmission.parameter(beta.max,beta.min,phi,t)
    # w<-getwaterperp(t)
    # #l<-getlatrinesperp(t)
    # beta<-get.transmisison.parameter.given.watsan(beta.max,beta.min,r.water,w,r.latrines,l)
    # Derivatives 
    beta<-beta.max
    dS <- -beta*S*I/N 
    dE <-  beta*S*I/N - gamma*E  
    dI <- gamma*E - rho*I  
    dR <- rho*I 
    dCumInc<-gamma*E
    list(c(dS,dE,dI,dR,dCumInc))
    
  })  # end with
} # end HepEmod


getLLHepEmod.mh.combined.no.watsan<-function(params, x1,x2,x3, plot=FALSE,verbose=FALSE, SEEIR=FALSE, SEIIR=TRUE){
  #The version is designed for use with Metropolis Hastings algorithm
  # assumes no effect of watsan
  # x1, x2, x3 are vectors of state variables (S, E, I, R) for each camp at start time (assume this vector is known )
  # if SEEIR ==TRUE thihs uses the dll for the SEEIR model which assumes two equal comparments (by default with same rate of leaving
  # so latent period has an erlang distribution) 
  # if SEIIR ==TRUE thihs uses the dll for the SEEIR model which assumes two equal comparments (by default with same rate of leaving
  # so latent period has an erlang distribution) 
 
  
  with(as.list(c(params, x1,x2,x3)), {
    # R01 is R0  in camp 1 
    # R02 is R0  in camp 2 
    # R03 is R0  in camp 3 
    # time.first.infected1  is time of first case in year 1 in camp 1 (scale is years)
    # time.first.infected2  is time of first case in year 1 in camp 2 (scale is years)
    # time.first.infected3  is time of first case in year 1 in camp 3 (scale is years)
    # gamma is rate of progression to exposed to infectious
    # rho is rate of progression from infectious to recovered (log(1/infectious period).
    # proportion.symptomatic is  proporotion of cases seen which are symptomatic
    # weeks.from.infection.symptoms - delay till symptoms show (assumed constant)
    # disp is the dispersion parameter for negative binomial likelihood. If NA use a Poisson likelihood
    wk.first.infected1<-round(time.first.infected1*52)
    wk.first.infected2<-round(time.first.infected2*52)
    wk.first.infected3<-round(time.first.infected3*52)
    wk<-1/52
    # output first time point then regular weekly timesteps 
    times1<- seq(ceiling(time.first.infected1*52)/52, 157/52, by = wk) 
    times2<- seq(ceiling(time.first.infected2*52)/52, 157/52, by = wk) 
    times3<- seq(ceiling(time.first.infected3*52)/52, 157/52, by = wk)
    if(time.first.infected1!=times1[1]) times1<-c(time.first.infected1,times1)
    if(time.first.infected2!=times2[1]) times2<-c(time.first.infected2,times2)
    if(time.first.infected3!=times3[1]) times3<-c(time.first.infected3,times3)

    parameters1<-c(
      beta=R01 *rho,
      gamma=gamma,
      rho=rho)
    
    parameters2<-c(
      beta=R02 *rho,
      gamma=gamma,
      rho=rho)
    
    parameters3<-c(
      beta=R03 *rho,
      gamma=gamma,
      rho=rho)
    
    if(SEEIR==FALSE & SEIIR==FALSE){
     out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
     out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
     out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
    } else if(SEEIR==TRUE & SEIIR==FALSE) {
      parameters1<-c(
        beta=R01 *rho,
        gamma1=gamma,
        gamma2=gamma,
        rho=rho)
      
      parameters2<-c(
        beta=R02 *rho,
        gamma1=gamma,
        gamma2=gamma,
        rho=rho)
      
      parameters3<-c(
        beta=R03 *rho,
        gamma1=gamma,
        gamma2=gamma,
        rho=rho)
      out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
   
     
    } else if(SEEIR==FALSE & SEIIR==TRUE) { # in this case set parameter for model where beta1 =beta2 and rho1 = rho2
        parameters1<-c(
        beta1=R01 *rho/2,
        beta2=R01 *rho/2,
        gamma=gamma,
        rho1=rho,
        rho2=rho)
      
      parameters2<-c(
        beta1=R02 *rho/2,
        beta2=R02 *rho/2, 
        gamma=gamma,
        rho1=rho,
        rho2=rho)
      
      parameters3<-c(
        beta1=R03 *rho/2,
        beta2=R03 *rho/2, 
        gamma=gamma,
        rho1=rho,
        rho2=rho)
        
      out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
    	
    } else {
    	stop("Something has gone wrong. This should not happen.")
    }
    
    cuminc1<-out1[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
    cuminc2<-out2[,6][-1]
    cuminc3<-out3[,6][-1]
    
    n1<-length(cuminc1)
    n2<-length(cuminc2)
    n3<-length(cuminc3)
    
    weeklyinc1<-c(cuminc1[1],cuminc1[2:n1]-cuminc1[1:(n1-1)]) # weekly incidence
    weeklyinc2<-c(cuminc2[1],cuminc2[2:n2]-cuminc2[1:(n2-1)]) # weekly incidence
    weeklyinc3<-c(cuminc3[1],cuminc3[2:n3]-cuminc3[1:(n3-1)]) # weekly incidence
    

    obs1<-OBSERVED.CASES1[(LATEST.START.WEEK1+1):length(OBSERVED.CASES1)] # only use obs after LATEST.START.WEEK in LL
    obs2<-OBSERVED.CASES2[(LATEST.START.WEEK2+1):length(OBSERVED.CASES2)] # only use obs after LATEST.START.WEEK in LL
    obs3<-OBSERVED.CASES3[(LATEST.START.WEEK3+1):length(OBSERVED.CASES3)] # only use obs after LATEST.START.WEEK in LL
    wk.inc1<-weeklyinc1[1:(n1-weeks.from.infection.symptoms)] # this gives total incident infectious caes - a proportion of whcih will be symptomatic & reported after a delay of weeks.from.infection.symptoms
    wk.inc2<-weeklyinc2[1:(n2-weeks.from.infection.symptoms)]
    wk.inc3<-weeklyinc3[1:(n3-weeks.from.infection.symptoms)]
    
    # now discard values of incidence before the  first obs we consider in the likelihood
    numbertochop1<-length(wk.inc1)-length(obs1)
    numbertochop2<-length(wk.inc2)-length(obs2)
    numbertochop3<-length(wk.inc3)-length(obs3)
    
    wk.inc1<-wk.inc1[(numbertochop1+1):length(wk.inc1)]
    wk.inc2<-wk.inc2[(numbertochop2+1):length(wk.inc2)]
    wk.inc3<-wk.inc3[(numbertochop3+1):length(wk.inc3)]
    
    # now obs and wk.inc vectors should be the same length and both start at LATEST.START.WEEK
    if(verbose){
      print(c(numbertochop1,numbertochop2,numbertochop3))
      print(obs1/wk.inc1)  # should be less than one for a feasible value
      print(obs2/wk.inc2)
      print(obs3/wk.inc3)
    }  
    disp<-params$disp
    if(is.na(disp)){ # then use Poisson likelihood
      LLterms1<-dpois(obs1, wk.inc1 * proportion.symptomatic,log=TRUE)
      LLterms2<-dpois(obs2, wk.inc2 * proportion.symptomatic,log=TRUE)
      LLterms3<-dpois(obs3, wk.inc3 * proportion.symptomatic,log=TRUE)
    } else {  # then use negbin likelihood
      LLterms1<-dnbinom(obs1,mu= wk.inc1 * proportion.symptomatic,size=disp, log=TRUE)
      LLterms2<-dnbinom(obs2,mu= wk.inc2 * proportion.symptomatic,size=disp,log=TRUE)
      LLterms3<-dnbinom(obs3,mu= wk.inc3 * proportion.symptomatic,size=disp,log=TRUE)
    }  

    LL<-sum(LLterms1)+sum(LLterms2)+sum(LLterms3) 
    
    if(plot==TRUE){
      expected.cases1<-proportion.symptomatic*wk.inc1
      expected.cases2<-proportion.symptomatic*wk.inc2
      expected.cases3<-proportion.symptomatic*wk.inc3
      ymax=max(obs1,obs2,obs3,expected.cases1,expected.cases2,expected.cases3)
      plot(obs1,type='p',col="red",ylim=c(0,ymax))
      lines(expected.cases1,type='l',col="red")
      points(obs2,col="blue",ylim=c(0,ymax))
      lines(expected.cases2,type='l',col="blue")
      points(obs3,col="green",ylim=c(0,ymax))
      lines(expected.cases3,type='l',col="green")
    }    
    return(list(LL=LL))
  })   
}




LATEST.START.WEEK1<-50 # last week in year one when epidemic could start. (use 50 for camp1, 43 for camp 2 and 47 for camp 3) to allow for delay to first cas
LATEST.START.WEEK2<-43
LATEST.START.WEEK3<-47
#LATEST.START.WEEK3<-49
# should be weeks.from.infection.symptoms before first case is seen

OBSERVED.CASES1<-Uganda[,2]
OBSERVED.CASES2<-Uganda[,3]
OBSERVED.CASES3<-Uganda[,4]
watsantimes<-c(71/52,87/52,81/52)  # times of watsan intervention start in time units of year (see UgandaWaterPerPerson.csv)
Uganda.wpp<-read.csv("UgandaWaterPerPerson.csv")
Uganda.san<-read.csv("Uganda_Latrines.csv")



N1<-16689  # Population size of Agoro
N2<-10442  # Population size of MadiOpei
N3<-10555  # Population size of Paloga
wks.to.symptoms<-1  # i.e. delay in weeks from onset of infectiousness to onset of symptoms

# initial values
init.latent.period<-2/52  # in units of 1 year (values between 2 and 7 weeks are plausible apparently)
init.infectious.period<-2/52  # values of about 2 weeks are plausible
init.propSymptoms<-1/7   # Assumed value for proportion symptomatic
init.R0<-4

E0<-1
I0<-0
R0<-0
CumInc0<-0
xstart1<-c(S=N1-1,E=E0,I=I0, R=R0, CumInc=CumInc0)
xstart2<-c(S=N2-1,E=E0,I=I0, R=R0, CumInc=CumInc0)
xstart3<-c(S=N3-1,E=E0,I=I0, R=R0, CumInc=CumInc0)

###### Metropolis-Hastings sampling
# we want to make inference about the following parameters: 
# R01, R02, R03 - the basic reproduction numbers in each of the three camps
# gamma - rate of progression from latently infected to infectious (1/latent.period)
# rho - rate of recovery for those who are  infectious (1/infectious.period)
# propSymptoms - proportion of people who have symptoms
# time.firstcase1 time.firstcase2 time.firstcase3  - time of first case in each of the three camps

# 1. Define Prior distributions
# we give R0 priors on [1,10] 
# we give gamma a diffuse gamma prior 
# we give rho a diffuse gamma prior 
# time.firstcase1 priors on [1,LATEST.START.WEEK1/52]  (time in units of years) 
# we give  propSymptoms a flat beta prior beta(1,1) which takes values between 0 and 1

# Uninformative priors below
R0.prior.density<-function(x){ dunif(x, min=1, max=20)}  #note we can either parameterise for R0 or for beta, but not both...
#beta.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)} #as R0=beta/rho
beta.prior.density<-function(x){ dunif(x, 0.0001,500)} 
gamma.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
rho.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
time.first.infected1.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK1/52)}
time.first.infected2.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK2/52)}
time.first.infected3.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK3/52)}
propSymptoms.prior.beta.shape1<-1  # parmeters for propSymptoms.prior
propSymptoms.prior.beta.shape2<-1  # parmeters for propSymptoms.prior
propSymptoms.prior.density<-function(x){ dbeta(x, propSymptoms.prior.beta.shape1,propSymptoms.prior.beta.shape2)}
disp.prior.density<-function(x){ dunif(x, min=0.1, max=10000)}

# Informative priors below
 
 R0.prior.density<-function(x){ dunif(x, min=1, max=30)}  #note we can either parameterise for R0 or for beta, but not both...
 #beta.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)} #as R0=beta/rho
# beta.prior.density<-function(x){ dunif(x, 0.0001,500)} #as


 gamma.prior.density<-function(x){ dgamma(x, shape=1,rate=0.0932)} #rate is 34/365. Conjugate based on one observation with a mean latent period of 34 days

# rho.prior.density<-function(x){ dgamma(x, shape=11,rate=1.145)}
 rho.prior.density<-function(x){ dgamma(x, shape=22,rate=1.145)} # twice rate so half duration in each I comparament - but two Is

# #above  based on Takahashi et al 2007: 11 patients with a mean duration of infectiousness of 38 days (assuming 1 week infectious before symptoms)
# #conjugate gamma shape and rate are shape=a+n,rate=r +n*x if prior has shape a and rate r
# #and if n observation with a mean of x. Time units are in years here


# 2. Define functions to propose new values for parameters
# Adjust variance of the proposals to tune the algorithm (ideally we want to accept about a fifth of proposals for each variable)
# Note that these functions may propose values outside possible ranges, but since prior distributions
# will have zero probability at these values they have no chance of being accepted
generate.proposal.R0<-function(currentR0){return(currentR0 + rnorm(1,0,.1))}
generate.proposal.beta<-function(currentbeta){return(currentbeta + rnorm(1,0,5))}
generate.proposal.gamma<-function(currentgamma){return(currentgamma + rnorm(1,0,.4))}
generate.proposal.rho<-function(currentrho){return(currentrho + rnorm(1,0,.4))}
generate.proposal.time.first.infected<-function(currenttime){return(currenttime + rnorm(1,0,.01))}
generate.proposal.propSymptoms<-function(currentprop){return(currentprop + rnorm(1,0,.01))} #   used if using Poisson likelihood (if binom likilihood use a Gibbs update for this one)
generate.proposal.disp<-function(currentdisp){return(currentdisp + rnorm(1,0,2))}
generate.proposal.watsan.effect<-function(currentw){return(currentw + rnorm(1,0,50))}

# initial values
parameters<-list(
  R01=4, 
  R02=5, 
  R03=6, 
  gamma=10.4,
  rho=20,
  time.first.infected1=0.87,
  time.first.infected2=0.7,
  time.first.infected3=0.83,
  proportion.symptomatic=init.propSymptoms,
  weeks.from.infection.symptoms=1,
  disp=1  # dispersion for negative binomial likelihood. Set to NA to use a Poisson likelihood
  #  w1=0, #watsan effect in camp1  - comment out if not using watsan model
  #  w2=0, #watsan effect in camp2  - comment out if not using watsan model
  #w3=0# watsan effect in camp3     - comment out if not using watsan model
  )  
#  parameters$disp<-NA # setting to NA forces a Poisson likelihood
parameters<-unlist(parameters)
parameters<-as.data.frame(t(parameters))
#thin<-10

# can do about 20000 runs in one hour - we need to do about 3 million runs in total so really a week long runs are needed.

thin<-100
N<-3000000

last.posterior.samples<-posterior.samples

m<-dim(last.posterior.samples)[1]
parameters<-last.posterior.samples[m,]

proposed.parameters<-parameters

# Below the SEIIR option, if set to true, causes likelihood to be based on SEIIR model, where time is in two 
# I compartments is equal (giving an erlang distributed infectious period - though if using this should use different priors 
# for the rho parameters )
current.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

#store results from the posterior sample with one column for each parameter & one for LL
posterior.samples<-as.data.frame(matrix(c(as.numeric(parameters),current.loglikelihood$LL),nrow=1))
colnames(posterior.samples)<-c(names(as.vector(parameters)),"LL")
beta1moves.accepted<-0    # note that we either update R0s or betas but not both
beta2moves.accepted<-0    # note that we either update R0s or betas but not both
beta3moves.accepted<-0    # note that we either update R0s or betas but not both
R01moves.accepted<-0    # note that we either update R0s or betas but not both
R02moves.accepted<-0    # note that we either update R0s or betas but not both
R03moves.accepted<-0    # note that we either update R0s or betas but not both
beta.moves.accepted<-0
gamma.moves.accepted<-0
rho.moves.accepted<-0
time.infected1.moves.accepted<-0
time.infected2.moves.accepted<-0
time.infected3.moves.accepted<-0
propSymptoms.moves.accepted<-0
disp.moves.accepted<-0
w1.moves.accepted<-0
w2.moves.accepted<-0
w3.moves.accepted<-0

for(i in 1:N){
  #  updates to beta1, beta1, and beta3 (not blocked together)
  #proposed.parameters$R01<-generate.proposal.R0(parameters$R01) # proposed update
  proposed.beta1<-
    (parameters$R01*parameters$rho) # proposed update
  proposed.parameters$R01<-proposed.beta1/parameters$rho
  #current.prior.density<-R0.prior.density(parameters$R01)
  current.prior.density<-beta.prior.density(parameters$R01*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta1)
  
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      beta1moves.accepted<-beta1moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  #proposed.parameters$R02<-generate.proposal.R0(parameters$R02) # proposed update
  proposed.beta2<-generate.proposal.beta(parameters$R02*parameters$rho) # proposed update
  proposed.parameters$R02<-proposed.beta2/ parameters$rho
  #current.prior.density<-R0.prior.density(parameters$R01)
  current.prior.density<-beta.prior.density(parameters$R02*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta2)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      beta2moves.accepted<-beta2moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  
  proposed.beta3<-generate.proposal.beta(parameters$R03*parameters$rho) # proposed update
  proposed.parameters$R03<-proposed.beta3/ parameters$rho
  #current.prior.density<-R0.prior.density(parameters$R01)
  current.prior.density<-beta.prior.density(parameters$R03*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta3)
  
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      beta3moves.accepted<-beta3moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  #  updates to gamma
  proposed.parameters$gamma<-generate.proposal.gamma(parameters$gamma) # proposed update
  current.prior.density<-gamma.prior.density(parameters$gamma)
  proposal.prior.density<-gamma.prior.density(proposed.parameters$gamma)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      gamma.moves.accepted<-gamma.moves.accepted+1
    } else proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  #  updates to rho
  proposed.parameters$rho<-generate.proposal.rho(parameters$rho) # proposed update
  proposed.parameters$R01<-parameters$R01*parameters$rho/ proposed.parameters$rho
  proposed.parameters$R02<-parameters$R02*parameters$rho/ proposed.parameters$rho
  proposed.parameters$R03<-parameters$R03*parameters$rho/ proposed.parameters$rho
  
  current.prior.density<-rho.prior.density(parameters$rho)
  proposal.prior.density<-rho.prior.density(proposed.parameters$rho)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE) 
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      rho.moves.accepted<-rho.moves.accepted+1
    }  else proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  #  updates to time.first.infected1, time.first.infected2, and time.first.infected3 (no longer blocked together)
  proposed.parameters$time.first.infected1<-generate.proposal.time.first.infected(parameters$time.first.infected1) # proposed update
  current.prior.density<-time.first.infected1.prior.density(parameters$time.first.infected1)
  proposal.prior.density<-time.first.infected1.prior.density(proposed.parameters$time.first.infected1)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      time.infected1.moves.accepted<-time.infected1.moves.accepted+1
    }  else  proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters 
  
  
  proposed.parameters$time.first.infected2<-generate.proposal.time.first.infected(parameters$time.first.infected2) # proposed update
  current.prior.density<-time.first.infected2.prior.density(parameters$time.first.infected2)
  proposal.prior.density<-time.first.infected2.prior.density(proposed.parameters$time.first.infected2)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      time.infected2.moves.accepted<-time.infected2.moves.accepted+1
    }  else  proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  proposed.parameters$time.first.infected3<-generate.proposal.time.first.infected(parameters$time.first.infected3) # proposed update
  current.prior.density<-time.first.infected3.prior.density(parameters$time.first.infected3)
  proposal.prior.density<-time.first.infected3.prior.density(proposed.parameters$time.first.infected3)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)  
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density) 
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      time.infected3.moves.accepted<-time.infected3.moves.accepted+1
    }  else  proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  #  updates to proportion.symptomatic 
  # (this can be done with a Gibbs step, as we have a beta prior and binomial likelihood (and don't change number of trials))
  #  If the prior for the probability of success (in this case reporting of a case) is beta(a, b), 
  #  and we have n trials with x successes, then the posterior for the probability of success will 
  # be beta(a+x, b+n-x). We can use this fact to use a Gibbs step in updating
  
  proposed.parameters$proportion.symptomatic<-generate.proposal.propSymptoms(parameters$proportion.symptomatic) # proposed update
  current.prior.density<-propSymptoms.prior.density(parameters$proportion.symptomatic)
  proposal.prior.density<-propSymptoms.prior.density(proposed.parameters$proportion.symptomatic)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      propSymptoms.moves.accepted<-propSymptoms.moves.accepted+1
    } else proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  # stuff below for use with binomial liklihood when we can use a conjugate prior
  #post.beta.shape1<-propSymptoms.prior.beta.shape1 + current.loglikelihood$n.successes
  #post.beta.shape2<-propSymptoms.prior.beta.shape2 + current.loglikelihood$n.trials- current.loglikelihood$n.successes
  #parameters$proportion.symptomatic=rbeta(1,post.beta.shape1,post.beta.shape2) # Gibbs update, so always accepted
  #proposed.parameters<-parameters
  # however, we still need to recalculate current likelihood as this is needed for future MH steps
  #  current.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=FALSE)
  
  #  updates to disp (dispersion parameter for negative binomial likeihood [if it is NA use Poisson likelihood insted])
  if(!is.na(parameters$disp)){
    proposed.parameters$disp<-generate.proposal.disp(parameters$disp) # proposed update
    current.prior.density<-disp.prior.density(parameters$disp)
    proposal.prior.density<-disp.prior.density(proposed.parameters$disp)
    if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
      proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE) 
      acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
      if(runif(1)<acceptance.prob){  
        parameters <- proposed.parameters  
        current.loglikelihood<- proposal.loglikelihood
        disp.moves.accepted<-disp.moves.accepted+1
      }  else proposed.parameters<-parameters # otherwise current parameters stay the same
    } else  proposed.parameters<-parameters
  }
  
  
  
  # now store results of this iteration if i is a multiple of the parameter "thin"
  if(i%% thin ==0) {
    new.posterior.sample<-as.data.frame(matrix(c(as.numeric(parameters),current.loglikelihood$LL),nrow=1))
    colnames(new.posterior.sample)<-c(names(as.vector(parameters)),"LL")
    posterior.samples<-rbind(posterior.samples, new.posterior.sample)
    print(c(i,c(as.numeric(parameters),current.loglikelihood$LL)))
  }
} #  end of main MCMC loops



#  Summzarize acceptance probs
print(c("Accept ratio for beta1", beta1moves.accepted/N))
print(c("Accept ratio for beta2", beta2moves.accepted/N))
print(c("Accept ratio for beta3", beta3moves.accepted/N))
#print(c("Accept ratio for R01", R01moves.accepted/N))
#print(c("Accept ratio for R02", R02moves.accepted/N))
#print(c("Accept ratio for R03", R03moves.accepted/N))
print(c("Accept ratio for gamma", gamma.moves.accepted/N))
print(c("Accept ratio for rho", rho.moves.accepted/N))
print(c("Accept ratio for time.infected 1", time.infected1.moves.accepted/N))
print(c("Accept ratio for time.infected 2", time.infected2.moves.accepted/N))
print(c("Accept ratio for time.infected 3", time.infected3.moves.accepted/N))
print(c("Accept ratio for prop symptomatic 3", propSymptoms.moves.accepted/N))
print(c("Accept ratio for disp", disp.moves.accepted/N))
print(c("Accept ratio for w1", w1.moves.accepted/N))
print(c("Accept ratio for w2", w2.moves.accepted/N))
print(c("Accept ratio for w3", w3.moves.accepted/N))
#save.image()
save.image("model4 SEIIR  same infectiousness no watsan MH1.9 inf priors 3M thin 100.RData")

# plot selected runs 
pdf("model4 SEIIR  - selected fits from 3M (thin 100) chain.pdf")
par(mfrow=c(3,2),mar=c(4,4,1,1))
parameters<-posterior.samples[1,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

#
parameters<-posterior.samples[3000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

parameters<-posterior.samples[6000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

parameters<-posterior.samples[9000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

parameters<-posterior.samples[12000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)

parameters<-posterior.samples[15000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIIR=TRUE)
# 
#parameters<-posterior.samples[10000,]
#getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE)
# 
#parameters<-posterior.samples[15000,]
#getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE)
# 
#parameters<-posterior.samples[25000,]
#getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE)
# 
dev.off()

#posterior.samples<-rbind(modv1.5_1.5mrunsthin100.NBLLflatpriors.posterior.samples[,1:11],posterior.samples[,1:11])

# plot chains
pdf("model4 SEIIR same infectiousness no watsan MH1.9 inf priors 3M thin 100.pdf")
par(mfrow=c(3,4),mar=c(4,4,1,1))
rows.to.plot<-seq(1,length(posterior.samples$R01),1)
plot(posterior.samples$R01[rows.to.plot],type='l',main="R01")
plot(posterior.samples$R02[rows.to.plot],type='l',main="R02")
plot(posterior.samples$R03[rows.to.plot],type='l',main="R03")
plot(posterior.samples$time.first.infected1[rows.to.plot],type='l',main="Time first infection 1")
plot(posterior.samples$time.first.infected2[rows.to.plot],type='l',main="Time first infection 2")
plot(posterior.samples$time.first.infected3[rows.to.plot],type='l',main="Time first infection 3")
plot(posterior.samples$gamma[rows.to.plot],type='l',main="gamma")
plot(posterior.samples$rho[rows.to.plot],type='l',main="rho")
plot(posterior.samples$proportion.symptomatic[rows.to.plot],type='l',main="Proportion symptomatic")
#plot(posterior.samples$w1[rows.to.plot],type='l',main="watsan camp 1")
#plot(posterior.samples$w2[rows.to.plot],type='l',main="watsan camp 2")
#plot(posterior.samples$w3[rows.to.plot],type='l',main="watsan camp 3")

dev.off()

# plot correlation
pdf("model4 SEIIR  inf prior  NBLL 3M correlation.pdf")
par(mfrow=c(3,3),mar=c(4,4,1,1))
rows.to.plot<-seq(1,length(posterior.samples$R01),1)
plot(posterior.samples$R01[rows.to.plot],posterior.samples$gamma[rows.to.plot],xlab="R01",ylab="gamma",type='p',main="R01 v gamma",pch=".")
plot(posterior.samples$R01[rows.to.plot],posterior.samples$rho[rows.to.plot],xlab="R01",ylab="rho",type='p',main="R01 v rho",pch=".")
plot(posterior.samples$R01[rows.to.plot]*posterior.samples$rho[rows.to.plot],posterior.samples$rho[rows.to.plot],xlab="beta1",ylab="rho",type='p',main="beta1 v rho",pch=".")
plot(posterior.samples$R02[rows.to.plot]*posterior.samples$rho[rows.to.plot],posterior.samples$rho[rows.to.plot],xlab="beta2",ylab="rho",type='p',main="beta2 v rho",pch=".")
plot(posterior.samples$R03[rows.to.plot]*posterior.samples$rho[rows.to.plot],posterior.samples$rho[rows.to.plot],xlab="beta3",ylab="rho",type='p',main="beta3 v rho",pch=".")

plot(posterior.samples$R01[rows.to.plot],posterior.samples$time.first.infected1[rows.to.plot],xlab="R01",ylab="time.first.infected1",type='p',main="R01 v time.first.infected1",pch=".")
plot(posterior.samples$rho[rows.to.plot],posterior.samples$gamma[rows.to.plot],xlab="rho",ylab="gamma",type='p',main="rho v gamma",pch=".")
plot(posterior.samples$rho[rows.to.plot],posterior.samples$time.first.infected1[rows.to.plot],xlab="rho",ylab="time.first.infected1",type='p',main="rho v time.first.infected1",pch=".")
dev.off()

# calc stats  - mean and 95% CIs for key quantities
# first chuck out first 1000 samples (i.e. 100,000 iterations) as a burnin 
posterior.samples<-posterior.samples[1001: dim(posterior.samples)[1],]


# calc DIC
range<-10000:29000
dbar<- -2*mean(posterior.samples$LL[range])
#post.means<-data.frame(t(mean(posterior.samples)))
post.means<-data.frame(t(apply(posterior.samples[range,],2,mean)))[1:13]

# 1. Calculate DIC using the Spiegelhalter approach
dhat<-  -2*getLLHepEmod.mh.combined.no.watsan(post.means, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=TRUE)$LL

pD<-dbar-dhat
DIC<- dbar+pD
dbar
DIC
pD

# 2. Calculate DIC using the Gelman approach (invariant to reparameterisation so might be preferred)

p_v<- var(-2*posterior.samples$LL[range])
DIC2<- dbar+ p_v
DIC2


# 1. latent period

# mean

mean(365/posterior.samples$gamma)
quantile( 365/posterior.samples$gamma, c(.025,.5,.975))

# 2. infectious period

# mean

mean(2* 365/posterior.samples$rho)
quantile(2* 365/posterior.samples$rho, c(.025,.5,.975))

# 3.proportion of infections which are symptomatic and reported
mean(posterior.samples$proportion.symptomatic)
quantile(posterior.samples$proportion.symptomatic, c(.025,.5,.975))


# 4.R0s for the three camps
mean(posterior.samples$R01)
quantile(posterior.samples$R01, c(.025,.5,.975))
mean(posterior.samples$R02)
quantile(posterior.samples$R02, c(.025,.5,.975))
mean(posterior.samples$R03)
quantile(posterior.samples$R03, c(.025,.5,.975))


# 5.Time of first infection for the three camps

365*mean(posterior.samples$time.first.infected1)
365*quantile(posterior.samples$time.first.infected1, c(.025,.5,.975))

365*mean(posterior.samples$time.first.infected2)
365*quantile(posterior.samples$time.first.infected2, c(.025,.5,.975))

365*mean(posterior.samples$time.first.infected3)
365*quantile(posterior.samples$time.first.infected3, c(.025,.5,.975))




# 6.Watsan effect (water per person) on transmission parameter (-ve means a reduction in transmission with more wpp)

mean(posterior.samples$w1)
quantile(posterior.samples$w1, c(.025,.5,.975))

mean(posterior.samples$w2)
quantile(posterior.samples$w2, c(.025,.5,.975))

mean(posterior.samples$w3)
quantile(posterior.samples$w3, c(.025,.5,.975))


# 7. dispersion parameter of negbin dist

mean(posterior.samples$disp)
quantile(posterior.samples$disp, c(.025,.5,.975))
# Plots of predicted v observed at each of the camps - start at week 27 2007 finish at week 26 2009

# First sample P runs for plotting randomly selected without replacement
P<-100
n.samples<-length(posterior.samples[,1])
runs.to.plot<-  c(1:n.samples)[order(runif(n.samples))][1:P]
predicted1<-NULL
predicted2<-NULL
predicted3<-NULL
q.025prob.observed1<-NULL
q.975prob.observed1<-NULL
q.025prob.observed2<-NULL
q.975prob.observed2<-NULL
q.025prob.observed3<-NULL
q.975prob.observed3<-NULL

for(i in runs.to.plot){
  weeks.from.infection.symptoms<-1
  parameters<-posterior.samples[i,]
  parameters1<-c(
    beta.max=parameters$R01 *parameters$rho,
    beta.min=parameters$R01 *parameters$rho, #assume no seasonality  
    phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
    gamma1=parameters$gamma,
    gamma2=parameters$gamma,
    rho=parameters$rho
    #  r.water=0.000001,
    # r.latrines=0.000001
    )
  parameters2<-c(
    beta.max=parameters$R02 *parameters$rho,
    beta.min=parameters$R02 *parameters$rho, #assume no seasonality  
    phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
    gamma1=parameters$gamma,
    gamma2=parameters$gamma,
    rho=parameters$rho
    # r.water=0.000001,
    # r.latrines=0.000001
    )
  parameters3<-c(
    beta.max=parameters$R03 *parameters$rho,
    beta.min=parameters$R03 *parameters$rho, #assume no seasonality  
    phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs)
    gamma1=parameters$gamma,
    gamma2=parameters$gamma,
    rho=parameters$rho
    # r.water=0.000001,
    # r.latrines=0.000001
    )  
  
  time.first.infected1<-parameters$time.first.infected1
  time.first.infected2<-parameters$time.first.infected2
  time.first.infected3<-parameters$time.first.infected3
  
  times1<- seq(ceiling(time.first.infected1*52)/52, 157/52, by = wk)
  times2<- seq(ceiling(time.first.infected2*52)/52, 157/52, by = wk)
  times3<- seq(ceiling(time.first.infected3*52)/52, 157/52, by = wk)
  
  
  if(time.first.infected1!=times1[1]) times1<-c(time.first.infected1,times1)
  if(time.first.infected2!=times2[1]) times2<-c(time.first.infected2,times2)
  if(time.first.infected3!=times3[1]) times3<-c(time.first.infected3,times3)
  
  # below should be same as out1 out2 and out3 calculated using the compiled code (and they ares)
  #out1a<-ode(y=xstart1,times=times1, func=HepEmod.nowatsan,parameters1)
  #out2a<-ode(y=xstart2,times=times2, func=HepEmod.nowatsan,parameters2)
  #out3a<-ode(y=xstart3,times=times3, func=HepEmod.nowatsan,parameters3)
  
  params1<-c(beta=parameters1[1],gamma=parameters1[4], rho=parameters1[5])
  params2<-c(beta=parameters2[1],gamma=parameters2[4], rho=parameters2[5])
  params3<-c(beta=parameters3[1],gamma=parameters3[4], rho=parameters3[5])
  
  # out1<-ode(y=xstart1,times=times1, func="derivs",parms=params1,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
  # out2<-ode(y=xstart2,times=times2, func="derivs",parms=params2,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
  #out3<-ode(y=xstart3,times=times3, func="derivs",parms=params3,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
  
  
  out1<-ode(y=xstart1,times=times1, func="derivs",parms=params1,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum")
  out2<-ode(y=xstart2,times=times2, func="derivs",parms=params2,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum")
  out3<-ode(y=xstart3,times=times3, func="derivs",parms=params3,jacfun="jac",dllname="SEEIR", initfunc="initmod",nout=1,outnames="Sum")
  
  
  
  cuminc1<-out1[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
  cuminc2<-out2[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
  cuminc3<-out3[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
  
  obs1<-OBSERVED.CASES1[27:130] # data from week 27 2007 to week 26 2009 
  obs2<-OBSERVED.CASES2[27:130] # data from week 27 2007 to week 26 2009 
  obs3<-OBSERVED.CASES3[27:130] # data from week 27 2007 to week 26 2009 
  
  n1<-length(cuminc1)
  n2<-length(cuminc2)
  n3<-length(cuminc3)
  
  # first remove weeks 131 to 157 from cuminc i.e. the last 27 ob
  cuminc1<-cuminc1[1:(n1-27)]
  cuminc2<-cuminc2[1:(n2-27)]
  cuminc3<-cuminc3[1:(n3-27)]
  
  n1<-length(cuminc1)
  n2<-length(cuminc2)
  n3<-length(cuminc3)
  
  weeklyinc1<-c(cuminc1[1],cuminc1[2:n1]-cuminc1[1:(n1-1)]) # weekly incidence
  weeklyinc2<-c(cuminc2[1],cuminc2[2:n2]-cuminc2[1:(n2-1)]) # weekly incidence
  weeklyinc3<-c(cuminc3[1],cuminc3[2:n3]-cuminc3[1:(n3-1)]) # weekly incidence
  
  n1<-length(cuminc1)
  n2<-length(cuminc2)
  n3<-length(cuminc3)
  
  wk.inc1<-weeklyinc1[1:(n1-weeks.from.infection.symptoms)] # this gives total incident infectious caes - a proportion of whcih will be symptomatic & reported after a delay of weeks.from.infection.symptoms
  wk.inc2<-weeklyinc2[1:(n2-weeks.from.infection.symptoms)] # this gives total incident infectious caes - a proportion of whcih will be symptomatic & reported after a delay of weeks.from.infection.symptoms
  wk.inc3<-weeklyinc3[1:(n3-weeks.from.infection.symptoms)] # this gives total incident infectious caes - a proportion of whcih will be symptomatic & reported after a delay of weeks.from.infection.symptoms
  
  n1<-length(wk.inc1) 
  n2<-length(wk.inc2) 
  n3<-length(wk.inc3) 
  
  
  firstwk<-ceiling(time.first.infected1*52) + weeks.from.infection.symptoms # first week where we have cuminc data
  if(firstwk>27){
    wk.inc1 <- c(rep(0,104-n1), wk.inc1)
  }
  if(firstwk<27){
    wk.inc1 <- wk.inc1[(28-firstwk):n1]   
  }
  
  firstwk<-ceiling(time.first.infected2*52) + weeks.from.infection.symptoms # first week where we have cuminc data
  if(firstwk>27){
    wk.inc2 <- c(rep(0,104-n2), wk.inc2)
  }
  if(firstwk<27){
    wk.inc2 <- wk.inc2[(28-firstwk):n2]   
  }
  
  firstwk<-ceiling(time.first.infected3*52) + weeks.from.infection.symptoms # first week where we have cuminc data
  if(firstwk>27){
    wk.inc3 <- c(rep(0,104-n3), wk.inc3)
  }
  if(firstwk<27){
    wk.inc3 <- wk.inc3[(28-firstwk):n3]   
  }
  
  predicted1<- rbind(predicted1,parameters$proportion.symptomatic*wk.inc1)
  predicted2<- rbind(predicted2,parameters$proportion.symptomatic*wk.inc2)
  predicted3<- rbind(predicted3,parameters$proportion.symptomatic*wk.inc3)

  # negative binomial likelihood
  
  q.025prob.observed1<-rbind(q.025prob.observed1,qnbinom(.025,mu=round(wk.inc1)*parameters$proportion.symptomatic,size=parameters$disp))
  q.975prob.observed1<-rbind(q.975prob.observed1,qnbinom(.975,mu=round(wk.inc1)*parameters$proportion.symptomatic,size=parameters$disp))
  q.025prob.observed2<-rbind(q.025prob.observed2,qnbinom(.025,mu=round(wk.inc2)*parameters$proportion.symptomatic,size=parameters$disp))
  q.975prob.observed2<-rbind(q.975prob.observed2,qnbinom(.975,mu=round(wk.inc2)*parameters$proportion.symptomatic,size=parameters$disp))
  q.025prob.observed3<-rbind(q.025prob.observed3,qnbinom(.025,mu=round(wk.inc3)*parameters$proportion.symptomatic,size=parameters$disp))
  q.975prob.observed3<-rbind(q.975prob.observed3,qnbinom(.975,mu=round(wk.inc3)*parameters$proportion.symptomatic,size=parameters$disp))
  
  # poisson likelihood
  
  #   q.025prob.observed1<-rbind(q.025prob.observed1,qpois(.025,lambda=round(wk.inc1)*parameters$proportion.symptomatic))
  #   q.975prob.observed1<-rbind(q.975prob.observed1,qpois(.975,lambda=round(wk.inc1)*parameters$proportion.symptomatic))
  #   q.025prob.observed2<-rbind(q.025prob.observed2,qpois(.025,lambda=round(wk.inc2)*parameters$proportion.symptomatic))
  #   q.975prob.observed2<-rbind(q.975prob.observed2,qpois(.975,lambda=round(wk.inc2)*parameters$proportion.symptomatic))
  #   q.025prob.observed3<-rbind(q.025prob.observed3,qpois(.025,lambda=round(wk.inc3)*parameters$proportion.symptomatic))
  #   q.975prob.observed3<-rbind(q.975prob.observed3,qpois(.975,lambda=round(wk.inc3)*parameters$proportion.symptomatic))
  #   
}
mean1<-apply(predicted1, 2, mean)
mean2<-apply(predicted2, 2, mean)
mean3<-apply(predicted3, 2, mean)

quants1<-apply(predicted1, 2, quantile,c(.025,.5,.975))
quants.ob.lower1<-apply(q.025prob.observed1, 2, mean)
quants.ob.upper1<-apply(q.975prob.observed1, 2, mean)

quants2<-apply(predicted2, 2, quantile,c(.025,.5,.975))
quants.ob.lower2<-apply(q.025prob.observed2, 2, mean)
quants.ob.upper2<-apply(q.975prob.observed2, 2, mean)

quants3<-apply(predicted3, 2, quantile,c(.025,.5,.975))
quants.ob.lower3<-apply(q.025prob.observed3, 2, mean)
quants.ob.upper3<-apply(q.975prob.observed3, 2, mean)

pdf("model4 SEIIR  Predicted versus observed NB 3m thin 100.pdf",height=6,width=5,pointsize=10)
par(mfrow=c(3,1),mar=c(4,4,2,4))

col1<-"red"
col1a<-"pink"
col1b<-"grey90"
col2<-"palegreen4"
col2a<-"palegreen1"
col3<-"blue"
col3a<-"lightcyan"

# camp 1
plot(27:130, obs1, type='p',col=col1,ylim=c(0,250),xlab="", ylab="", main="Agoro")
# need to think what we want to plot below..
polygon(c(27:130,130:27),c(quants.ob.lower1,rev(quants.ob.upper1)),col=col1b,border=col1b)

polygon(c(27:130,130:27),c(quants1[1,],rev(quants1[3,])),col=col1a,border=col1a)
lines(27:130, mean1,col=col1)
points(27:130, obs1,col=col1)

# camp 2
plot(27:130, obs2, type='p',col=col2,ylim=c(0,200),xlab=" ", ylab="Cases", main="MadiOpei")

polygon(c(27:130,130:27),c(quants.ob.lower2,rev(quants.ob.upper2)),col=col1b,border=col1b)

polygon(c(27:130,130:27),c(quants2[1,],rev(quants2[3,])),col=col2a,border=col2a)
lines(27:130, mean2,col=col2)
points(27:130, obs2,col=col2)

# camp 3
plot(27:130, obs3, type='p',col=col3,ylim=c(0,200),xlab="Week number", ylab="", main="Paloga")

polygon(c(27:130,130:27),c(quants.ob.lower3,rev(quants.ob.upper3)),col=col1b,border=col1b)

polygon(c(27:130,130:27),c(quants3[1,],rev(quants3[3,])),col=col3a,border=col3a)
lines(27:130, mean3,col=col3)
points(27:130, obs3,col=col3)
dev.off()

