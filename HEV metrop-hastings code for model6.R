
# HEV metrop-hastings code for model6
#( SEIWR model with moderately informative priors allowing transmisison via shared and saturating water sources)
# Ben Cooper 2016

# Note that code here is set up for prior pi1 e1. Just edit these lines for different prior assumptions.
# see #priors for SEIRW model  and definitions for pi.prior.density and mu.prior.density
# Fit the model to all three camps together using Metropolis-Hastings Algorithm
# estimating common natural history parameters and the proportion symptomaticc
# Model 6 does this with  an SEIWR where W corresponds to contaminated water

# camp specific R0s, camp-specific time of first case but estimating 
# gamma and proportion sympomatic which are the same for all three camps.
# Assume data follow Negative binomial distribution around predicted means 
# and use compiled C code to run the ODE model 
# ************************************************************************************************



library(deSolve)

# first compiled and load the C code for the SEIR model

system("R CMD SHLIB SEEIR.c")
dyn.load("SEEIR.so")

system("R CMD SHLIB SEIIR.c")
dyn.load("SEIIR.so")

system("R CMD SHLIB SEIWR.c")
dyn.load("SEIWR.so")

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
    # Return values 
    list(c(dS,dE,dI,dR,dCumInc))
    
  })  # end with
} # end HepEmod



getLLHepEmod.mh.combined.no.watsan<-function(params, x1,x2,x3, plot=FALSE,verbose=FALSE, SEEIR=FALSE, SEIIR=FALSE, SEIRW=TRUE){
  #The version is designed for use with Metropolis Hastings algorithm
  # assumes no effect of watsan
  # x1, x2, x3 are vectors of state variables (S, E, I, R and possibly W  for each camp at start time (assume this vector is known )
  # if SEEIR ==TRUE thihs uses the dll for the SEEIR model which assumes two equal comparments (by default with same rate of leaving
  # so latent period has an erlang distribution) 
  # if SEIIR ==TRUE thihs uses the dll for the SEEIR model which assumes two equal comparments (by default with same rate of leaving
  # so latent period has an erlang distribution) 
  # if SEIRW == TRUE this uses the DLL for the SEIRW model accounting for contaminated Water 
 
  
  with(as.list(c(params, x1,x2,x3)), {
    # R01 is R0  in camp 1 
    # R02 is R0  in camp 2 
    # R03 is R0  in camp 3 
    # time.first.infected1  is time of first case in year 1 in camp 1 (scale is years)
    # time.first.infected2  is time of first case in year 1 in camp 2 (scale is years)
    # time.first.infected3  is time of first case in year 1 in camp 3 (scale is years)
    # gamma is rate of progression to exposed to infectious
    # rho is rate of progression from infectious to recovered (log(1/infectious period)
    # lambda is rate of shedding from I to W
    # mu is rate of decay from W (i.e. rate contamination is lost)
    
    # proportion.symptomatic is  proporotion of cases seen which are symptomatic
    # proportion.of.transmission.early.in.infection - is proportion of total  transmission due to the first I compartment
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
    
   
      
    if(SEEIR==FALSE & SEIIR==FALSE & SEIRW == FALSE){
     out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
     out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum")
     out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
    } else if(SEEIR==TRUE & SEIIR==FALSE & SEIRW==FALSE) {
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
   
     
    } else if(SEEIR==FALSE & SEIIR==TRUE &  SEIRW == FALSE) { # in this case set parameter for model where beta1 =beta2 and rho1 = rho2
        parameters1<-c(
        beta1=proportion.of.transmission.early.in.infection*R01 *rho,
        beta2=(1-proportion.of.transmission.early.in.infection)*R01 *rho,
        gamma=gamma,
        rho1=rho,
        rho2=rho)
      
      parameters2<-c(
        beta1=proportion.of.transmission.early.in.infection*R02 *rho,
        beta2=(1-proportion.of.transmission.early.in.infection)*R02 *rho, 
        gamma=gamma,
        rho1=rho,
        rho2=rho)
      
      parameters3<-c(
        beta1=proportion.of.transmission.early.in.infection*R03 *rho,
        beta2=(1-proportion.of.transmission.early.in.infection)*R03 *rho, 
        gamma=gamma,
        rho1=rho,
        rho2=rho)
        
      out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEIIR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
    } else if(SEEIR==FALSE & SEIIR==FALSE  &  SEIRW == TRUE) { 
    	#env contamination model
    	# where R01D is R0 component in pop 1 due to direct transmission
    	# and RO1E is R0 component in pop 1 due to indirect (environmental) transmission etc 

      parameters1<-c(
        betaD=R01D *rho,  # direct transmission 
        gamma=gamma,
        rho=rho,
        lambda=lambda,
        mu=mu,
        betaE=R01E*rho*mu/lambda)
      
      parameters2<-c(
        betaD=R02D *rho,  # direct transmission 
        gamma=gamma,
        rho=rho,
        lambda=lambda,
        mu=mu,
        betaE=R02E*rho*mu/lambda)
      
      
      parameters3<-c(
        betaD=R03D *rho,  # direct transmission 
        gamma=gamma,
        rho=rho,
        lambda=lambda,
        mu=mu,
        betaE=R03E*rho*mu/lambda)    
 
      out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIWR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out2<-ode(y=x2,times=times2, func="derivs",parms=parameters2,jacfun="jac",dllname="SEIWR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      out3<-ode(y=x3,times=times3, func="derivs",parms=parameters3,jacfun="jac",dllname="SEIWR", initfunc="initmod",nout=1,outnames="Sum",maxsteps=1000000)
      
      
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
    
    # now use weeklyinc to construct a likelihood  
    # Could use a binomial likelihood noting that dbinom returns a 0 if number of trials is fewer than 
    # number successes (so returned log likelihood will be -Inf)
    # below OBSERVED.CASES is the actual data (specified as global variable)
    # and the likelihood is calculated starting at LATEST.START.WEEK (last time virus could be introduced into camp)
    # However using a poisson or negbin likelihood is more robust (and negbin is likely to fit much better)
    
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
    #dbinomLLterms1<-dbinom(obs1,pmax(obs1,round(wk.inc1)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms2<-dbinom(obs2,pmax(obs2,round(wk.inc2)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms3<-dbinom(obs3,pmax(obs3,round(wk.inc3)), proportion.symptomatic,log=TRUE)
    # pmax term used because trials should be >= num successes othereise LL is -Inf
    # pmax term is there to allow us to find good initial values  - should then check that converged
    # values have weekly inc greater than observed (otherwise not feasible)
    if(verbose){
      print(c(numbertochop1,numbertochop2,numbertochop3))
      print(obs1/wk.inc1)  # should be less than one for a feasible value
      print(obs2/wk.inc2)
      print(obs3/wk.inc3)
    }  
    #binomLL<- sum(dbinomLLterms1) +sum(dbinomLLterms2) +sum(dbinomLLterms3)  
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
    # n.trials=sum(wk.inc1)+sum(wk.inc2)+sum(wk.inc3) # number of binomial trials (needed for Gibss updates)
    #n.successes=sum(obs1)+sum(obs2)+sum(obs3)
    #  return(list(LL=binomLL,n.trials=n.trials,n.successes=n.successes))
    return(list(LL=LL))
  })   
}



getLLHepEmod.mh.combined.stepwatsan<-function(params, x1,x2,x3, watsantimes,plot=FALSE,verbose=FALSE){
  #The version is designed for use with Metropolis Hastings algorithm
  # does not currently support SEEIR or SEIIR models
  # Treats watsan intervention as a step function, taking 0 or 1 values
  # x1, x2, x3 are vectors of state variables (S, E, I, R) for each camp at start time (assume this vector is known )
  # watsantimes is a vector hodling the time of the watsan interventions at the three camps (in units of years)
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
    # w1 is the watsan effect in camp1 where beta=beta0*exp(w1*watsan)  where watsan is a measure of watsan intervention 
    # w2 is the watsan effect in camp2 
    # w3 is the watsan effect in camp3 
    
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
    forcings1<-cbind(times1,as.integer(times1>=watsantimes[1]))
    forcings2<-cbind(times2,as.integer(times2>=watsantimes[2]))
    forcings3<-cbind(times3,as.integer(times3>=watsantimes[3]))
    
    
    parameters1<-c(
      beta=R01 *rho,
      gamma=gamma,
      rho=rho,
      w=w1)
    
    parameters2<-c(
      beta=R02 *rho,
      gamma=gamma,
      rho=rho,
      w=w2)
    
    parameters3<-c(
      beta=R03 *rho,
      gamma=gamma,
      rho=rho,
      w=w3)
    
    out1<-ode(y=x1,times=times1, func="derivs",parms=as.numeric(parameters1),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings1, initfunc="initmod",nout=1,outnames="Sum")
    out2<-ode(y=x2,times=times2, func="derivs",parms=as.numeric(parameters2),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings2, initfunc="initmod",nout=1,outnames="Sum")
    out3<-ode(y=x3,times=times3, func="derivs",parms=as.numeric(parameters3),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings3, initfunc="initmod",nout=1,outnames="Sum",maxsteps=10000)
    
    
    cuminc1<-out1[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
    cuminc2<-out2[,6][-1]
    cuminc3<-out3[,6][-1]
    
    n1<-length(cuminc1)
    n2<-length(cuminc2)
    n3<-length(cuminc3)
    
    weeklyinc1<-c(cuminc1[1],cuminc1[2:n1]-cuminc1[1:(n1-1)]) # weekly incidence
    weeklyinc2<-c(cuminc2[1],cuminc2[2:n2]-cuminc2[1:(n2-1)]) # weekly incidence
    weeklyinc3<-c(cuminc3[1],cuminc3[2:n3]-cuminc3[1:(n3-1)]) # weekly incidence
    
    # now use weeklyinc to construct a likelihood  
    # Use a binomial likelihood noting that dbinom returns a 0 if number of trials is fewer than 
    # number successes (so returned log likelihood will be -Inf)
    # below OBSERVED.CASES is the actual data (specified as global variable)
    # and the likelihood is calculated starting at LATEST.START.WEEK (last time virus could be introduced into camp)
    
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
    #dbinomLLterms1<-dbinom(obs1,pmax(obs1,round(wk.inc1)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms2<-dbinom(obs2,pmax(obs2,round(wk.inc2)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms3<-dbinom(obs3,pmax(obs3,round(wk.inc3)), proportion.symptomatic,log=TRUE)
    # pmax term used because trials should be >= num successes othereise LL is -Inf
    # pmax term is there to allow us to find good initial values  - should than check that converged
    # values have weekly inc greater than observed (otherwise not feasibl)
    if(verbose){
      print(c(numbertochop1,numbertochop2,numbertochop3))
      print(obs1/wk.inc1)  # should be less than one for a feasible value
      print(obs2/wk.inc2)
      print(obs3/wk.inc3)
    }  
    #binomLL<- sum(dbinomLLterms1) +sum(dbinomLLterms2) +sum(dbinomLLterms3)  
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
    # n.trials=sum(wk.inc1)+sum(wk.inc2)+sum(wk.inc3) # number of binomial trials (needed for Gibss updates)
    #n.successes=sum(obs1)+sum(obs2)+sum(obs3)
    #  return(list(LL=binomLL,n.trials=n.trials,n.successes=n.successes))
    return(list(LL=LL))
  })   
}

getLLHepEmod.mh.combined.wpp.watsan<-function(params, x1,x2,x3, watsan.wpp,plot=FALSE,verbose=FALSE){
  #The version is designed for use with Metropolis Hastings algorithm
  # does not currently support SEEIR or SEIIR models
  # Treats watsan intervention through water person data (wpp) which is interpolated 
  # x1, x2, x3 are vectors of state variables (S, E, I, R) for each camp at start time (assume this vector is known )
  # watsantimes is a vector hodling the time of the watsan interventions at the three camps (in units of years)
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
    # w1 is the watsan effect in camp1 where beta=beta0*exp(w1*watsan)  where watsan is a measure of watsan intervention 
    # w2 is the watsan effect in camp2 
    # w3 is the watsan effect in camp3 
    
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
    wpptimes<-seq(1/52, length(watsan.wpp[,1])/52, 1/52)
    wpp1<-approx(wpptimes, watsan.wpp[,2], xout=times1,rule = 2)$y
    wpp2<-approx(wpptimes, watsan.wpp[,3], xout=times2,rule = 2)$y
    wpp3<-approx(wpptimes, watsan.wpp[,4], xout=times3,rule = 2)$y
    
    forcings1<-cbind(times1,wpp1)
    forcings2<-cbind(times2,wpp2)
    forcings3<-cbind(times3,wpp3)
    
 
    parameters1<-c(
      beta=R01 *rho,
      gamma=gamma,
      rho=rho,
      w=w1)
    
    parameters2<-c(
      beta=R02 *rho,
      gamma=gamma,
      rho=rho,
      w=w2)
    
    parameters3<-c(
      beta=R03 *rho,
      gamma=gamma,
      rho=rho,
      w=w3)
    
    out1<-ode(y=x1,times=times1, func="derivs",parms=as.numeric(parameters1),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings1, initfunc="initmod",nout=1,outnames="Sum")
    out2<-ode(y=x2,times=times2, func="derivs",parms=as.numeric(parameters2),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings2, initfunc="initmod",nout=1,outnames="Sum")
    out3<-ode(y=x3,times=times3, func="derivs",parms=as.numeric(parameters3),jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings3, initfunc="initmod",nout=1,outnames="Sum")
    
    
    cuminc1<-out1[,6][-1] # the [-1] drops the first time point for cuminc, so we start counting at exact times of each new week
    cuminc2<-out2[,6][-1]
    cuminc3<-out3[,6][-1]
    
    n1<-length(cuminc1)
    n2<-length(cuminc2)
    n3<-length(cuminc3)
    
    weeklyinc1<-c(cuminc1[1],cuminc1[2:n1]-cuminc1[1:(n1-1)]) # weekly incidence
    weeklyinc2<-c(cuminc2[1],cuminc2[2:n2]-cuminc2[1:(n2-1)]) # weekly incidence
    weeklyinc3<-c(cuminc3[1],cuminc3[2:n3]-cuminc3[1:(n3-1)]) # weekly incidence
    
    # now use weeklyinc to construct a likelihood  
    # Use a binomial likelihood noting that dbinom returns a 0 if number of trials is fewer than 
    # number successes (so returned log likelihood will be -Inf)
    # below OBSERVED.CASES is the actual data (specified as global variable)
    # and the likelihood is calculated starting at LATEST.START.WEEK (last time virus could be introduced into camp)
    
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
    #dbinomLLterms1<-dbinom(obs1,pmax(obs1,round(wk.inc1)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms2<-dbinom(obs2,pmax(obs2,round(wk.inc2)), proportion.symptomatic,log=TRUE)
    #dbinomLLterms3<-dbinom(obs3,pmax(obs3,round(wk.inc3)), proportion.symptomatic,log=TRUE)
    # pmax term used because trials should be >= num successes othereise LL is -Inf
    # pmax term is there to allow us to find good initial values  - should than check that converged
    # values have weekly inc greater than observed (otherwise not feasibl)
    if(verbose){
      print(c(numbertochop1,numbertochop2,numbertochop3))
      print(obs1/wk.inc1)  # should be less than one for a feasible value
      print(obs2/wk.inc2)
      print(obs3/wk.inc3)
    }  
    #binomLL<- sum(dbinomLLterms1) +sum(dbinomLLterms2) +sum(dbinomLLterms3)  
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
    # n.trials=sum(wk.inc1)+sum(wk.inc2)+sum(wk.inc3) # number of binomial trials (needed for Gibss updates)
    #n.successes=sum(obs1)+sum(obs2)+sum(obs3)
    #  return(list(LL=binomLL,n.trials=n.trials,n.successes=n.successes))
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
#fixed
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
W0<-0
CumInc0<-0
xstart1<-c(S=N1-1,E=E0,I=I0, R=R0,CumInc=CumInc0,W=W0 )
xstart2<-c(S=N2-1,E=E0,I=I0, R=R0,CumInc=CumInc0,W=W0)
xstart3<-c(S=N3-1,E=E0,I=I0, R=R0,CumInc=CumInc0,W=W0)

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

# 
# # Uninformative priors below
R0.prior.density<-function(x){ dunif(x, min=1, max=20)}  #note we can either parameterise for R0 or for beta, but not both...
#beta.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)} #as R0=beta/rho
beta.prior.density<-function(x){ dunif(x, 0.0001,500)} #as
# betaE.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
betaE.prior.density<-function(x){dunif(x, 0.0001,500)}
#gamma.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
gamma.prior.density<-function(x){dunif(x, 0.0001,500)}
#rho.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
rho.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
#lambda.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
lambda.prior.density<-function(x){dunif(x, 0.0001,1000)}
lambda.prior.density<-function(x){dunif(x, 0.0001,1)} # with lower values in this range saturation is quite limited, but with other upper values it does limit transmission


#mu.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
mu.prior.density<-function(x){dunif(x, 0.0001,500)}

time.first.infected1.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK1/52)}
time.first.infected2.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK2/52)}
time.first.infected3.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK3/52)}
propSymptoms.prior.beta.shape1<-1  # parmeters for propSymptoms.prior
propSymptoms.prior.beta.shape2<-1  # parmeters for propSymptoms.prior
propSymptoms.prior.density<-function(x){ dbeta(x, propSymptoms.prior.beta.shape1,propSymptoms.prior.beta.shape2)}
proportion.of.transmission.early.in.infection.density<-function(x){ dbeta(x, 1,1)} 
disp.prior.density<-function(x){ dunif(x, min=0.1, max=10000)}

#pi.prior.density<-function(x){dbeta(x, 1,1)}
#pi is the proportion of R0 due to environmental spread. 
#pi.prior.density<-function(x){dbeta(x, 81,20)}   # prior expresses quite strong belief that majority of transmission is through environment (about 80%)
#pi.prior.density<-function(x){dbeta(x, 21,80)}   # prior expresses quite strong belief that majority of transmission is direct(about 80%)
pi.prior.density<-function(x){dbeta(x, 1,1)}  #  flat prior

# So if a and b are beta priors and there are x successess from n trials
# from beta conjuage properties for binomial prob the posteriors are
#  a_post= a +x and b_post =b + n -x 



# Informative priors below
# 
# R0.prior.density<-function(x){ dunif(x, min=1, max=30)}  #note we can either parameterise for R0 or for beta, but not both...
# #beta.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)} #as R0=beta/rho
# beta.prior.density<-function(x){ dunif(x, 0.0001,500)} #as

# pi is the proportion of R0 due to environmental spread. 
#pi.prior.density<-function(x){dbeta(x, 81,20)}   # prior expresses quite strong belief that majority of transmission is through environment (about 80%)
#pi.prior.density<-function(x){dbeta(x, 21,80)}   # prior expresses quite strong belief that majority of transmission is direct(about 80%)

gamma.prior.density<-function(x){ dgamma(x, shape=1,rate=0.0932)} 
#rate is 34/365. Conjugate based on one observation with a mean latent period of 34 days
#Chauhan A, Jameel S, Dilawari JB, Chawla YK, Kaur U, Ganguly NK. Lancet 1993; 341: 149–50.

rho.prior.density<-function(x){ dgamma(x, shape=11,rate=1.145)}
 # above  based on Takahashi et al 2007: 11 patients with a mean duration of infectiousness of 38 days (assuming 1 week infectious before symptoms)
#  this comes from using a gamma(0,1)
#lambda.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}



# #above  based on Takahashi et al 2007: 11 patients with a mean duration of infectiousness of 38 days (assuming 1 week infectious before symptoms)
# #conjugate gamma shape and rate are shape=a+n,rate=r +n*x if prior has shape a and rate r
# #and if n observation with a mean of x. Time units are in years here
# time.first.infected1.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK1/52)}
# time.first.infected2.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK2/52)}
# time.first.infected3.prior.density<-function(x){ dunif(x, min=0, max=LATEST.START.WEEK3/52)}
# propSymptoms.prior.beta.shape1<-2  # parmeters for propSymptoms.prior (corresponds to week evidence that 1 in 8 positive)
# propSymptoms.prior.beta.shape2<-8  # parmeters for propSymptoms.prior
# propSymptoms.prior.density<-function(x){ dbeta(x, propSymptoms.prior.beta.shape1,propSymptoms.prior.beta.shape2)}
# disp.prior.density<-function(x){ dunif(x, min=0.1, max=10000)}
# watsan.effect.prior.density<-function(x){ dnorm(x, mean=0, sd=100)}


#priors for SEIRW model   
pi.prior.density<-function(x){dbeta(x, 0.03,3.013)}  #  (mean= 0.01, 95% interval=  5.202893e-54 1.247330e-01 )
mu.prior.density<-function(x){ dgamma(x, shape=5.245,scale=34.795)} # rate of decay from environment 
# corresponding to: e1  (mean= 365/2 , 95% CIs: 61.29262 368.65418 )  [so mean duration and 9%% is 2 (1,6) days ]



# 2. Define functions to propose new values for parameters
# Adjust variance of the proposals to tune the algorithm (ideally we want to accept about a fifth of proposals for each variable)
# Note that these functions may propose values outside possible ranges, but since prior distributions
# will have zero probability at these values they have no chance of being accepted
generate.proposal.R0<-function(currentR0){return(currentR0 + rnorm(1,0,.1))}
generate.proposal.beta<-function(currentbeta){return(currentbeta + rnorm(1,0,5))}
generate.proposal.gamma<-function(currentgamma){return(currentgamma + rnorm(1,0,.4))}
generate.proposal.rho<-function(currentrho){return(currentrho + rnorm(1,0,.4))}
generate.proposal.lambda<-function(currentlambda){return(currentlambda + rnorm(1,0,.05))}

generate.proposal.mu<-function(currentmu){return(currentmu + rnorm(1,0,5))}
generate.proposal.time.first.infected<-function(currenttime){return(currenttime + rnorm(1,0,.01))}
generate.proposal.propSymptoms<-function(currentprop){return(currentprop + rnorm(1,0,.01))} #   used if using Poisson likelihood (if binom likilihood use a Gibbs update for this one)
generate.proposal.proportion.of.transmission.early.in.infection <-function(currentprop){return(currentprop + rnorm(1,0,.01))} # only for SEIIR model
generate.proposal.disp<-function(currentdisp){return(currentdisp + rnorm(1,0,2))}
generate.proposal.watsan.effect<-function(currentw){return(currentw + rnorm(1,0,50))}
generate.proposal.pi<-function(currentpi){return(currentpi + rnorm(1,0,.05))}

# initial values
pi1<-.5
pi2<-.5
pi3<-.5


rho=20
lambda=0.1
gamma=10
mu=5
R01D<-4
R02D<-5
R03D<-6

betaD1=R01D*gamma
betaD2=R02D*gamma
betaD3=R03D*gamma

R01E<-lambda*pi1*betaD1/(rho * mu *(gamma * lambda - gamma*lambda*pi1))
R02E<-lambda*pi2*betaD2/(rho * mu *(gamma * lambda - gamma*lambda*pi2))
R03E<-lambda*pi3*betaD3/(rho * mu *(gamma * lambda - gamma*lambda*pi3))


parameters<-list(
  R01D=4,  # R0 component due to direct transmission
  R02D=5, 
  R03D=6,
  R01E=R01E, # R0 component due to evnironmental transmission
  R02E=R02E, 
  R03E=R03E, 
  gamma=10.4,
  rho=20,
  lambda=0.1,
  mu=5,
  proportion.of.transmission.early.in.infection=0.5, #only relevant in SEIIR model
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
# Note that we can parameterise the model different ways havong different "fundamental" parameters which we give priors too. 
# Depending on what we want to put priors on we can change these fundabmental prameters, but don't need to change this list "parameters"
# defined above - and it is convenient to keep it the same so we don't need to change the interface to the functions that calculate the likelihiood.
# However, if changing the parametersiation we do need to be careful when implementing the MCMC sampler to make sure when parameters xyz changes all
# parameters defined in terms of xyz change accordingly.
# In the current version (metrop-hastings1.12R) pi (proportion of transmission due to env route) is given a prior and considered to be fundamental
# while R01E etc is an implied parameters. 

  currentpi1<-R01E/(R01E+R01D)
  currentpi2<-R02E/(R02E+R02D)
  currentpi3<-R03E/(R03E+R03D)

#  parameters$disp<-NA # setting to NA forces a Poisson likelihood
parameters<-unlist(parameters)
parameters<-as.data.frame(t(parameters))

parameters[1:17]<-posterior.samples[50001,1:17]
# uncomment above line to start with same parameters left in the preivios run

#thin<-10

# with version 1.12 we can do about 75000 runs in one hour 
# So 1 million overnight should be fine

thin<-100

N<-5000000
last.posterior.samples<-posterior.samples

#parameters<-modv1.670to1670k.posterior.samples[10001,1:10]
m<-dim(last.posterior.samples)[1]
#parameters<-last.posterior.samples[m,1:11]
parameters<-last.posterior.samples[m,]
l<-length(parameters)
parameters<-parameters[1:(l-1)]
#parameters$w1<-0
#parameters$w2<-0
#parameters$w3<-0
#parameters$time.first.infected3<-0.9
#parameters$disp<-5
# parameters<-posterior.samples.150k.negbin[m,1:11]


proposed.parameters<-parameters


#parameters$disp<-1
#current.loglikelihood<-getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3,Uganda.wpp, plot=TRUE)
#parameters$disp<-1
#current.loglikelihood<-getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3,Uganda.wpp, plot=TRUE)


# Below the SEIIR option, if set to true, causes likelihood to be based on SEIIR model, where time is in two 
# I compartments is equal (giving an erlang distributed infectious period - though if using this should use different priors 
# for the rho parameters )
current.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW=TRUE)

#store results from the posterior sample with one column for each parameter & one for LL
posterior.samples<-as.data.frame(matrix(c(as.numeric(parameters),current.loglikelihood$LL),nrow=1))
colnames(posterior.samples)<-c(names(as.vector(parameters)),"LL")
betaD1moves.accepted<-0    # note that we either update R0s or betas but not both
betaD2moves.accepted<-0    # note that we either update R0s or betas but not both
betaD3moves.accepted<-0    # note that we either update R0s or betas but not both
pi1.moves.accepted<-0 # pi1 is R01E/(R01E+R01D) i.e. proportion of total R0 due to the environmental route  in population 1
pi2.moves.accepted<-0 
pi3.moves.accepted<-0 


#sbetaE1moves.accepted<-0    # note that we either update R0s or betas but not both
#sbetaE2moves.accepted<-0    # note that we either update R0s or betas but not both
#sbetaE3moves.accepted<-0    # note that we either update R0s or betas but not both
#sR01Dmoves.accepted<-0    # note that we either update R0s or betas but not both
#sR02Dmoves.accepted<-0    # note that we either update R0s or betas but not both
#sR03Dmoves.accepted<-0    # note that we either update R0s or betas but not both
#sR01Emoves.accepted<-0    # note that we either update R0s or betas but not both
#sR02Emoves.accepted<-0    # note that we either update R0s or betas but not both
#sR03Emoves.accepted<-0    # note that we either update R0s or betas but not both
#beta.moves.accepted<-0
gamma.moves.accepted<-0
lambda.moves.accepted<-0
mu.moves.accepted<-0
rho.moves.accepted<-0
proportion.of.transmission.early.in.infection.moves.accepted<-0
time.infected1.moves.accepted<-0
time.infected2.moves.accepted<-0
time.infected3.moves.accepted<-0
propSymptoms.moves.accepted<-0
disp.moves.accepted<-0
w1.moves.accepted<-0
w2.moves.accepted<-0
w3.moves.accepted<-0

# In this version parameters in the model are pi1 (==R01E/(R01E+R01D), pi2, pi3 beta1D, beta2D , beta3d, gamma, rho, lambda, mu
# and the times first infected. 
# We have priors for all of these. Since updating any parameter on it's own would change pi - for any updates we recaculate beta?E terms 

for(i in 1:N){
  #  updates to beta1D, beta2D, and beta3d (not blocked together)
  #proposed.parameters$R01<-generate.proposal.R0(parameters$R01) # proposed update
  proposed.beta1D<-generate.proposal.beta(parameters$R01D*parameters$rho) # proposed update
  proposed.parameters$R01D<-(proposed.beta1D/parameters$rho) 
  # also need to update the R01E turn to maintain the same pi1 value (just calc these by solving for R01E in terms of pi etc)
  #old line (wrong - though it would be OK if gamma changed to rho): proposed.parameters$R01E<-parameters$lambda*currentpi1*proposed.beta1D/(parameters$gamma * parameters$lambda - parameters$gamma*parameters$lambda*currentpi1)
  # but line below is simpler (and correct)
  proposed.parameters$R01E<-proposed.parameters$R01D*currentpi1/(1-currentpi1)
    
  #current.prior.density for beta1D (no term for beta1E  as the is modelled using pi term )
  current.prior.density<-beta.prior.density(parameters$R01D*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta1D)
  
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=FALSE, SEEIR=FALSE,SEIRW=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      betaD1moves.accepted<-betaD1moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  #proposed.parameters$R02<-generate.proposal.R0(parameters$R02) # proposed update
  proposed.beta2D<-generate.proposal.beta(parameters$R02D*parameters$rho) # proposed update
  proposed.parameters$R02D<-(proposed.beta2D/parameters$rho) 
  # also need to update the R01E turn to maintain the same pi2 value 
  #proposed.parameters$R02E<-parameters$lambda*currentpi2*proposed.beta2D/(parameters$gamma * parameters$lambda - parameters$gamma*parameters$lambda*currentpi2)
  proposed.parameters$R02E<-proposed.parameters$R02D*currentpi2/(1-currentpi2)
  
  
  #current.prior.density<-R0.prior.density(parameters$R01)
  current.prior.density<-beta.prior.density(parameters$R02D*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta2D)
  
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      betaD2moves.accepted<-betaD2moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  
  proposed.beta3D<-generate.proposal.beta(parameters$R03D*parameters$rho) # proposed update
  proposed.parameters$R03D<-(proposed.beta3D/parameters$rho)
  
  #proposed.parameters$R03E<-parameters$lambda*currentpi3*proposed.beta3D/(parameters$gamma * parameters$lambda - parameters$gamma*parameters$lambda*currentpi3)
  proposed.parameters$R03E<-proposed.parameters$R03D*currentpi3/(1-currentpi3)
  
  
  #current.prior.density<-R0.prior.density(parameters$R01)
  current.prior.density<-beta.prior.density(parameters$R03D*parameters$rho)
  proposal.prior.density<-beta.prior.density(proposed.beta3D)
    
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      betaD3moves.accepted<-betaD3moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  }  
  
  #   update pi1 pi2 and pi3 and find the corresponding R01E R02E R03E 
  proposed.pi1<- generate.proposal.pi(currentpi1)
  proposed.parameters$R01E<- proposed.pi1 * proposed.parameters$R01D/(1- proposed.pi1)
  current.prior.density<-pi.prior.density(currentpi1)
  proposal.prior.density<-pi.prior.density(proposed.pi1)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero 
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      currentpi1<-proposed.pi1
      current.loglikelihood<- proposal.loglikelihood
      pi1.moves.accepted<-pi1.moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
   } else {
    proposed.parameters<-parameters
  } 
  
  proposed.pi2<- generate.proposal.pi(currentpi2)
  proposed.parameters$R02E<- proposed.pi2 * proposed.parameters$R02D/(1- proposed.pi2)
  current.prior.density<-pi.prior.density(currentpi2)
  proposal.prior.density<-pi.prior.density(proposed.pi2)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero 
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      currentpi2<-proposed.pi2
      current.loglikelihood<- proposal.loglikelihood
      pi2.moves.accepted<-pi2.moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  } 

  proposed.pi3<- generate.proposal.pi(currentpi3)
  proposed.parameters$R03E<- proposed.pi3 * proposed.parameters$R03D/(1- proposed.pi3)
  current.prior.density<-pi.prior.density(currentpi3)
  proposal.prior.density<-pi.prior.density(proposed.pi3)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero 
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    #print(c("R0 accept prob",acceptance.prob ))
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      currentpi3<-proposed.pi3
      current.loglikelihood<- proposal.loglikelihood
      pi3.moves.accepted<-pi3.moves.accepted+1
    }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else {
    proposed.parameters<-parameters
  } 


# for all updates to parameters that would change pi1 pi2 or p3 (where pi1 - R01E/(R01E + R01D)) we modify corresponding betaxE
# to maintain the same pi1 pi2 and pi3, and also update R01D and R01D values. As in this parameterisation the parameters do not include betaE terms but instead these are calculateed 
# to keep pi1 etc the same. Also applies to other parameters that would effect pi1. Point of this is so that we can put informative priors on the pi terms
 
# *** updates to betaxE below commented out because now we update pi1 instead (which as a side effect update betaxE terms)  ***
#   # note formula for R01E below accounts for environmental reservoir with flow in rate lambda and flow out rate mu, with units of person equivalents in terms of infectivity
#   proposed.beta1E<-generate.proposal.beta(parameters$R01E*parameters$rho*parameters$mu/parameters$lambda) # proposed update
#   proposed.parameters$R01E<-proposed.beta1E*parameters$lambda/(parameters$rho*parameters$mu) 
#   #current.prior.density<-R0.prior.density(parameters$R01)
#   current.prior.density<-betaE.prior.density(parameters$R01E*parameters$rho*parameters$mu/parameters$lambda)
#   proposal.prior.density<-betaE.prior.density(proposed.beta1E)
#   
#   if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
#     proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=FALSE, SEEIR=FALSE,SEIRW=TRUE)
#     
#     acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
#     #print(c("R0 accept prob",acceptance.prob ))
#     if(runif(1)<acceptance.prob){  
#       parameters <- proposed.parameters  
#       current.loglikelihood<- proposal.loglikelihood
#       betaE1moves.accepted<-betaE1moves.accepted+1
#     }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
#   } else {
#     proposed.parameters<-parameters
#   }  
#   
#   proposed.beta2E<-generate.proposal.beta(parameters$R02E*parameters$rho*parameters$mu/parameters$lambda) # proposed update
#   proposed.parameters$R02E<-proposed.beta2E*parameters$lambda/(parameters$rho*parameters$mu) 
#   #current.prior.density<-R0.prior.density(parameters$R01)
#   current.prior.density<-betaE.prior.density(parameters$R02E*parameters$rho*parameters$mu/parameters$lambda)
#   proposal.prior.density<-betaE.prior.density(proposed.beta2E)
#   
#   if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
#     proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=FALSE, SEEIR=FALSE,SEIRW=TRUE)
#     
#     acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
#     #print(c("R0 accept prob",acceptance.prob ))
#     if(runif(1)<acceptance.prob){  
#       parameters <- proposed.parameters  
#       current.loglikelihood<- proposal.loglikelihood
#       betaE2moves.accepted<-betaE2moves.accepted+1
#     }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
#   } else {
#     proposed.parameters<-parameters
#   } 
#   
#   proposed.beta3E<-generate.proposal.beta(parameters$R03E*parameters$rho*parameters$mu/parameters$lambda) # proposed update
#   proposed.parameters$R03E<-proposed.beta3E*parameters$lambda/(parameters$rho*parameters$mu) 
#   #current.prior.density<-R0.prior.density(parameters$R01)
#   current.prior.density<-betaE.prior.density(parameters$R03E*parameters$rho*parameters$mu/parameters$lambda)
#   proposal.prior.density<-betaE.prior.density(proposed.beta3E)
#   
#   if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
#     proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIIR=FALSE, SEEIR=FALSE,SEIRW=TRUE)
#     
#     acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
#     #print(c("R0 accept prob",acceptance.prob ))
#     if(runif(1)<acceptance.prob){  
#       parameters <- proposed.parameters  
#       current.loglikelihood<- proposal.loglikelihood
#       betaE3moves.accepted<-betaE3moves.accepted+1
#     }  else  proposed.parameters<-parameters  # otherwise current parameters stay the same
#   } else {
#     proposed.parameters<-parameters
#   } 
  
  #  updates to gamma  
  proposed.parameters$gamma<-generate.proposal.gamma(parameters$gamma) # proposed update
  current.prior.density<-gamma.prior.density(parameters$gamma)
  proposal.prior.density<-gamma.prior.density(proposed.parameters$gamma)


  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)  
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      gamma.moves.accepted<-gamma.moves.accepted+1
    } else proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  #  updates to rho - note that these will  affect R0D and R0E  terms but not the pi terms.


  proposed.parameters$rho<-generate.proposal.rho(parameters$rho) # proposed update
 #note that changing rho changes R0D and R0E as well...so need to update that
  proposed.parameters$R01D <- parameters$R01D*parameters$rho/proposed.parameters$rho 
  proposed.parameters$R02D <- parameters$R02D*parameters$rho/proposed.parameters$rho 
  proposed.parameters$R03D <- parameters$R03D*parameters$rho/proposed.parameters$rho 
  proposed.parameters$R01E <- parameters$R01E*parameters$rho/proposed.parameters$rho 
  proposed.parameters$R02E <- parameters$R02E*parameters$rho/proposed.parameters$rho 
  proposed.parameters$R03E <- parameters$R03E*parameters$rho/proposed.parameters$rho

  current.prior.density<-rho.prior.density(parameters$rho)
  proposal.prior.density<-rho.prior.density(proposed.parameters$rho)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE) 
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      rho.moves.accepted<-rho.moves.accepted+1
    }  else proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
    #  updates to lambda - this will affect the R0E terms so to keep the same R0E (and hence pi) values we update BetaE 
#   (but since betaE is implicit we don't have to do anything - but betaE implicitly changes and R0E stays the same)
  proposed.parameters$lambda<-generate.proposal.lambda(parameters$lambda) # proposed update
 #note that changing lambda changes R0E as well...so need to update that

  current.prior.density<-lambda.prior.density(parameters$lambda)
  proposal.prior.density<-lambda.prior.density(proposed.parameters$lambda)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE) 
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood <- proposal.loglikelihood
      lambda.moves.accepted <- lambda.moves.accepted+1
    }  else proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  

    #  updates to mu -  this will affect the R0E terms so to keep the same R0E (and hence pi) values we update BetaE 
#   (but since betaE is implicit we don't have to do anything - but betaE implicitly changes and R0E stays the same)
  proposed.parameters$mu<-generate.proposal.mu(parameters$mu) # proposed update

  current.prior.density<-mu.prior.density(parameters$mu)
  proposal.prior.density<-mu.prior.density(proposed.parameters$mu)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE) 
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood <- proposal.loglikelihood
      mu.moves.accepted <- mu.moves.accepted+1
    }  else proposed.parameters<-parameters # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  
  
  #  updates to time.first.infected1, time.first.infected2, and time.first.infected3 (no longer blocked together)
  proposed.parameters$time.first.infected1<-generate.proposal.time.first.infected(parameters$time.first.infected1) # proposed update
  current.prior.density<-time.first.infected1.prior.density(parameters$time.first.infected1)
  proposal.prior.density<-time.first.infected1.prior.density(proposed.parameters$time.first.infected1)
  if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
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
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
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
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)  
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
    proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
    acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    if(runif(1)<acceptance.prob){  
      parameters <- proposed.parameters  
      current.loglikelihood<- proposal.loglikelihood
      propSymptoms.moves.accepted<-propSymptoms.moves.accepted+1
    } else proposed.parameters<-parameters  # otherwise current parameters stay the same
  } else  proposed.parameters<-parameters
  
  
 
    #  updates to  proportion.of.transmission.early.in.infection (only relevant for SEIIR model)

# #   
  # proposed.parameters$proportion.of.transmission.early.in.infection <-generate.proposal.proportion.of.transmission.early.in.infection(parameters$proportion.of.transmission.early.in.infection) # proposed update
  # current.prior.density<- proportion.of.transmission.early.in.infection.density(parameters$proportion.of.transmission.early.in.infection)
  # proposal.prior.density<-proportion.of.transmission.early.in.infection.density(proposed.parameters$proportion.of.transmission.early.in.infection)
  # if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
    # proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE)
    
    # acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
    
    # if(runif(1)<acceptance.prob){  
      # parameters <- proposed.parameters  
      # current.loglikelihood<- proposal.loglikelihood
      # proportion.of.transmission.early.in.infection.moves.accepted<-proportion.of.transmission.early.in.infection.moves.accepted+1
    # } else proposed.parameters<-parameters  # otherwise current parameters stay the same
  # } else  proposed.parameters<-parameters
  
  
  
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
      proposal.loglikelihood<-getLLHepEmod.mh.combined.no.watsan(proposed.parameters, xstart1,xstart2,xstart3, plot=FALSE, SEIRW=TRUE) 
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
    temp.pi1<-parameters$R01E/(parameters$R01D + parameters$R01E)
    temp.pi2<-parameters$R02E/(parameters$R02D + parameters$R02E)
    temp.pi3<-parameters$R03E/(parameters$R03D + parameters$R03E)
    print(c(i,c(as.numeric(parameters),temp.pi1, temp.pi2, temp.pi3, current.loglikelihood$LL)))
  }
} #  end of main MCMC loops



#  Summzarize accepteance probs

print(c("Accept ratio for pi1", pi1.moves.accepted/N))
print(c("Accept ratio for pi2", pi2.moves.accepted/N))
print(c("Accept ratio for pi3", pi3.moves.accepted/N))

print(c("Accept ratio for betaD1", betaD1moves.accepted/N))
print(c("Accept ratio for betaD2", betaD2moves.accepted/N))
print(c("Accept ratio for betaD3", betaD3moves.accepted/N))
#print(c("Accept ratio for betaE1", betaE1moves.accepted/N))
#print(c("Accept ratio for betaE2", betaE2moves.accepted/N))
#print(c("Accept ratio for betaE3", betaE3moves.accepted/N))

#print(c("Accept ratio for R01", R01moves.accepted/N))
#print(c("Accept ratio for R02", R02moves.accepted/N))
#print(c("Accept ratio for R03", R03moves.accepted/N))
print(c("Accept ratio for gamma", gamma.moves.accepted/N))
print(c("Accept ratio for rho", rho.moves.accepted/N))
print(c("Accept ratio for lambda", lambda.moves.accepted/N))
print(c("Accept ratio for mu", mu.moves.accepted/N))

print(c("Accept ratio for proportion.of.transmission.early.in.infection.moves.accepted", proportion.of.transmission.early.in.infection.moves.accepted/N))
print(c("Accept ratio for time.infected 1", time.infected1.moves.accepted/N))
print(c("Accept ratio for time.infected 2", time.infected2.moves.accepted/N))
print(c("Accept ratio for time.infected 3", time.infected3.moves.accepted/N))
print(c("Accept ratio for prop symptomatic 3", propSymptoms.moves.accepted/N))
print(c("Accept ratio for disp", disp.moves.accepted/N))
print(c("Accept ratio for w1", w1.moves.accepted/N))
print(c("Accept ratio for w2", w2.moves.accepted/N))
print(c("Accept ratio for w3", w3.moves.accepted/N))
#save.image()
save.image("model7 SEIRW MH1.12  5M thin 100 - pi1e1 priors .RData") 

# plot selected runs 
pdf("modv1 - selected fits from 2M (thin 100) chain pi1e1 prior.pdf")
par(mfrow=c(3,2),mar=c(4,4,1,1))
parameters<-posterior.samples[1,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW=TRUE)

#
parameters<-posterior.samples[2000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW =TRUE)

parameters<-posterior.samples[5000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW =TRUE)

parameters<-posterior.samples[9000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW =TRUE)

parameters<-posterior.samples[12000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW =TRUE)

parameters<-posterior.samples[15000,]
#getLLHepEmod.mh.combined.wpp.watsan(parameters, xstart1,xstart2,xstart3, Uganda.wpp, plot=TRUE)
getLLHepEmod.mh.combined.no.watsan(parameters, xstart1,xstart2,xstart3, plot=TRUE, SEIRW =TRUE)
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
pdf("modv1.11 SEIRW prior pi1e1 NBLL  2M runs thin 100.pdf")
par(mfrow=c(3,5),mar=c(4,4,1,1))
rows.to.plot<-seq(1,length(posterior.samples$R01D),1)
plot(posterior.samples$R01D[rows.to.plot],type='l',main="R01D")
plot(posterior.samples$R02D[rows.to.plot],type='l',main="R02D")
plot(posterior.samples$R03D[rows.to.plot],type='l',main="R03D")
plot(posterior.samples$R01E[rows.to.plot],type='l',main="R01E")
plot(posterior.samples$R02E[rows.to.plot],type='l',main="R02E")
plot(posterior.samples$R03E[rows.to.plot],type='l',main="R03E")

plot(posterior.samples$time.first.infected1[rows.to.plot],type='l',main="Time first infection 1")
plot(posterior.samples$time.first.infected2[rows.to.plot],type='l',main="Time first infection 2")
plot(posterior.samples$time.first.infected3[rows.to.plot],type='l',main="Time first infection 3")
plot(posterior.samples$gamma[rows.to.plot],type='l',main="gamma")
plot(posterior.samples$rho[rows.to.plot],type='l',main="rho")
plot(posterior.samples$proportion.symptomatic[rows.to.plot],type='l',main="Proportion symptomatic")
plot(posterior.samples$mu[rows.to.plot],type='l',main="mu")
plot(posterior.samples$lambda[rows.to.plot],type='l',main="lambda")


#plot(posterior.samples$w1[rows.to.plot],type='l',main="watsan camp 1")
#plot(posterior.samples$w2[rows.to.plot],type='l',main="watsan camp 2")
#plot(posterior.samples$w3[rows.to.plot],type='l',main="watsan camp 3")

dev.off()



# plot correlation
pdf("modv1.7 watsan inf prior  NBLL 3M correlation.pdf")
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

dbar<- -2*mean(posterior.samples$LL)
post.means<-data.frame(t(mean(posterior.samples)))
#dhat<-  -2*getLLHepEmod.mh.combined.wpp.watsan(post.means, xstart1,xstart2,xstart3,Uganda.wpp, plot=FALSE)$LL
dhat<-  -2*getLLHepEmod.mh.combined.no.watsan(post.means, xstart1,xstart2,xstart3, plot=FALSE, SEEIR=TRUE)$LL

pD<-dbar-dhat
DIC<- dbar+pD
dbar
DIC
pD

# 1. latent period

# mean

mean(365/posterior.samples$gamma)
quantile( 365/posterior.samples$gamma, c(.025,.5,.975))

# 2. infectious period

# mean

mean(365/posterior.samples$rho)
quantile( 365/posterior.samples$rho, c(.025,.5,.975))

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
  # binomial likelihood
  #   q.025prob.observed1<-rbind(q.025prob.observed1,qbinom(.025,round(wk.inc1),parameters$proportion.symptomatic))
  #   q.975prob.observed1<-rbind(q.975prob.observed1,qbinom(.975,round(wk.inc1),parameters$proportion.symptomatic))
  #   q.025prob.observed2<-rbind(q.025prob.observed2,qbinom(.025,round(wk.inc2),parameters$proportion.symptomatic))
  #   q.975prob.observed2<-rbind(q.975prob.observed2,qbinom(.975,round(wk.inc2),parameters$proportion.symptomatic))
  #   q.025prob.observed3<-rbind(q.025prob.observed3,qbinom(.025,round(wk.inc3),parameters$proportion.symptomatic))
  #   q.975prob.observed3<-rbind(q.975prob.observed3,qbinom(.975,round(wk.inc3),parameters$proportion.symptomatic))
  
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

pdf("Predicted versus observed NB 1.25m thin 100.pdf",height=6,width=5,pointsize=10)
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
#out1<-ode(y=x1,times=times1, func="derivs",parms=parameters1,jacfun="jac",dllname="SEIRforced", initforc = "forcc", forcings=forcings1, initfunc="initmod",nout=1,outnames="Sum")

