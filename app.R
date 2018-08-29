# HEV Shiny app for evaluating HEV vaccination options in response to epidemics in emergency situations
# Ben Cooper, August 2018 
# This app is intended to accompany this paper https://www.biorxiv.org/content/early/2017/11/16/219154 ()

library("shiny")
library("deSolve")
library("ggplot2")
library("scales")


HepEwithVAxmod.shiny <- function(t, x, parms) { 
  # This is the standard SEIR model. Vaccination is not in ODEs but specified by events in the call to ode 
  #v3  allows for a different efficacy after second dose of vaccine - but in this case 
  # vaccination is implemented as an event in the above function immunisation.event (avoids numerical problems when vaccination is a rate)
  # collects output data on the number of infections in the vaccinated and unvaccinated populations  
  # Also to allow for the fact that a proportion of the population are not elligible for vaccination we divde susceptibles into 
  # 2 groups SV (susceptible vaccine candidates) and SNV susceptible non vaccine candidates
  # we also assume that the time it takes to generate immunity is fixed (rather than exponentially distributed) so 
  # delay from the vaccination time to immunity is dealt with when specify the time of the immunisation (should account for delay from vaccination)
  
  
  # t is time, 
  # x is the vector of state variables (SNV,SV VS, VE,  EV, ENV, IV, INV, RV, RNV) # where SV is susceptible and unvaccinated  but vaccinatable, SNV is  Susceptible unvaccinated and unvaccinatable, VS is unusccesfully vaccinated (and still scuscptible) 
  # and VE is successfully vaccinated and waiting for immunity, EV and ENV are respectively latently infected and vaccinatable and latently infectied and non-vaccinatable. Similarly for RV and RNV (both infec)
  # we also track i) the number of people vaccinated with the first dose of vaccine V1
  # the number of poeple vaccinated with the second dose, V2, and cumulative number of HepE in fections, CumInc
  
  
  #  Note that if we have an imperfect vaccine  a real vaccination program wouldn't vaccinate people more than once in a each vaccination episode 
  # (where an episode is a  short period lasting a few days when a camp is vaccinated) - so to allow for this fact if vaccinated unsuccesfully 
  # vaccinated people move to a compartment VS (vaccinated and still susceptible). If susccessfully vaccinated they move to VE where they may still be infected. 
  # but after a while progress to R from VE as immunity developes. 
  # We assume that patients in the VE compartment can still be infected, and that vaccine has no effect on anyone who is in the E, I or R state when vaccinated. 
  # We want to model the use of up to two doses of vaccine (typically about a month apart). When doing this we assume that for both the first and second vaccination
  # episodes everyone in the target age group is equally likely to be vaccinated (and therefore that someone missing the first vaccine is no more or less likely to miss the
  # next vaccine). Those who are vaccinated effectiveley but not yet immune (VE) can still be infected.
  # Finally, in this version of the model we assume that each dose of vaccine to someone susceptible (S or VS) either gives full protection (with some fixed probabaility that is 
  # the same for the first or second dose) or no protection.
  # For vaccination, we assume that vaccine only has any effect if used in those who are in S state or VS (vaccinated but still suceptible due to vaccine failure)
  
  with(as.list(c(parms, x)), {

    # N is current population size
    # gamma is rate of prgression to exposed to infectious
    # rho is rate of progression from infectious to recovered.
    # vacc.eff1 is vaccine efficacy of first dose
    # vacc.eff2 is vaccine efficacy of second dose
    #  proportion.vaccinatable is the proportion of the population in the vaccine candidate population (e.g. proportion over 15 and not pregnant)
    #  proportion.to.vaccinate is the proportion of the vaccine candidate population to be vaccinated at each vaccination episode
   
    
    N<-SNV + SV +ENV + EV+INV +IV +RNV +RV +VS    

    beta<-beta.min

    # vacc.eff2given1 is probability of immunisation from vaccine 2 given lack of immunisaiton from 1. 
    # derived from the probablity of successfull immunisation after 1 or 2 doses vacc.eff1 and  vacc.eff2
    # using (1-vacc.eff2)=(1-vacc.eff1)(1-vacc.eff2given1)
    
    # Derivatives 
    
    dSNV <- -beta*SNV*(INV +IV)/N    
    dSV <- -beta*SV*(INV +IV)/N    
    dENV <-  beta*SNV*(INV +IV)/N - gamma*ENV  
    dEV <-  beta*(VS +SV)*(INV +IV)/N - gamma*EV   #vaccinatable group (no necessarily vaccinated)
    dINV <- gamma*(ENV) - rho*INV  
    dIV <- gamma*(EV) - rho*IV  
    dRNV <- rho*INV 
    dRV <- rho*IV 
    dVS<- - beta*VS*(INV +IV)/N  
    dV1<-  0
    dV2<-  0
    dCumIncNV<-gamma*ENV  #cumulative incidence of infection in not vaccinatable group
    dCumIncV<-gamma*EV    #cumulative incidence of infection in vaccinatable group
    
    # Return values 
    list(c(dSNV,dSV,dENV,dEV,dINV, dIV,dRNV, dRV, dVS, dV1, dV2,dCumIncNV, dCumIncV ))
    
  })  # end with
} # end HepEwithVAxmod.shiny


immunisation.event<-function(t, y,parms){
  #  update state variable for HepEwithVAxmod model
  #  Here vacc.eff1 is prob of immunity from vaccine after first dose. 
  #  and vacc.eff2 is prob of immunity after having received two dose
  #  Assume that vaccination has no effected on someone on the E, I or R classes so don't need to split these into vaccinated and not
  with(as.list(c(y,parms)),{
    vacc.efficacy.given.prior.vacc<-max(0,(vacc.eff2-vacc.eff1)/(1-vacc.eff1) )# i.e prob of getting immunity given susceptible and previously vaccinated
    newVS <- VS*(1- proportion.to.vaccinate* vacc.efficacy.given.prior.vacc) + SV * (1- vacc.eff1)* proportion.to.vaccinate         # susceptibles in vaccinatable population , who have been vaccinated but are still susceptible
    newRV<-RV+ vacc.eff1*SV*proportion.to.vaccinate + vacc.efficacy.given.prior.vacc*proportion.to.vaccinate*VS
    newSV<-SV * (1- proportion.to.vaccinate) # susceptibles in vaccinatable population , not yet vaccinated 
    oldV2<-V2
    V2<- V2 + (V1-V2)*proportion.to.vaccinate # number having recieved at least two additional doses (following initial dose) 
    V1<- V1+ (SV +EV +IV +RV+ VS - V1)* proportion.to.vaccinate # number having receIved at least one additional dose (following initial dose) = no. vaccinatable who have not been vaccinated with 2nd dose already * proportion of those to vaccinate + number who had previously received at least one dose. 
    VS<-newVS
    SV<-newSV
    RV<-newRV
    return(c(SNV,SV, ENV,EV,INV, IV,RNV, RV, VS, V1, V2,CumIncNV, CumIncV ))
  })
}


shinyUI<-navbarPage(
  title = 'HEV vaccination policy evaluation tool',
  tabPanel('Model',  
           pageWithSidebar(
             
             #  Application title
             headerPanel(""),
             
             # Sidebar with sliders that demonstrate various available options
             sidebarPanel(
               
               tags$h4("Model parameters"),
               
               sliderInput("N", "Population size:", 
                           min=10^3, max=10^5, value=17000, step=10^3),
               
               sliderInput("E0", "Number latently infected at time 0:", 
                           min=0, max=10^2, value=1, step=1),
               sliderInput("R.at.time.zero", "Percentage of population immune at time 0", 
                           min=0, max=100, value=0, step=1), 
               sliderInput("R0", "Basic reproduction number", 
                           min=0, max=20, value=6.5, step=.1),
               
               sliderInput("One.over.rho", "Mean infectious period (days):", 
                           min=10, max=100, value=36, step=1),
               
               sliderInput("One.over.gamma", "Mean latent period (days):", 
                           min=10, max=100, value=34, step=1), 
               
               sliderInput("percent.symptomatic", "Percentage of infections which are symptomatic & reported", 
                           min=1, max=100, value=12.5, step=.5), 
               
               sliderInput("percent.sympt", "Percentage of infections which are symptomatic", 
                           min=0, max=100, value=20, step=1), 
               
               sliderInput("prob.death.pregnant", "Probability of death given symptomatic infection if pregnant", 
                           min=0, max=1, value=.2, step=0.01) ,
               sliderInput("prob.death.not.pregnant", "Probability of death given symptomatic infection if not pregnant", 
                           min=0, max=1, value=.02, step=0.01) ,
               
               sliderInput("Percent.age.15to65", "Percentage of population aged 15-65", 
                           min=0, max=100, value=49, step=1), 
               sliderInput("Percent.pregnant", "Percentage of population who are pregnant", 
                           min=0, max=10, value=3, step=.5) 
               
             ),
             
             # Show a table summarizing the values entered
             mainPanel(
               #  h5("What to plot:"),
               checkboxGroupInput("Indicators", "",
                                  c("Susceptible", 
                                    "Latently infected", 
                                    "Infectious",
                                    "Immune",
                                    "Observed cases",
                                    "Cumulative infections",
                                    "Vaccinated with 2 doses",
                                    "Vaccinated with 3 doses"
                                  ),
                                  selected=c(
                                    "Susceptible", 
                                    "Latently infected", 
                                    "Infectious",
                                    "Immune",
                                    "Observed cases",
                                    "Cumulative infections",
                                    "Vaccinated with 2 doses",
                                    "Vaccinated with 3 doses"
                                  ),
                                  inline=TRUE),
               plotOutput("graph1"),
               h4("Who to vaccinate"),
               checkboxGroupInput("WhoToVax", "",
                                  c("Ages 16-65 yrs, excluding pregnant women"="Age16to65",
                                    "Pregnant women"="Preg",
                                    "Everyone else"="EveryoneElse"
                                  ),
                                  selected=c(
                                    "Age16to65"
                                  ),
                                  inline=TRUE
               ),
               
               fluidRow(
                 
                 column(4,
                        sliderInput("vax.coverage", "Vaccine coverage in target groups (%)", 
                                    min=0, max=100, value=90, step=1)
                 ),
                 
                 column(4,
                        sliderInput("cases.before.vaccine", "Number of reported cases before vaccine is used", 
                                    min=0, max=100, value=50, step=1) 
                 ),
                 column(4,
                        sliderInput("t.to.immunity", "Days from effective vaccination to immunity", 
                                    min=0, max=30, value=14, step=1)
                 ), 
                 column(4,
                        sliderInput("vacc.eff1", "Vaccine effectiveness after two doses (%) ", 
                                    min=-0, max=100, value=80, step=1)
                 ), 
                 column(4,
                        sliderInput("vacc.eff2", "Vaccine effectiveness after three doses (%) ", 
                                    min=-0, max=100, value=93, step=1)
                 ), 
                 column(4,
                        sliderInput("t.dose2.to.3", "Days from second to third dose", 
                                    min=0, max=200, value=182, step=1)  )
               ), 
               h6("Note that vaccine has only been evaluated in those aged 16-65 who are not pregnant. Equal efficacy in other groups is assumed."),
               
               plotOutput("graph2"),
               
               h4("Outcomes without vaccination"),
               
               tableOutput("datatable.novax"),
               
               h4("Outcomes with vaccination"),
               
               tableOutput("datatable.vax")
               
             )
           )
           
           ),

  tabPanel('Information',        
           mainPanel( 
             h3("Hepatitis E background information"),
             p("Hepatitis E Virus (HEV) is the leading cause of acute viral hepatitis globally. Symptomatic infection is associated with case fatality rates of ~20% in pregnant women 
                and ~2% in others [1]; it is estimated to account for ~10,000 annual pregnancy-related deaths in southern Asia alone [2]."),
             
             p("Recently, large and well-documented outbreaks with high mortality have occurred in internally-displaced persons (IDP) camps in Sudan, Uganda and South Sudan. 
               Controlling these outbreaks is difficult, but there  is a safe and effective HEV vaccine that is licensed in China for those aged 16-65 years who are not pregnant [3]. 
              In 2015 the Strategic Advisory Group of Experts on immunization could not recommend the routine use of the vaccine for population sub-groups including children 
                aged less than 16 years and pregnant women but emphasized that the use of the vaccine during outbreaks of hepatitis E should be considered [4]."), 

              p(" The purpose of this app is  to provide a tool to help evaluate potential  strategies for using the vaccine  in such outbreak settings."),
            
             h3("What this app does"),
             
             p("This app allows for  analysis using Model 1a from Cooper et al., 2018 [5]. This is the basic SEIR model which was shown to provide the best fit to outbreak data 
               from Uganda. An SEIR model divides up the population into those who are susceptible to infection (S), exposed but not yet infectious (E), infectious (I), and immune and therefore reisistant to infection (R).
               The app allows users to explore the effects of different assumptions about parameter values, and different assumptions about the timing, coverage, and effectiveness of the vaccine."),
        
             h3("Assumptions"),
             p("The app uses baseline parameter values for  HEV transmission estimated in the modelling paper [5] by fitting the model to outbreak data from Uganda, though users can freely adjust all parameter values."),
             p("Baseline assumptions about vaccine effectiveness after 2 and 3 doses are based on analysis of data from Zhu et al. [3], as described in [5].
                Since these data do not allow evaluation of the effectiveness of a single dose, we make the conservative assumption that a single dose has no effect.
               Note also that the vaccine is currently only licensed for those aged 16-65 years who are not pregnant. There are, however, reasons for thinking it might also be effective 
                in pregnant women. The app therefore provides options for  vaccinating pregant women and other age groups under the assumption that the vaccine effectiveness is the same across all groups."),
             p("The recommended vaccine dose schedule is for three doses, speparated by intervals of one and five  months; these are the default intervals in the app. The app allows users to vary the time between the second and third doses,  though the impact of doing this on effectiveness is not known."),
             p("Full details of further assumptions, consideration of additional models, results of model fitting and a fuller discussion of caveats and limitations can be found in the paper describing analysis of HEV outbreak data using this and other models [5]."),

             h3("References"),
             a(href="https://aasldpubs.onlinelibrary.wiley.com/doi/full/10.1002/hep.25505", "[1] Rein D, Stevens G, Theaker J, Wittenborn J, Wiersma S. The global burden of hepatitis E virus genotypes 1 and 2 in 2005. Hepatology. 2012;55: 988–997."),
               br(),     
             a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3437697/", "[2] Labrique AB, Sikder SS, Krain LJ, West KP Jr, Christian P, Rashid M, et al. Hepatitis E, a vaccine-preventable cause of maternal deaths. 
              Emerg Infect Diseases. 2012;18: 1401–1404"), 
             br(),  
             a(href="http://www.natap.org/2010/newsUpdates/hepvaccine.pdf", "[3] Zhu FC, Zhang J, Zhang XF, Zhou C, Wang ZZ, Huang SJ, et al. Efficacy and safety of a recombinant hepatitis E vaccine in healthy adults: 
               a large-scale, randomised, double-blind placebo-controlled, phase 3 trial. Lancet. 2010; 376: 895–902."),
             br(),  
             a(href="http://www.who.int/wer/2015/wer9018.pdf?ua=1", "[4] World Health Organization. Hepatitis E vaccine: WHO position paper, May 2015.Weekly Epidemiological Record 2015; 90,185–200."),
             br(),  
             a(href="https://www.biorxiv.org/content/biorxiv/early/2017/11/16/219154.full.pdf", "[5] Cooper BS, White LJ, Siddiqui MR. Reactive and pre-emptive vaccination strategies to control hepatitis E infection in emergency and refugee settings: a modelling study. BioRxiv doi: https://doi.org/10.1101/219154"),
             br(),  
             a(href="https://aasldpubs.onlinelibrary.wiley.com/doi/pdf/10.1002/hep.25522", "[6] Wu T, Zhu FC, Huang SJ, Zhang XF, Wang ZZ, Zhang J, et al. Safety of the hepatitis E vaccine for pregnant women: a preliminary analysis. Hepatology. 2012;55: 2038."),
             h3("Contact information"),
             p("This app was created by Ben Cooper, Nuffield Dept of Medicine, University of Oxford. It is based on modelling work that was done in collaboration with Lisa White, University of Oxford 
                and Ruby Siddiqui, Médecins Sans Frontières."),
             p("For further information about either this app or the underlying modelling work please contact:"), 
             tags$a("Ben Cooper", href="https://www.ndm.ox.ac.uk/principal-investigators/researcher/ben-cooper")
            )
           )
 
)

              


#shinyServer(function(input, output) {
server <- function(input, output) {
  dataInput <- reactive({
    # Model Parameters:
    N  <- input$N   # Total Population
    R0  <- input$R0   # R0
    rho<-365/input$One.over.rho
    gamma<-365/input$One.over.gamma
    proportion.symptomatic<-input$percent.symptomatic/100
    times<-seq(0,3,1/52)
    cases.b4.vacc<-input$cases.before.vaccine 
    time.dose.3<-input$t.dose2.to.3/365
    #immunisation.times<-c(1/12,6/12)  # time of second and third doses relative to first (in years) 
    immunisation.times<-c(1/12,time.dose.3) 
    
    immun.delay1<-input$t.to.immunity/365  # delay from vaccination to immunity in years
    proportion.pregnant<-input$Percent.pregnant/100
    proportion.15to65<-input$Percent.age.15to65/100
    proportion.vaccinatable<-0
    if(c("Age16to65") %in% input$WhoToVax ){
      proportion.vaccinatable<-proportion.vaccinatable+ proportion.15to65 - proportion.pregnant
    }
    pregnant.vaccinated<-FALSE
    if(c("Preg") %in% input$WhoToVax ){
      proportion.vaccinatable<-proportion.vaccinatable+  proportion.pregnant
      pregnant.vaccinated<-TRUE
    }
    
    if(c("EveryoneElse") %in% input$WhoToVax ){
      #  proportion.vaccinatable<-proportion.vaccinatable+ 1-(proportion.15to65 + proportion.pregnant)
      proportion.vaccinatable<-proportion.vaccinatable+ 1-proportion.15to65 
    }
    
    proportion.of.these.to.vax<-input$vax.coverage/100
    
    VE1<-input$vacc.eff1 
    VE2<-input$vacc.eff2
    E0<-input$E0    # numbers initially latently infected 
    R.at.time.zero<-input$R.at.time.zero*N/100 # numbers initially immune
    
    parameters<-c(
      beta.min= R0*rho , #  currently no seasonality
      beta.max= R0*rho,
      phi=0,  # phase angle of seasonal forcing in radians (specifies when peak occurs). Currently no seasonality.
      gamma=gamma,
      rho=rho,
      vacc.eff1=VE1/100,   # vaccine efficacy:  proportion of people successfully immunized by two dose
      vacc.eff2=VE2/100,   # vaccine efficacy: proportion of people successfully immunized by three doses
      #proportion.vaccinatable= proportion.vaccinatable,
      #proportion.to.vaccinate= proportion.of.these.to.vax
      proportion.vaccinatable=proportion.vaccinatable,
      proportion.to.vaccinate= proportion.vaccinatable * proportion.of.these.to.vax
    )
    
    S0<-N - E0 - R.at.time.zero
    I0<-0
    R00<-R.at.time.zero #don't confuse wiht R0 the basic reproduction number
    VS0<-0
    V10<-0
    V20<-0
    xstart<-c(
      SNV=S0*(1-proportion.vaccinatable),  # susceptibles not eligible to be vaccinated given current policy
      SV=S0*proportion.vaccinatable,       # susceptibles  eligible to be vaccinated
      ENV=E0,
      EV=0,
      INV=I0 * (1- proportion.vaccinatable ),
      IV=I0 * proportion.vaccinatable ,
      RNV=R00 * (1- proportion.vaccinatable ),
      RV=R00 * proportion.vaccinatable ,
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
    no.vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.shiny, params.novax)
    
    #  Now run the model with vaccination
    if(cases.b4.vacc>0) {
      #   then   determine when first day when there  are  cases.b4.vacc - and delay immunisation times by this much
      no.vax.cum.inc<-(no.vax.out[,13]+no.vax.out[,14])*proportion.symptomatic
      n.time.points<-length(no.vax.cum.inc)
      first.time.gt.Xcases<-no.vax.out[min(n.time.points,1+findInterval(c(-Inf, cases.b4.vacc), no.vax.cum.inc)[2]),1]
      
      new.immunisation.times <-times[findInterval(immunisation.times+ first.time.gt.Xcases +immun.delay1,times)]
      vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.shiny,parameters,events=list(func=immunisation.event,time= new.immunisation.times))
      
    } else {
      #print(parameters)
      immunisation.times<-c(0,0.02) # to model pre-emptive vaccination assume 2nd dose is given an time 0, and third v shortly after
      immunisation.times <-times[findInterval(immunisation.times + immun.delay1,times)] #this ensures immunisation time is at included time point
      vax.out<-ode(y=xstart,times=times, func=HepEwithVAxmod.shiny,parameters,events=list(func=immunisation.event,time= immunisation.times))
    }
    
    list(no.vax.out=no.vax.out, vax.out=vax.out)
    
    
    # Now make long tables for plotting and wide tables for outputting
    S.novax<-no.vax.out[,2]+no.vax.out[,3]
    E.novax<-no.vax.out[,4]+no.vax.out[,5]
    I.novax<-no.vax.out[,6]+no.vax.out[,7]
    R.novax<-no.vax.out[,8]+no.vax.out[,9]
    O.novax<-I.novax*input$percent.symptomatic/100
    C.novax<-no.vax.out[,13]+no.vax.out[,14]
    V1.novax<-no.vax.out[,11]  
    V2.novax<-no.vax.out[,12]  
    C.novax.NV<-no.vax.out[,13] #cumulative infections in non-vaccinatable population
    C.novax.V<-no.vax.out[,14] #cumulative infections in vaccinatable population
    
    
    S.vax<-vax.out[,2]+vax.out[,3]
    E.vax<-vax.out[,4]+vax.out[,5]
    I.vax<-vax.out[,6]+vax.out[,7]
    R.vax<-vax.out[,8]+vax.out[,9]
    O.vax<-I.vax*input$percent.symptomatic/100
    C.vax<-vax.out[,13]+vax.out[,14]
    V1.vax<-vax.out[,11]  
    V2.vax<-vax.out[,12]  
    C.vax.NV<-vax.out[,13] #cumulative infections in non-vaccinatable population
    C.vax.V<-vax.out[,14] #cumulative infections in vaccinatable population
    
    np<-length(times)
    
    long.novax <- data.frame(
      Year=rep(times,8), 
      Numbers.novax = c(S.novax, E.novax, I.novax, R.novax, O.novax, C.novax,V1.novax,V2.novax), 
      Indicator=rep(c("Susceptible", 
                      "Latently infected", 
                      "Infectious",
                      "Immune",
                      "Observed cases",
                      "Cumulative infections",
                      "Vaccinated with 2 doses",
                      "Vaccinated with 3 doses" ), 
                    each=np))
    wide.novax <- cbind(times, S.novax, E.novax, I.novax, R.novax, O.novax, C.novax,V1.novax,V2.novax,C.novax.NV,C.novax.V)
    
    
    long.vax <- data.frame(
      Year=rep(times,8), 
      Numbers.vax = c(S.vax, E.vax, I.vax, R.vax, O.vax, C.vax,V1.vax,V2.vax), 
      Indicator=rep(c("Susceptible", 
                      "Latently infected", 
                      "Infectious",
                      "Immune",
                      "Observed cases",
                      "Cumulative infections",
                      "Vaccinated with 2 doses",
                      "Vaccinated with 3 doses" ), 
                    each=np))
    wide.vax <- cbind(times, S.vax, E.vax, I.vax, R.vax, O.vax, C.vax,V1.vax,V2.vax,C.vax.NV,C.vax.V)
    
    list(long.novax=long.novax, wide.novax=wide.novax,long.vax=long.vax, wide.vax=wide.vax)
    
  })
  
  output$graph1 <- renderPlot({
    long <- dataInput()[["long.novax"]]
    p <- ggplot(long[long$Indicator %in% input$Indicators,], 
                aes(x=Year, y=Numbers.novax, group=Indicator))    
    p <- p + 
      geom_line(aes(colour = Indicator), size=1, alpha=.75) + 
      ggtitle("No vaccination")+
      theme(text=element_text(size=15)) +
      scale_x_continuous(name="Year")+ 
      scale_y_continuous(labels = comma, name="")
    print(p)
  })
  
  
  output$graph2 <- renderPlot({
    long <- dataInput()[["long.vax"]]
    p <- ggplot(long[long$Indicator %in% input$Indicators,], 
                aes(x=Year, y=Numbers.vax, group=Indicator))    
    p <- p + 
      geom_line(aes(colour = Indicator), size=1, alpha=1.0) + 
      ggtitle("Vaccination")+
   theme(text=element_text(size=15)) +   
      scale_x_continuous(name="Year")+ 
      scale_y_continuous(labels = comma, name="")
    print(p)
  })
  

  output$datatable.novax <- 
    renderTable({
      table.novax <- dataInput()[["wide.novax"]]   
      total.infections<-max(table.novax[,c(7)])
      # total.infections.V<-max(table.novax[,c(10)]) #in vaccinatable pouplation
      # total.infections.NV<-max(table.novax[,c(11)]) #in non-vaccinatable pouplation
      
      proportion.symptomatic<-input$percent.symptomatic/100
      proportion.pregnant<-input$Percent.pregnant/100
      proportion.15to65<-input$Percent.age.15to65/100
      S0<-input$N - input$E0
      relative.risk.of.symptoms.given.infected.if.pregnant<-1 #in absence of any information to the contrary 
      prob.death.pregnant<-input$prob.death.pregnant
      prob.death.not.pregnant<-input$prob.death.not.pregnant
      Number.infected.pregnant<-total.infections*proportion.pregnant
      prob.pregnant.infected<-total.infections/S0
      Number.infected.not.pregnant<-total.infections-Number.infected.pregnant
      prob.symptoms.if.not.pregnant<-total.infections*proportion.symptomatic/(Number.infected.not.pregnant+ Number.infected.pregnant* relative.risk.of.symptoms.given.infected.if.pregnant)
      #checked. where relative.risk.of.symptoms.given.infected.if.pregnant = prob of symtpoms given infected and pregnant/prob of symtpoms given infected and not pregnant
      
      prob.symptoms.if.pregnant<- relative.risk.of.symptoms.given.infected.if.pregnant*prob.symptoms.if.not.pregnant
      
      Number.cases.pregnant = prob.symptoms.if.pregnant * Number.infected.pregnant
      Number.cases.not.pregnant<-prob.symptoms.if.not.pregnant * Number.infected.not.pregnant
      Number.deaths.pregnant<-prob.death.pregnant*Number.cases.pregnant 
      Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.cases.not.pregnant
      Number.cases<-Number.cases.pregnant + Number.cases.not.pregnant
      Number.deaths<-Number.deaths.pregnant+ Number.deaths.not.pregnant           
      table.data<-data.frame(Number.cases, Number.deaths.not.pregnant,  Number.cases.pregnant,Number.deaths.pregnant )
      colnames(table.data)<-c("Total symptomatic cases", "Total deaths", "Pregnant cases", "Pregnant deaths")
      table.data
    },digits=0)
  
  
  output$datatable.vax <- 
    renderTable({
      table.vax <- dataInput()[["wide.vax"]]   
      total.infections<-max(table.vax[,c(7)])
      total.infections.V<-max(table.vax[,c(11)]) #in vaccinatable pouplation
      total.infections.NV<-max(table.vax[,c(10)]) #in non-vaccinatable pouplation
      proportion.symptomatic<-input$percent.symptomatic/100
      proportion.pregnant<-input$Percent.pregnant/100
      proportion.vaccinatable<-0
      proportion.15to65<-input$Percent.age.15to65/100
      if(c("Age16to65") %in% input$WhoToVax ){  # group is 16-65 but not pregnant
        proportion.vaccinatable<-proportion.vaccinatable+ proportion.15to65 - proportion.pregnant
      }
      pregnant.vaccinated<-FALSE
      if(c("Preg") %in% input$WhoToVax ){
        proportion.vaccinatable<-proportion.vaccinatable+  proportion.pregnant
        pregnant.vaccinated<-TRUE
      }
      
      if(c("EveryoneElse") %in% input$WhoToVax ){
        proportion.vaccinatable<-proportion.vaccinatable+ 1-(proportion.15to65 )
        
      }
      
      S0<-input$N - input$E0
      relative.risk.of.symptoms.given.infected.if.pregnant<-1
      prob.death.pregnant<-input$prob.death.pregnant
      prob.death.not.pregnant<-input$prob.death.not.pregnant
      
      if(pregnant.vaccinated) {
        Number.infected.pregnant<-total.infections.V*proportion.pregnant/ proportion.vaccinatable
        prob.pregnant.infected<-total.infections.V/(S0*proportion.vaccinatable)
        proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable-proportion.pregnant
        proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)
        
      } else {
        Number.infected.pregnant<-total.infections.NV*proportion.pregnant/ (1-proportion.vaccinatable)
        prob.pregnant.infected<-total.infections.NV/(S0*(1-proportion.vaccinatable))
        proportion.vaccinatable.and.not.pregnant<-proportion.vaccinatable
        proportion.not.vaccinatable.and.not.pregnant<-(1-proportion.vaccinatable)-proportion.pregnant
      }
      
      Number.infected.not.pregnant<-total.infections-Number.infected.pregnant
      prob.symptoms.if.not.pregnant<-total.infections*proportion.symptomatic/(Number.infected.not.pregnant+ Number.infected.pregnant* relative.risk.of.symptoms.given.infected.if.pregnant)
      prob.symptoms.if.pregnant<- relative.risk.of.symptoms.given.infected.if.pregnant*prob.symptoms.if.not.pregnant
      Number.cases.pregnant = prob.symptoms.if.pregnant * Number.infected.pregnant
      Number.cases.not.pregnant<-prob.symptoms.if.not.pregnant * Number.infected.not.pregnant
      Number.deaths.pregnant<-prob.death.pregnant*Number.cases.pregnant 
      Number.deaths.not.pregnant<-prob.death.not.pregnant*Number.cases.not.pregnant
      Number.cases<-Number.cases.pregnant + Number.cases.not.pregnant
      Number.deaths<-Number.deaths.pregnant+ Number.deaths.not.pregnant           
      table.data<-data.frame(Number.cases, Number.deaths.not.pregnant,  Number.cases.pregnant,Number.deaths.pregnant)
      colnames(table.data)<-c("Total symptomatic cases", "Total deaths", "Pregnant cases", "Pregnant deaths")
      table.data
    },digits=0)
}


shinyApp(shinyUI, server)