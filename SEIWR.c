
//  seir(w) model
//
//  Created by BSC on 19/12/14.
//  Copyright 2014. All rights reserved.
//
// Implements deterministic SEIR(W) model to be used with R
// Here the W (for water) corresponds to an environmental reservoir.


// People may be infected directly from W (and persistence in W can be quite long)
// and Infectious people also shed into W
// Infections from environment follow a Holling type II functional response (linear increase with low level contamination but plateauing)
// Because units of contamination in w are not clear and not all parameters would be identifiable 
// but we can rescale W units so that rate of infection from environment is proportional to w/(1+w)
// This means we require three parameters to model env reservoir - rate of shedding to env (lambda), rate of loss from env Nu), and rate of transmission from env (beta2)
// units of W are person equivalents. So W measures the person equivalents of environmental contamination and we have same rate of 
// of transmission from person equivalents to suseptibles as from infected people to susceptibles. This should make it identifiable. 
// To compile from R: system("R CMD SHLIB SEIR.c")
// To load from R (in Unix): dyn.load("SEIWR.so")
// To call from R use deSolve package and see details in help on using compiled code.


#include </Library/Frameworks/R.framework/Versions/2.12/Resources/include/R.h> 

static double parms[6];

#define beta parms[0]
#define gamma parms[1]
#define rho parms[2]
#define lambda parms[3]
#define nu parms[4]
#define beta2 parms[5]

//  note that mu is used later so use the name nu instead

// beta is transmission rate from I
// gamma is rate of progression from E to I
// rho is rate of recovery from I
// lambda is rate of shedding from I to W
// nu is rate of decay from W (i.e. rate contamination is lost)
// beta2 is rate of transmission from the environment

/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=6;
    odeparms(&N, parms);
}

/* derivatives and one output variable 
y[0] is S (susceptibles)
y[1] is E (latently infected)
y[2] is I (infectious)
y[3] is R (recovered and immune)
y[4] is CumInc (cumulative incidence)
y[5] is W (environmental contaminatin (in water)) 
*/
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1");
    double N;
    N= y[0]+y[1]+y[2]+y[3];
    
    ydot[0] = -beta*y[0]*y[2]/N -beta2*y[0]*y[5]/(N*(1+y[5]));
    ydot[1] =  beta*y[0]*y[2]/N + beta2*y[0]*y[5]/(N*(1+y[5])) - gamma*y[1];  
    ydot[2] = gamma*y[1] - rho*y[2];  
    ydot[3] =  rho*y[2]; 
    ydot[4] = gamma*y[1];
    ydot[5]=lambda*y[2] - nu*y[5];
    yout[0]=y[0]+y[1]+y[2]+y[3];
}

/* Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    double N;
    
    N= y[0]+y[1]+y[2]+y[3];
    pd[0] =-beta*y[2]/N -beta2*y[5]/(N*(1+y[5]));
    pd[1] = beta*y[2]/N + beta2*y[5]/(N*(1+y[5]));
    pd[2] =0;
    pd[3] =0;
    pd[4]=0;
    pd[5]=0;
    pd[(*nrowpd)]=0;
    pd[(*nrowpd)+1]=-gamma;
    pd[(*nrowpd)+2]=gamma;
    pd[(*nrowpd)+3]=0;
    pd[(*nrowpd)+4]=gamma;
    pd[(*nrowpd)+5]=0;
    pd[2*(*nrowpd)]=-beta*y[0]/N;
    pd[2*(*nrowpd)+1]=beta*y[0]/N;
    pd[2*(*nrowpd)+2]=-rho;
    pd[2*(*nrowpd)+3]=rho;
    pd[2*(*nrowpd)+4]=0;
    pd[2*(*nrowpd)+5]=lambda;
    pd[3*(*nrowpd)]=0;
    pd[3*(*nrowpd)+1]=0;  
    pd[3*(*nrowpd)+2]=0;  
    pd[3*(*nrowpd)+3]=0;  
    pd[3*(*nrowpd)+4]=0; 
    pd[3*(*nrowpd)+5]=0; 
    pd[4*(*nrowpd)]=0;
    pd[4*(*nrowpd)+1]=0;  
    pd[4*(*nrowpd)+2]=0;  
    pd[4*(*nrowpd)+3]=0;  
    pd[4*(*nrowpd)+4]=0; 
    pd[4*(*nrowpd)+5]=0; 
    pd[5*(*nrowpd)]=beta2*y[0]*(y[5]/(N*(1+y[5])*(1+y[5])) - 1/(N*(1+y[5])));
    pd[5*(*nrowpd)+1]=-beta2*y[0]*(y[5]/(N*(1+y[5])*(1+y[5])) - 1/(N*(1+y[5])));   
    pd[5*(*nrowpd)+2]=0;  
    pd[5*(*nrowpd)+3]=0;  
    pd[5*(*nrowpd)+4]=0; 
    pd[5*(*nrowpd)+5]=-nu; 
}


