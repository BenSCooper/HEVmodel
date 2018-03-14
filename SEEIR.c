
//  seeir model (extension of seir to allow latent period to be the sum of exponential rvs.
//
//  Created by BSC on  16/5/14.
//  Copyright 2014. All rights reserved.
//
// Implements deterministic SEEIR model to be used with R 
// To compile from R: system("R CMD SHLIB SEEIR.c")
// To load from R (in Unix): dyn.load("SEEIR.so")
// To call from R use deSolve package and see details in help on using compiled code.


#include </Library/Frameworks/R.framework/Versions/2.12/Resources/include/R.h> 

static double parms[4];

#define beta parms[0]
#define gamma1 parms[1]
#define gamma2 parms[2]
#define rho parms[3]

/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=4;
    odeparms(&N, parms);
}

/* derivatives and one output variable 
y[0] is S (susceptibles)
y[1] is E1 (latently infected compartment 1)
y[2] is E2 (latently infected compartment 2)
y[3] is I (infectious)
y[4] is R (recovered and immune)
y[5] is CumInc (cumulative incidence)
*/
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1");
    double N;
    N= y[0]+y[1]+y[2]+y[3]+y[4];
    
    ydot[0] = -beta*y[0]*y[3]/N;
    ydot[1] =  beta*y[0]*y[3]/N - gamma1*y[1];
    ydot[2] =  gamma1*y[1] - gamma2*y[2];  
    ydot[3] = gamma2*y[2] - rho*y[3];  
    ydot[4] =  rho*y[3]; 
    ydot[5] = gamma1*y[1];
    
    yout[0]=y[0]+y[1]+y[2]+y[3]+y[4];
}

/* Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    double N;
    
    N= y[0]+y[1]+y[2]+y[3] + y[4];
    pd[0]=-beta*y[3]/N;
    pd[1]=beta*y[3]/N;
    pd[2]=0;
    pd[3]=0;
    pd[4]=0;
    pd[5]=0;

    pd[(*nrowpd)]=0;
    pd[(*nrowpd)+1]=-gamma1;
    pd[(*nrowpd)+2]=gamma1;
    pd[(*nrowpd)+3]=0;
    pd[(*nrowpd)+4]=0;
    pd[(*nrowpd)+5]=gamma1;

    pd[2*(*nrowpd)]=0;
    pd[2*(*nrowpd)+1]=0;
    pd[2*(*nrowpd)+2]=-gamma2;
    pd[2*(*nrowpd)+3]=gamma2;
    pd[2*(*nrowpd)+4]=0;
    pd[2*(*nrowpd)+5]=0;

    pd[3*(*nrowpd)]=-beta*y[0]/N;
    pd[3*(*nrowpd)+1]=beta*y[0]/N;  
    pd[3*(*nrowpd)+2]=0;  
    pd[3*(*nrowpd)+3]=-rho;  
    pd[3*(*nrowpd)+4]=rho; 
    pd[3*(*nrowpd)+5]=0; 

    pd[4*(*nrowpd)]=0;
    pd[4*(*nrowpd)+1]=0;  
    pd[4*(*nrowpd)+2]=0;  
    pd[4*(*nrowpd)+3]=0;  
    pd[4*(*nrowpd)+4]=0; 
    pd[4*(*nrowpd)+5]=0; 

    pd[5*(*nrowpd)]=0;
    pd[5*(*nrowpd)+1]=0;  
    pd[5*(*nrowpd)+2]=0;  
    pd[5*(*nrowpd)+3]=0;  
    pd[5*(*nrowpd)+4]=0; 	
    pd[5*(*nrowpd)+5]=0; 
}


