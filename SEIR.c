
//  seir model
//
//  Created by BSC on 23/3/14.
//  Copyright 2014. All rights reserved.
//
// Implements deterministic SEIR model to be used with R 
// To compile from R: system("R CMD SHLIB testSEIR.c")
// To load from R (in Unix): dyn.load("testSEIR.so")
// To call from R use deSolve package and see details in help on using compiled code.


#include </Library/Frameworks/R.framework/Versions/2.12/Resources/include/R.h> 

static double parms[3];

#define beta parms[0]
#define gamma parms[1]
#define rho parms[2]

/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=3;
    odeparms(&N, parms);
}

/* derivatives and one output variable 
y[0] is S (susceptibles)
y[1] is E (latently infected)
y[2] is I (infectious)
y[3] is R (recovered and immune)
y[4] is CumInc (cumulative incidence)
*/
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1");
    double N;
    N= y[0]+y[1]+y[2]+y[3];
    
    ydot[0] = -beta*y[0]*y[2]/N;
    ydot[1] =  beta*y[0]*y[2]/N - gamma*y[1];  
    ydot[2] = gamma*y[1] - rho*y[2];  
    ydot[3] =  rho*y[2]; 
    ydot[4] = gamma*y[1];
    
    yout[0]=y[0]+y[1]+y[2]+y[3]+y[4];
}

/* Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    double N;
    
    N= y[0]+y[1]+y[2]+y[3];
    pd[0]=-beta*y[2]/N;
    pd[1]=beta*y[2]/N;
    pd[2]=0;
    pd[3]=0;
    pd[4]=0;
    pd[(*nrowpd)]=0;
    pd[(*nrowpd)+1]=-gamma;
    pd[(*nrowpd)+2]=gamma;
    pd[(*nrowpd)+3]=0;
    pd[(*nrowpd)+4]=gamma;
    pd[2*(*nrowpd)]=-beta*y[0]/N;
    pd[2*(*nrowpd)+1]=beta*y[0]/N;
    pd[2*(*nrowpd)+2]=-rho;
    pd[2*(*nrowpd)+3]=rho;
    pd[2*(*nrowpd)+4]=0;
    pd[3*(*nrowpd)]=0;
    pd[3*(*nrowpd)+1]=0;  
    pd[3*(*nrowpd)+2]=0;  
    pd[3*(*nrowpd)+3]=0;  
    pd[3*(*nrowpd)+4]=0; 
    pd[4*(*nrowpd)]=0;
    pd[4*(*nrowpd)+1]=0;  
    pd[4*(*nrowpd)+2]=0;  
    pd[4*(*nrowpd)+3]=0;  
    pd[4*(*nrowpd)+4]=0; 
}


