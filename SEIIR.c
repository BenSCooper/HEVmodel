
//  seir model
//
//  Created by BSC on 9/5/14.
//  Copyright 2014. All rights reserved.
//
// Implements deterministic SEIIR model to be used with R 
// (two I comparments with different transmissibilities and mean sojourn times)
// To compile from R: system("R CMD SHLIB SEIIR.c")
// To load from R (in Unix): dyn.load("SEIIR.so")
// To call from R use deSolve package and see details in help on using compiled code.


#include </Library/Frameworks/R.framework/Versions/2.12/Resources/include/R.h> 

static double parms[5];

#define beta1 parms[0]
#define beta2 parms[1]
#define gamma parms[2]
#define rho1 parms[3]
#define rho2 parms[4]

/* initializer */
void initmod(void (* odeparms)(int * , double *))
{
    int N=5;
    odeparms(&N, parms);
}

/* derivatives and one output variable 
y[0] is S (susceptibles)
y[1] is E (latently infected)
y[2] is I1 (infectious compartment 1)
y[3] is I2 (infectious compartment 2)
y[4] is R (recovered and immune)
y[5] is CumInc (cumulative incidence)
*/
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if(ip[0]<1) error("nout should be at least 1");
    double N;
    N= y[0]+y[1]+y[2]+y[3]+y[4];
    
    ydot[0] = -beta1*y[0]*y[2]/N  -beta2*y[0]*y[3]/N;
    ydot[1] =  beta1*y[0]*y[2]/N +beta2*y[0]*y[3]/N - gamma*y[1];  
    ydot[2] = gamma*y[1] - rho1*y[2];
    ydot[3] = rho1*y[2] - rho2*y[3];  
    ydot[4] =  rho2*y[3]; 
    ydot[5] = gamma*y[1] ;
    
    yout[0]=y[0]+y[1]+y[2]+y[3]+y[4];
}

/* Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    double N;
    
    N= y[0]+y[1]+y[2]+y[3]+y[4];
    pd[0]=-beta1*y[2]/N -beta2*y[3]/N ;
    pd[1]=beta1*y[2]/N + beta2*y[3]/N;
    pd[2]=0;
    pd[3]=0;
    pd[4]=0;
    pd[5]=0;

    pd[(*nrowpd)]=0;
    pd[(*nrowpd)+1]=-gamma;
    pd[(*nrowpd)+2]=gamma;
    pd[(*nrowpd)+3]=0;
    pd[(*nrowpd)+4]=0;
    pd[(*nrowpd)+5]=gamma;
  
    pd[2*(*nrowpd)]=-beta1*y[0]/N;
    pd[2*(*nrowpd)+1]=beta1*y[0]/N;
    pd[2*(*nrowpd)+2]=-rho1;
    pd[2*(*nrowpd)+3]=rho1;
    pd[2*(*nrowpd)+4]=0;
    pd[2*(*nrowpd)+5]=0;
    
    pd[3*(*nrowpd)]=-beta2*y[0]/N;
    pd[3*(*nrowpd)+1]=beta2*y[0]/N;  
    pd[3*(*nrowpd)+2]=0;  
    pd[3*(*nrowpd)+3]=-rho2;  
    pd[3*(*nrowpd)+4]=rho2; 
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


