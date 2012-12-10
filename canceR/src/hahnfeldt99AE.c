#include <R.h>
static double parms[8];
#define lambda1 parms[0]
#define lambda2 parms[1]
#define b parms[2]
#define d parms[3]
#define eA parms[4]
#define clrA parms[5]
#define eE parms[6]
#define clrE parms[7]

/* initializer  */
void initPhilAE(void (* odeparms)(int *, double *))
{
    int N=8;
    odeparms(&N, parms);
}

/* Derivatives */
void derivsPhilAE (int *neq, double *t, double *y, double *ydot)
{   double V,K,gA,gE;
    V=y[0];K=y[1];gA=y[2];gE=y[3];
    ydot[0] = -lambda1*V*log(V/K);
    ydot[1] = -lambda2*K+b*V-d*K*pow(V,.667)-eA*K*gA-eE*K*gE;
    ydot[2] = -clrA*gA;
    ydot[3] = -clrE*gE;
}

