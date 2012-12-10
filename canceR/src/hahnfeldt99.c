#include <R.h>
static double parms[6];
#define lambda1 parms[0]
#define lambda2 parms[1]
#define b parms[2]
#define d parms[3]
#define e parms[4]
#define clr parms[5]

/* initializer  */
void initPhil(void (* odeparms)(int *, double *))
{
    int N=6;
    odeparms(&N, parms);
}

/* Derivatives */
void derivsPhil (int *neq, double *t, double *y, double *ydot)
{   double V,K,g;
    V=y[0];K=y[1];g=y[2];
    ydot[0] = -lambda1*V*log(V/K);
    ydot[1] = -lambda2*K+b*V-d*K*pow(V,.667)-e*K*g;
    ydot[2] = -clr*g;
}

