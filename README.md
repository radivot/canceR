canceR
======

An R package housing dynamical system models of non-myeloid cancers.
To install it use:  
```
devtools::install_github("radivot/canceR",subdir="canceR")
```

There is currently only one model in this package, the dynamic  carrying capacity model of Hahnfeldt et al. 

# Tumor Development under Angiogenic Signaling: A Dynamical Theory of Tumor Growth, Treatment Response, and Postvascular Dormancy

[Hahnfeldt et al  *Cancer Research* **59** 4770-4775 (1999)](https://www.ncbi.nlm.nih.gov/pubmed/10519381) provide a model that includes as  state variables, the tumor volume V, its dynamic carrying capacity K, and the concentrations gA and gE of angiostatin and endostatin. In R their model is 
```
library(canceR)
hahnfeldt99AE<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dV = -lambda1*V*log(V/K)
				dK = -lambda2*K+b*V-d*K*V^.667-eA*K*gA-eE*K*gE
				dgA = -clrA*gA
				dgE = -clrE*gE
				return(list(c(dV,dK,dgA,dgE)))
			})
}
```

and in equivalent but faster C it is 
```
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
```

Using Metrum Research Group's mrgsolve, such C code is automatically generated and compiled using this neat R code.
```
library(mrgsolve)
code='
$PARAM lambda1=.192,lambda2=0,b=5.85,d=0.00873,eA=0.15,eE=0.66,clrA=.38,clrE=1.7
$INIT V=180,K= 625,gA=0,gE=0 
$ODE 
dxdt_V = -lambda1*V*log(V/K);
dxdt_K = -lambda2*K+b*V-d*K*pow(V,.667)-eA*K*gA-eE*K*gE;
dxdt_gA = -clrA*gA;
dxdt_gE = -clrE*gE;
'
mod <- mread("phil", "~/tmp", code)
mod%>%mrgsim(end = 30, delta = 0.1)%>%plot(xlab="Days") #no drugs
e=ev(time=10,amt=100,cmt=3)+ev(time=20,amt=100,cmt=4)
e # adds 100 mg/kg to Angiostatin (state 3) at 10 days and endostatin (state 4) at 20 days
mod%>%ev(e)%>%mrgsim(end = 30, delta = 0.1)%>%plot(xlab="Days")
```
which generates these plots.

![](docs/noDrug.png)
![](docs/day10n20.png)


The data provided with this  model is seen by this code 
```
g=NULL              
g$dC=tibble(Trt="Control",ID=1,time=c(0,4,7,10,13,16,19), V=c(192,728,2319,3590,5852,8591,10438),gA=0,gE=0)
g$dA=tibble(Trt="Angiostatin",ID=2, time=c(0,4,7,10,13), V=c(176,246,265,136,207),gA=20,gE=0)
g$dE=tibble(Trt="Endostatin",ID=3, time=c(0,4,7,10 ), V=c(186,157,135,76 ),gA=0,gE=20)
g$dEL=tibble(Trt="Elow",ID=4, time=c(0,4,7,10,13,16), V=c(291,429,859,2420,3120,4100),gA=0,gE=4) # endo at 4mg/kg
g$dAE=tibble(Trt="A+E",ID=5, time=c(0,4,7,10,13,16,19,22,25), V=c(309,180,88,79,29,6,0,0,0),gA=20,gE=20)
(d=bind_rows(g)%>%mutate(Trt=as_factor(Trt)))
sbb=theme(strip.background=element_blank())
ltp=theme(legend.position="top",legend.title=element_blank())
tc=function(sz) theme_classic(base_size=sz)
(g=d%>%ggplot(aes(x=time,y=V))+facet_wrap(~Trt)+geom_point()+scale_y_log10()+tc(13)+sbb+ltp+xlab("Days"))
ggsave("~/tmp/data.png",width=5,height=4)
```

![](docs/data.png)
which shows that Angiostatin and Endostatin alone at 20mg/kg are static, low doses of 
Endostatin (4mg/kg) is simimar to controls, and the combination is most effective. To add smooth curves of the model fit provided in the paper  we need to simulate the model with doses given daily out to 25 days to reach all of the data points. 

```
ev1=ev(Trt="Control",ID=1, time=0, cmt=1,amt=0)
ev2=ev(Trt="Angiostatin",ID=2, time=0,ii=1,addl=24, cmt=3,amt=20)
ev3=ev(Trt="Endostatin",ID=3, time=0,ii=1,addl=24, cmt=4,amt=20)
ev4=ev(Trt="Elow",ID=4, time=0,ii=1,addl=24, cmt=4,amt=4)
ev5a=ev(Trt="A+E",ID=5, time=0,ii=1,addl=24, cmt=3,amt=20)
ev5b=ev(Trt="A+E",ID=5, time=0,ii=1,addl=24, cmt=4,amt=20)
(evnt=ev1+ev2+ev3+ev4+ev5a+ev5b)
out=mod%>%ev(evnt)%>%mrgsim(end = 25, delta = 0.1) 
D=out@data%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
g+geom_line(data=D)
ggsave("~/tmp/smooth.png",width=5,height=4)
```

![](docs/smooth.png)


We see that the model provides a decent fit. 




