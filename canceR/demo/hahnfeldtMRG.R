rm(list=ls())
library(tidyverse)  
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
out=mod%>%mrgsim(end = 25, delta = 0.1) 
plot(out,xlab="Days")

e=ev(M="C",ID=1,time=10,amt=100,cmt=3)+ev(M="C",ID=1,time=20,amt=100,cmt=4)
out=mod%>%ev(e)%>%mrgsim(end = 25, delta = 0.1) 
plot(out,xlab="Days")


# dataset
g=NULL               # time in days
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

# (out=mod%>%mrgsim(end = 30, delta = 0.1))
# mod%>%data_set(d)%>%mrgsim(end = 30, delta = 0.1) # ignores delta

#driving events
ev1=ev(Trt="Control",ID=1, time=0, cmt=1,amt=0)
ev2=ev(Trt="Angiostatin",ID=2, time=0,ii=1,addl=100, cmt=3,amt=20)
ev3=ev(Trt="Endostatin",ID=3, time=0,ii=1,addl=100, cmt=4,amt=20)
ev4=ev(Trt="Elow",ID=4, time=0,ii=1,addl=100, cmt=4,amt=4)
ev5a=ev(Trt="A+E",ID=5, time=0,ii=1,addl=100, cmt=3,amt=20)
ev5b=ev(Trt="A+E",ID=5, time=0,ii=1,addl=100, cmt=4,amt=20)
(evnt=ev1+ev2+ev3+ev4+ev5a+ev5b)
out=mod%>%ev(evnt)%>%mrgsim(end = 25, delta = 1) 
D=out@data%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
g+geom_line(data=D)
ggsave("~/tmp/smooth.png",width=5,height=4)

(pars=c( lambda1=.192,b=5.85,d=0.00873,eA=0.15,eE=0.66))
PH=function(pars) {
  evnt@data=merge(evnt@data,data.frame(t(pars)))
  (mod <- mod %>% ev(evnt))
  (out=mod%>%mrgsim(end = 100, delta = 1) )
  D=as.data.frame(as_tibble(out))
  D[!duplicated(D[1:2]),]
}

D=PH(pars)%>%select(ID,time,V)%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
head(D)
sd=2
D$V=D$V+rnorm(dim(D)[1],sd=sd)
D%>%ggplot(aes(x=time,y=V))+facet_wrap(Trt~.,scales = "free")+geom_line(size=1)+tc(11)+sbb
ggsave("~/tmp/simulated.png",width=5,height=4)


library(FME)
D=D%>%select(time,V)
PHcost <- function (pars) {
  out=PH(pars)%>%select(time,V)
  modCost(model = out, obs = D)
  return(modCost(model = out, obs = D))
}

(Fit <- modFit(f = PHcost, p = 1.3*pars))
summary(Fit)

# ##### logs dont help
# PHcostLog <- function(lpars)  PHcost(c(exp(lpars)))
# (FitL <- modFit(f = PHcostLog, p = log(2*pars)))
# exp(coef(FitL))
# summary(FitL)

(pars=c( lambda1=.192,b=5.85,d=0.00873,eE=0.66))
(Fit <- modFit(f = PHcost, p = 1.3*pars))
summary(Fit)


(pars=c( lambda1=.192,b=5.85,d=0.00873))
(Fit <- modFit(f = PHcost, p = 1.3*pars))
summary(Fit)

(pars=c( lambda1=.192,b=5.85,d=0.00873))
(Fit <- modFit(f = PHcost, p = 2*pars))
summary(Fit)

######### stop here on going into readme front page


Sfun <- sensFun(PHcost, pars)
summary(Sfun)



str(Sfun)
Sfun=
plot(Sfun[Sfun$x!=0,], which = c("V"), xlab="Time (Days)",ylab="Sensitivity", lwd = 2,legpos="bottomright")



PHcostReal <- function (pars) {
  out=PH(pars)%>%select(time,V)
  return(modCost(model = out, obs = (d%>%select(time,V))))
}



(Fit <- modFit(f = PHcostReal, p = pars))
PHcostLog <- function(lpars)  PHcost(c(exp(lpars)))
(FitL <- modFit(f = PHcostLog, p = log(pars)))
coef(Fit)
exp(coef(FitL))
summary(Fit)

ini=PH(pars)%>%select(ID,time,V)%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
final=PH(pars = coef(Fit))%>%select(ID,time,V)%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
final=PH(pars = exp(coef(FitL)))%>%select(ID,time,V)%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))

g+geom_line(data=ini%>%filter(time<25))+geom_line(data=final%>%filter(time<25),col="red")


head(ini)
d%>%ggplot(aes(x))
pdf("~/Results/myelo/fit5pars1p5mrg.pdf",width=5,height=5)
plot(d[,1:2], xlab = "time", ylab = "Q",col="gray")
lines(ini$time, ini$Q, col="red")
lines(final$time, final$Q,col="blue")
legend("topright", c("data", "initial", "fitted"),
       lty = c(NA,1,1),col=c("gray","red","blue"), pch = c(1, NA, NA))
dev.off()
par(mfrow = c(1, 1))




Sfun <- sensFun(PHcost, pars)
summary(Sfun)

pdf("~/Results/myelo/FMEphil.pdf",width=5,height=5)
plot(Sfun, which = c("V"), xlab="Time (Days)",ylab="Sensitivity", lwd = 2,legpos="bottomright")
dev.off()

pdf("~/Results/myelo/sensPairsPhil.pdf",width=5,height=5)
pairs(Sfun, which = c("V"))
dev.off()

ident <- collin(Sfun)
ident
tail(ident,7)
plot(ident, log = "y")
