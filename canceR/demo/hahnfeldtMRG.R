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
d%>%ggplot(aes(x=time,y=V))+facet_wrap(~Trt)+geom_point()+scale_y_log10()+tc(13)+sbb+ltp+xlab("Days")
ggsave("~/tmp/data.png",width=5,height=4)

(out=mod%>%mrgsim(end = 30, delta = 0.1))
mod%>%data_set(d)%>%mrgsim(end = 30, delta = 0.1) # ignores delta

#driving events
ev1=ev(Trt="Control",ID=1, time=0, cmt=1,amt=0)
ev2=ev(Trt="Angiostatin",ID=2, time=0,ii=1,addl=24, cmt=3,amt=20)
ev3=ev(Trt="Endostatin",ID=3, time=0,ii=1,addl=24, cmt=4,amt=20)
ev4=ev(Trt="Elow",ID=4, time=0,ii=1,addl=24, cmt=4,amt=4)
ev5a=ev(Trt="A+E",ID=5, time=0,ii=1,addl=24, cmt=3,amt=20)
ev5b=ev(Trt="A+E",ID=5, time=0,ii=1,addl=24, cmt=4,amt=20)
(evnt=ev1+ev2+ev3+ev4+ev5a+ev5b)
out=mod%>%ev(evnt)%>%mrgsim(end = 25, delta = 0.1) 
D=out@data%>%mutate(Trt=levels(d$Trt)[ID])%>%mutate(Trt=as_factor(Trt))
d%>%ggplot(aes(x=time,y=V))+facet_wrap(~Trt)+geom_point()+scale_y_log10()+tc(13)+sbb+ltp+xlab("Days")+geom_line(data=D)
ggsave("~/tmp/smooth.png",width=5,height=4)


(pars=c( lambda1=.192,b=5.85,d=0.00873,eA=0.15,eE=0.66))
H=function(pars) {
  evnt@data=merge(evnt@data,data.frame(t(pars)))
  (mod <- mod %>% ev(evnt))
  (out=mod%>%mrgsim(end = 25, delta = 1) )
  as.data.frame(as_tibble(out))
}

D=H(pars)%>%select(ID,time,EV=V)
head(D,2)
D=left_join(d,D)
head(D,2)
D%>%ggplot(aes(x=time,y=EV))+facet_wrap(Trt~.,scales = "free")+geom_line(size=1)+tc(14)+sbb

library(FME)

Hcost <- function (pars) {
  D=H(pars)%>%select(ID,time,EV=V)
  head(D,2)
  D=left_join(d,D)
  (D$EV-D$V)^2
}

Sfun <- sensFun(Hcost, pars)
summary(Sfun)
plot(Sfun, which = c("V"), xlab="Time (Days)",ylab="Sensitivity", lwd = 2,legpos="bottomright")

pdf("~/Results/myelo/FMEsensHSCmrg.pdf",width=5,height=5)
dev.off()

pdf("~/Results/myelo/sensPairsMRG.pdf",width=5,height=5)
pairs(Sfun, which = c("Q"))
dev.off()

