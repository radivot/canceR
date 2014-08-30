# first growth with no treatment, from help page
library(canceR)
out1=lsoda(y=c(V=180,K= 625,g=0),times=seq(0,20,1),hahnfeldt99,  
		parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,e=0,clr=0),rtol=1e-4, atol= rep(1e-4,3)) 
plot(out1) 

# now do Angiostatin and Endostatin together, daily for 24 days
hahnfeldt99AE<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dV = -lambda1*V*log(V/K)
				dK = -lambda2*K+b*V-d*K*V^.667-eA*K*gA-eE*K*gE
				dgA = -clrA*gA
				dgE = -clrE*gE
				return(list(c(dV,dK,dgA,dgE)))
			})
}

library(deSolve)
y0<-c(V=180,K= 625,gA=20,gE=20) # angio and endo high dose combined
outs=NULL
eventdat <- data.frame(var = rep(c("gA", "gE"),25),
		time = rep(0:24,each=2) ,
		value = rep(20,50),
		method = c("add"))
eventdat

out=lsoda(y=y0,times=seq(0,30,.2),hahnfeldt99AE, parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,
					eA=0.15,clrA=0.38,eE=0.66,clrE=1.7), rtol=1e-4, atol= rep(1e-4,4),events = list(data = eventdat))
plot(out)


out1=lsoda(y=c(V=180,K= 625,g=0),times=seq(0,20,1),hahnfeldt99,  
		parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,e=0,clr=0),rtol=1e-4, atol= rep(1e-4,3)) 
plot(out1) 

# now check to see if DLL's are woring the same way
(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "canceR"),
					paste("canceR",.Platform$dynlib.ext,sep="")))
dyn.load(f)
#dlls <- getLoadedDLLs()
#getDLLRegisteredRoutines(dlls[["canceR"]])
#getDLLRegisteredRoutines("canceR")
#getNativeSymbolInfo("canceR")

out1=ode(y=c(V=180,K= 625,g=0),times=seq(0,20,1),  
		parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,e=0,clr=0),rtol=1e-4, atol= rep(1e-4,3),
		func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )
plot(out1) 

out=lsoda(y=y0,times=seq(0,30,.2),parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,
eA=0.15,clrA=0.38,eE=0.66,clrE=1.7), rtol=1e-4, atol= rep(1e-4,4),events = list(data = eventdat),
				func = "derivsPhilAE", dllname = "canceR",initfunc = "initPhilAE" )
plot(out)

######################################################################################
# Warning: the following code is old to the extent that it uses neither bbmle nor events
#########################################################################################
(fsrc=file.path(system.file("demo/hahnfeldtDemoFuncs.R", package = "canceR")))
source(fsrc)

g=NULL
g$dC=data.frame( days=c(0,4,7,10,13,16,19), size=c(192,728,2319,3590,5852,8591,10438))
g$dT=data.frame( days=c(0,4,7,10,13 ), size=c(170,381,875,1363,1750 ))
g$dA=data.frame( days=c(0,4,7,10,13), size=c(176,246,265,136,207))
g$dE=data.frame( days=c(0,4,7,10 ), size=c(186,157,135,76 ))
g$dEL=data.frame( days=c(0,4,7,10,13,16), size=c(291,429,859,2420,3120,4100)) # endo at 4mg/kg
g$dAE=data.frame( days=c(0,4,7,10,13,16,19,22,25), size=c(309,180,88,79,29,6,0,0,0))
g$nData=dim(g$dC)[1]+dim(g$dT)[1]+dim(g$dA)[1]+dim(g$dE)[1]+dim(g$dEL)[1]+dim(g$dAE)[1]
parICs=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,eT=1.3,clrT=10.1,eA=0.15,clrA=0.38,eE=0.66,clrE=1.7,V0=180,V0valid=300,K0=625)
g$params=data.frame(initial=parICs,final=parICs,opt=TRUE,CI95prct="not fitted")
g


g=simData(g)
# first take a look at the initial parameter value situation
plotg(g)  # mg/kg time courses
plotVK(g) # Tumor volume setpoint and process variable time courses
g$params

g$params[c("lambda2"),"opt"]=FALSE   
#g$params[c("lambda2","eT","clrT","eA","clrA","eE","clrE"),"opt"]=FALSE  # only fit control parameters and V0valid
p0=t(g$params[g$params[,"opt"],"initial",drop=F])
p0=unlist(as.data.frame(p0))
p0=sapply(p0,log)
g$optNames<-names(p0)
g

t0 = Sys.time();  
opt<-optim(p0,fopt,hessian=TRUE,model=g,control=list(trace=F,maxit=5000))
difftime(Sys.time(),t0,units="mins")[[1]]

opar<-opt$par
Hess=opt$hessian
value=opt$value
sg<-sqrt(value/(g$nData-length(p0)))
if (abs(det(Hess))>1e-8) {
	sig=sg*sqrt(diag(solve(Hess/2)))
	upper=signif(opar+1.96*sig,3)
	lower=signif(opar-1.96*sig,3)
	g$CI=exp(cbind(lower,upper))
}  
g$params[g$optNames,"final"]=signif(exp(opar),3)  
g$params[g$optNames,"CI95prct"]=paste("(",signif(g$CI[,"lower"],3),", ",signif(g$CI[,"upper"],3),")",sep="")  
g$params

# now look at final parameter estimate situation
g=simData(g)
plotg(g)  # mg/kg time courses
plotVK(g) # Tumor volume setpoint and process variable time courses
dyn.unload(f)
