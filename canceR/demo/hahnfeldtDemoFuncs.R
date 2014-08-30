updateSSE<-function(dX,oldSSE) {
#  oldSSE+sum((dX$EV-dX$size)^2)  #too much emphasis on control data
	oldSSE+sum((sqrt(dX$EV)-sqrt(dX$size))^2)  # sqrt transform: still not great since K0 still less than V0
#  oldSSE+sum((dX$EV-dX$size)^2)/mean(dX$size)^2 # normalize: this led to K0 that are too small to believe
#  oldSSE+sum((log(dX$EV)-log(dX$size))^2) # problem is that some volumes are zero
}

simData<-function(g){
	SSE=0
	Pframe=as.data.frame(t(g$params[,"final",drop=F]))
	attach(Pframe)
	y0<-c(V=V0,K= K0,g=0)
	out1=ode(y=y0,times=seq(0,20,1),parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,e=0,clr=0), 
			rtol=1e-4, atol= rep(1e-4,3),
			func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )

	g$outC=as.data.frame(out1)
	g$dC$EV=g$outC[match(g$dC$days,g$outC$time),"V"]
	SSE=updateSSE(g$dC,0)
	
	y0<-c(V=V0,K= K0,g=30)
	outs=NULL
	for (i in 0:6) {
		out1=lsoda(y=y0,times=seq(2*i,2*(i+1),.2), parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,e=eT,clr=clrT), 
				rtol=1e-4, atol= rep(1e-4,3),
				func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )
		y0=out1[nrow(out1),2:4]+c(0,0,30)
		outs=rbind(outs,out1)
	}
	g$outT=as.data.frame(outs)
	g$dT$EV=g$outT[match(g$dT$days,g$outT$time),"V"]
	SSE=updateSSE(g$dT,SSE)
	
	outs=NULL
	y0<-c(V=V0,K= K0,g=20)
	for (i in 0:13) {
		out1=lsoda(y=y0,times=seq(i,i+1,.2),parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,e=eA,clr=clrA), 
				rtol=1e-4, atol= rep(1e-4,3),
				func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )
		y0=out1[nrow(out1),2:4]+c(0,0,20)
		outs=rbind(outs,out1)
	}
	g$outA=as.data.frame(outs)
	g$dA$EV=g$outA[match(g$dA$days,g$outA$time),"V"]
	SSE=updateSSE(g$dA,SSE)
	
	y0<-c(V=V0,K= K0,g=20)
	outs=NULL
	for (i in 0:9) {
		out1=lsoda(y=y0,times=seq(i,i+1,.2), parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,e=eE,clr=clrE), 
				rtol=1e-4, atol= rep(1e-4,3),
				func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )
		y0=out1[nrow(out1),2:4]+c(0,0,20)
		outs=rbind(outs,out1)
	}
	g$outE=as.data.frame(outs)
	g$dE$EV=g$outE[match(g$dE$days,g$outE$time),"V"]
	SSE=updateSSE(g$dE,SSE)
	
	y0<-c(V=V0valid,K= K0,g=4)  # endo low dose not enough to control tumor
	outs=NULL
	for (i in 0:15) {
		out1=lsoda(y=y0,times=seq(i,i+1,.2), parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,e=eE,clr=clrE), 
				rtol=1e-4, atol= rep(1e-4,3),
				func = "derivsPhil", dllname = "canceR",initfunc = "initPhil" )
		y0=out1[nrow(out1),2:4]+c(0,0,4)
		outs=rbind(outs,out1)
	}
	g$outEL=as.data.frame(outs)
	g$dEL$EV=g$outEL[match(g$dEL$days,g$outEL$time),"V"]
	SSE=updateSSE(g$dEL,SSE)
	
	y0<-c(V=V0valid,K= K0,gA=20,gE=20) # angio and endo high dose combined
	outs=NULL
	for (i in 0:24) {
		out1=lsoda(y=y0,times=seq(i,i+1,.2), parms=c(lambda1=lambda1,lambda2=0,b=b,d=d,
#            eA=0.15,clrA=0.38,eE=0.66,clrE=1.7), rtol=1e-4, atol= rep(1e-4,4),dll="philModAE")
						eA=eA,clrA=clrA,eE=eE,clrE=clrE), rtol=1e-4, atol= rep(1e-4,4),
				func = "derivsPhilAE", dllname = "canceR",initfunc = "initPhilAE" )
		y0=out1[nrow(out1),2:5]+c(0,0,20,20)
		outs=rbind(outs,out1)
	}
	g$outAE=as.data.frame(outs)
	g$dAE$EV=g$outAE[match(g$dAE$days,g$outAE$time),"V"]
	g$SSE=updateSSE(g$dAE,SSE)
	
	detach(Pframe)
	g
} ### 



fopt <- function(pars,model) {
	pars=exp(pars)
	model$params[names(pars),"final"]=pars  
	model=simData(model)
	return(model$SSE)   
}

plotVK<-function(g) {
	windows()
	par(mfrow=c(2,3))
	with(g$outC, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="Control",ylim=c(0,12000));
				lines(time,K,type="l",col="red")} )
	points(g$dC)
	with(g$outT, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="TNP470 (30mg/kg/2day)",ylim=c(0,2500));
				lines(time,K,type="l",col="red")})
	points(g$dT)
	with(g$outA, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="Angiostatin (20mg/kg/day)",ylim=c(0,500));
				lines(time,K,type="l",col="red")})
	points(g$dA)
	with(g$outE, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="Endostatin (20mg/kg/day)",ylim=c(0,400));
				lines(time,K,type="l",col="red")})
	points(g$dE)
	with(g$outEL, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="Endostatin (4mg/kg/day)",ylim=c(0,5000));
				lines(time,K,type="l",col="red")})
	points(g$dEL)
	with(g$outAE, {plot(time,V,type="l",ylab="Tumor Size (uL)",xlab="Time (days)",main="Angio+Endo (20mg/kg/day)",ylim=c(0,350));
				lines(time,K,type="l",col="red")})
	points(g$dAE)
}


plotg<-function(g) {
	windows()
	par(mfrow=c(2,3))
	plot.new()
	with(g$outT, plot(time,g,type="l",ylab="drug (mg/kg)",xlab="Time (days)",main="TNP470 (30mg/kg/2day)"))
	with(g$outA, plot(time,g,type="l",ylab="drug (mg/kg)",xlab="Time (days)",main="Angiostatin (20mg/kg/day)"))
	with(g$outE, plot(time,g,type="l",ylab="drug (mg/kg)",xlab="Time (days)",main="Endostatin (20mg/kg/day)"))
	with(g$outEL, plot(time,g,type="l",ylab="drug (mg/kg)",xlab="Time (days)",main="Endostatin (4mg/kg/day)"))
	with(g$outAE, {plot(time,gA,type="l",ylab="drug (mg/kg)",xlab="Time (days)",main="Angio+Endo (20mg/kg/day)",ylim=c(0,65));
				lines(time,gE)})
}

