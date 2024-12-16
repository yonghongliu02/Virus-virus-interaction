
#################COVID&Influenza####################
sipr2 <- function(time, state, pars) {
  
  with(as.list(c(state, pars)), {

    beta_cov <-B1*( b1*(sin(2*pi*(tt-120)/365))+1)
    beta_inf <- B2*(b2*(sin(2*pi*(tt)/365))+1)
    dSS <- n-beta_cov* SS*(IS+II+IP+IR) / N -beta_inf* SS*(SI+II+PI+RI) / N+zeta_cov*RS+zeta_inf*SR-SS*theta
    dIS <-  beta_cov* SS*(IS+II+IP+IR) / N- beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N-gamma_cov*IS+zeta_inf*IR-IS*theta
    dPS <- gamma_cov*IS-rho*PS-beta_inf*(1-sigma2)*PS*(SI+II+PI+RI) / N+zeta_inf*PR-PS*theta
    dRS <- rho*PS-beta_inf* RS*(SI+II+PI+RI) / N-zeta_cov*RS+zeta_inf*RR-RS*theta
    #####
    dSI <- beta_inf* SS*(SI+II+PI+RI) / N-beta_cov*(1-sigma1)* SI*(IS+II+IP+IR) / N-gamma_inf*SI+zeta_cov*RI-SI*theta
    dII <- beta_cov*(1-sigma1)* SI*(IS+II+IP+IR) / N + beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N - gamma_cov*II - gamma_inf*II-II*theta
    dPI <- beta_inf*(1-sigma2)* PS*(SI+II+PI+RI) / N+gamma_cov*II-PI*gamma_inf-PI*rho-PI*theta
    dRI <- beta_inf*RS*(SI+II+PI+RI)/N+PI*rho-RI*gamma_inf-zeta_cov*RI-RI*theta
    #####
    dSP <- gamma_inf*SI-rho*SP-beta_cov*(1-sigma1)*SP*(IS+II+IP+IR) / N+zeta_cov*RP-SP*theta
    dIP <- II*gamma_inf+beta_cov*(1-sigma1)*SP*(IS+II+IP+IR) / N -gamma_cov*IP-rho*IP-IP*theta
    dPP <- IP*gamma_cov+PI*gamma_inf-2*rho*PP-PP*theta
    dRP <- RI*gamma_inf+PP*rho-RP*rho-zeta_cov*RP-RP*theta
    
    #####
    dSR <- rho*SP-beta_cov* SR*(IS+II+IP+IR) / N -zeta_inf*SR+zeta_cov*RR-SR*theta
    dIR <-beta_cov* SR*(IS+II+IP+IR) / N+IP*rho-IR*gamma_cov -zeta_inf*IR-IR*theta
    dPR <- gamma_cov*IR+PP*rho-PR*rho-zeta_inf*PR-PR*theta
    dRR <- (PR+RP)*rho-zeta_cov*RR-zeta_inf*RR-RR*theta
    #####
    dN <- dSS+ dIS+dPS+ dRS+ dSI+dII+dPI+dRI+dSP+dIP+dPP+dRP+dSR+dIR+dPR+dRR
    dV1 <- beta_cov* SS*(IS+II+IP+IR) / N+beta_cov*(1-sigma1)* SI*(IS+II+IP+IR) / N +beta_cov*(1-sigma1)*SP*(IS+II+IP+IR) / N +beta_cov* SR*(IS+II+IP+IR)/N
    dV2 <- beta_inf* SS*(SI+II+PI+RI) / N+beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N +beta_inf*(1-sigma2)*PS*(SI+II+PI+RI) / N +beta_inf* RS*(SI+II+PI+RI)/N
    dtt <- 1
    return(list(c(dSS, dIS,dPS, dRS, dSI,dII,dPI,dRI,dSP,dIP,dPP,dRP,dSR,dIR,dPR,dRR,dN,dV1,dV2,dtt)))
  })
}
library(deSolve)

N=20000000
is <- 1
si <- 1
ss <-N-1-1
Time <- 1:(365*13+1)
ps <- rs <- ii <- Pi <- ri<- sp <- ip <-pp <- rp <-  sr <-ir <- pr <-  rr <- 0
init <- c(SS=ss, IS=is,PS=ps, RS=rs, SI=si,II=ii,PI=Pi,RI=ri,SP=sp,IP=ip,PP=pp,RP=rp,SR=sr,IR=ir,PR=pr,RR=rr,N=N,V1=0,V2=0,tt=1)
sigma1_plan <- seq(0,1,0.2)
sigma2_plan <- 0
R0_plan <-c(1.2,1.8,2.4)
re <- data.frame(cases =NA,virus=NA,sigma1=NA,sigma2=NA,Rt=NA)

for(i in 1:length(sigma1_plan ) ){
  for(t in 1:length(sigma2_plan )){
    for(z in 1:length(R0_plan )){
      N=20000000
      is <- 0
      si <- 1
      ss <-N-1
      Time <- 1:(365*13+1)
      ps <- rs <- ii <- Pi <- ri<- sp <- ip <-pp <- rp <-  sr <-ir <- pr <-  rr <- 0
      init <- c(SS=ss, IS=is,PS=ps, RS=rs, SI=si,II=ii,PI=Pi,RI=ri,SP=sp,IP=ip,PP=pp,RP=rp,SR=sr,IR=ir,PR=pr,RR=rr,N=N,V1=0,V2=0,tt=1)
      pars <- c(b1=0.1,
                b2=0.1,
                B1=R0_plan[z]/5,
                B2=1.55/4,
                sigma1=sigma1_plan[i],
                sigma2=sigma2_plan[t],
                rho=1/56,
                gamma_cov=1/5,
                gamma_inf=1/4,
                zeta_cov=1/365*2,
                zeta_inf=1/365,
                theta=1/70/365,
                n=1/70/365*N
      )
      res.sir <- as.data.frame(lsoda(y = init, times = Time, func = sipr2, parms = pars ))
      res.sir <- res.sir[365*13+1,]
      init <- c(SS=res.sir$SS, IS=1,PS=res.sir$PS, RS=res.sir$PS, 
                SI=res.sir$SI,II=res.sir$II,PI=res.sir$PI,RI=res.sir$RI,
                SP=res.sir$SP,IP=res.sir$IR,PP=res.sir$PP,RP=res.sir$RP,
                SR=res.sir$SR,IR=res.sir$IR,PR=res.sir$PR,RR=res.sir$RR,N=N,V1=0,V2=res.sir$V2,tt=1)
      
      res.sir <- as.data.frame(lsoda(y = init, times = Time, func = sipr2, parms = pars ))
      
      V1 <- round(diff(res.sir$V1))
      V2 <- round(diff(res.sir$V2))
      V <- c(V1,V2)
      label1 <- sigma1_plan[i]#paste("相互作用系数",i,sep="=")
      label2 <- paste("Rt",R0_plan[z],sep="=")
      data <- data.frame(cases=V,
                         virus=c(rep("virus 1",length(V1)),rep("virus 2",length(V1))),
                         sigma1=rep(label1,length(V1)*2),
                         sigma2=rep(sigma2_plan[t],length(V1)*2),
                         Rt=rep(R0_plan[z],length(V1)*2) )
      re <- rbind(re,data)
      
    }
  }}
re <- re[-1,]
tim <- 1:(365*3)

re$time <- rep((1:(365*13)),length(sigma1_plan)*length(R0_plan)*2)
re2 <- re[which(re$time>(365*10)),]




#################RSV&Influenza###########
sipr3 <- function(time, state, pars) {
  
  with(as.list(c(state, pars)), {
    
    beta_rsv <-B1*( b1*(sin(2*pi*(tt)/365))+1)
    beta_inf <- B2*(b2*(sin(2*pi*(tt)/365))+1)
    dSS <- n-beta_rsv* SS*(IS+II+IP+IR) / N -beta_inf* SS*(SI+II+PI+RI) / N+zeta_rsv*RS+zeta_inf*SR-SS*theta
    dIS <-  beta_rsv* SS*(IS+II+IP+IR) / N- beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N-gamma_rsv*IS+zeta_inf*IR-IS*theta
    dPS <- gamma_rsv*IS-rho*PS-beta_inf*(1-sigma2)*PS*(SI+II+PI+RI) / N+zeta_inf*PR-PS*theta
    dRS <- rho*PS-beta_inf* RS*(SI+II+PI+RI) / N-zeta_rsv*RS+zeta_inf*RR-RS*theta
    #####
    dSI <- beta_inf* SS*(SI+II+PI+RI) / N-beta_rsv*(1-sigma1)* SI*(IS+II+IP+IR) / N-gamma_inf*SI+zeta_rsv*RI-SI*theta
    dII <- beta_rsv*(1-sigma1)* SI*(IS+II+IP+IR) / N + beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N - gamma_rsv*II - gamma_inf*II-II*theta
    dPI <- beta_inf*(1-sigma2)* PS*(SI+II+PI+RI) / N+gamma_rsv*II-PI*gamma_inf-PI*rho-PI*theta
    dRI <- beta_inf*RS*(SI+II+PI+RI)/N+PI*rho-RI*gamma_inf-zeta_rsv*RI-RI*theta
    #####
    dSP <- gamma_inf*SI-rho*SP-beta_rsv*(1-sigma1)*SP*(IS+II+IP+IR) / N+zeta_rsv*RP-SP*theta
    dIP <- II*gamma_inf+beta_rsv*(1-sigma1)*SP*(IS+II+IP+IR) / N -gamma_rsv*IP-rho*IP-IP*theta
    dPP <- IP*gamma_rsv+PI*gamma_inf-2*rho*PP-PP*theta
    dRP <- RI*gamma_inf+PP*rho-RP*rho-zeta_rsv*RP-RP*theta
    
    #####
    dSR <- rho*SP-beta_rsv* SR*(IS+II+IP+IR) / N -zeta_inf*SR+zeta_rsv*RR-SR*theta
    dIR <-beta_rsv* SR*(IS+II+IP+IR) / N+IP*rho-IR*gamma_rsv -zeta_inf*IR-IR*theta
    dPR <- gamma_rsv*IR+PP*rho-PR*rho-zeta_inf*PR-PR*theta
    dRR <- (PR+RP)*rho-zeta_rsv*RR-zeta_inf*RR-RR*theta
    #####
    dN <- dSS+ dIS+dPS+ dRS+ dSI+dII+dPI+dRI+dSP+dIP+dPP+dRP+dSR+dIR+dPR+dRR
    dV1 <- beta_rsv* SS*(IS+II+IP+IR) / N+beta_rsv*(1-sigma1)* SI*(IS+II+IP+IR) / N +beta_rsv*(1-sigma1)*SP*(IS+II+IP+IR) / N +beta_rsv* SR*(IS+II+IP+IR)/N
    dV2 <- beta_inf* SS*(SI+II+PI+RI) / N+beta_inf*(1-sigma2)* IS*(SI+II+PI+RI) / N +beta_inf*(1-sigma2)*PS*(SI+II+PI+RI) / N +beta_inf* RS*(SI+II+PI+RI)/N
    dtt <- 1
    return(list(c(dSS, dIS,dPS, dRS, dSI,dII,dPI,dRI,dSP,dIP,dPP,dRP,dSR,dIR,dPR,dRR,dN,dV1,dV2,dtt)))
  })
}
library(deSolve)
#################rsvID&Influenza
N=20000000
is <- 1
si <- 1
ss <-N-1-1
Time <- 1:(365*13+1)
ps <- rs <- ii <- Pi <- ri<- sp <- ip <-pp <- rp <-  sr <-ir <- pr <-  rr <- 0
init <- c(SS=ss, IS=is,PS=ps, RS=rs, SI=si,II=ii,PI=Pi,RI=ri,SP=sp,IP=ip,PP=pp,RP=rp,SR=sr,IR=ir,PR=pr,RR=rr,N=N,V1=0,V2=0,tt=1)
sigma1_plan <- seq(-2.5,0,0.5)
sigma2_plan <-seq(-2.5,0,0.5)

re <- data.frame(cases =NA,virus=NA,sigma1=NA,sigma2=NA)

for(i in 1:length(sigma1_plan ) ){
  for(t in 1:length(sigma2_plan )){
    
    N=20000000
    is <- 1
    si <- 1
    ss <-N-1-1
    Time <- 1:(365*13+1)
    ps <- rs <- ii <- Pi <- ri<- sp <- ip <-pp <- rp <-  sr <-ir <- pr <-  rr <- 0
    init <- c(SS=ss, IS=is,PS=ps, RS=rs, SI=si,II=ii,PI=Pi,RI=ri,SP=sp,IP=ip,PP=pp,RP=rp,SR=sr,IR=ir,PR=pr,RR=rr,N=N,V1=0,V2=0,tt=1)
    pars <- c(b1=0.1,
              b2=0.1,
              B1=2.5/9,
              B2=1.55/4,
              sigma1=sigma1_plan[i],
              sigma2=sigma2_plan[t],
              rho=1/56,
              gamma_rsv=1/9,
              gamma_inf=1/4,
              zeta_rsv=1/200,
              zeta_inf=1/365,
              theta=1/70/365,
              n=1/70/365*N
    )
    res.sir <- as.data.frame(lsoda(y = init, times = Time, func = sipr3, parms = pars ))
    
    V1 <- round(diff(res.sir$V1))
    V2 <- round(diff(res.sir$V2))
    V <- c(V1,V2)
    label1 <- sigma1_plan[i]#paste("相互作用系数",i,sep="=")
    label2 <- sigma2_plan[t]
    data <- data.frame(cases=V,
                       virus=c(rep("virus 1",length(V1)),rep("virus 2",length(V1))),
                       sigma1=rep(label1,length(V1)*2),
                       sigma2=rep(sigma2_plan[t],length(V1)*2))
    re <- rbind(re,data)
    
    
  }}
re <- re[-1,]
tim <- 1:(365*3)

re$time <- rep((1:(365*13)),length(sigma1_plan)*length(sigma2_plan)*2)
re3 <- re[which(re$time>(365*10)),]
#################plot#########

tiff(file = "fig.tif", width = 6000, height = 3000, res = 500,compression = "lzw")
par(mfrow=c(3,4),mar=c(4,4,1,1))
for(a in 1:3){
  A <- c("A","B","C")
  R0_plan <-c(1.2,1.8,2.4)
  i=R0_plan[a]
  plot(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==0&re2$Rt==i),"cases"],ylim=c(0,130000),type="l",xlab="",ylab="Cases",col=x[1],lwd=1,xaxt = "n",main=paste0("SARS-CoV-2: Rt=",i),
       cex.main=0.9,cex.lab=1  )
  x <- c("#a8ddb5","#7bccc4","#4eb3d3","#2b8cbe","#0868ac","#084081")
  mtext("Year", side = 1, line = 2,cex=0.7)
  lines(tim,re2[which(re2$virus=="virus 2"&re2$sigma1==0&re2$Rt==i),"cases"],lty=2,lwd=1)
  lines(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==0.2&re2$Rt==i),"cases"],col=x[2],lwd=1)
  lines(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==0.4&re2$Rt==i),"cases"],col=x[3],lwd=1)
  lines(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==(0.2+0.4)&re2$Rt==i),"cases"],col=x[4],lwd=1)
  lines(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==0.8&re2$Rt==i),"cases"],col=x[5],lwd=1)
  lines(tim,re2[which(re2$virus=="virus 1"&re2$sigma1==1.0&re2$Rt==i),"cases"],col=x[6],lwd=1)
  xmaks <- c("11","12","13","14")
  axis.Date(1, at = c(0,366,365*2+1,365*3+1), format = xmaks,cex.axis=1 )
  legend(-400,160000,legend=A[a],lwd=0.5,cex=1.2,col="white",
         bty="n",ncol=1, xpd = TRUE)
  
}
plot(0:100, 0:100, type="n", axes=FALSE, xlab="", ylab="")  
legend(-30,100,legend=c("IAV","SARS-CoV-2: σ =  0.0","SARS-CoV-2: σ = -0.2",
                        "SARS-CoV-2: σ = -0.4","SARS-CoV-2: σ = -0.6",
                        "SARS-CoV-2: σ = -0.8","SARS-CoV-2: σ = -1.0"),lwd=0.5,cex=1,
       col=c("black",x[1:3],x[4:6]),lty=c(2,1,1,1,1,1,1,1),bty="n",ncol=1, xpd = TRUE)

A <- c("D","E","F")
for(a in 1:3){  
  sigma1_plan <- c(-0.5,-1,-2)
  i=sigma1_plan[a]
  x <- c("#fdbb84","#fc8d59","#ef6548","#d7301f","#7f0000")
  plot(tim,re3[which(re3$virus=="virus 2"&re3$sigma2==0&re3$sigma1==i),"cases"],ylim=c(0,320000),type="l",xlab="",ylab="Cases",col=x[1],lwd=1,xaxt = "n",
       cex.main=0.9,cex.lab=1,main=paste0("RSV:σ =",-i))
  mtext("Year", side = 1, line = 2,cex=0.7)
  lines(tim,re3[which(re3$virus=="virus 1"&re3$sigma2==-0.5&re3$sigma1==i),"cases"],lty=2,lwd=1)
  lines(tim,re3[which(re3$virus=="virus 2"&re3$sigma2==-0.5&re3$sigma1==i),"cases"],col=x[2],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 2"&re3$sigma2==-1&re3$sigma1==i),"cases"],col=x[3],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 2"&re3$sigma2==-1.5&re3$sigma1==i),"cases"],col=x[4],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 2"&re3$sigma2==-2&re3$sigma1==i),"cases"],col=x[5],lwd=1)
  
  xmaks <- c("11","12","13","14")
  axis.Date(1, at = c(0,366,365*2+1,365*3+1), format = xmaks,cex.axis=1 )
  legend(-400,395000,legend=A[a],lwd=0.5,cex=1.2,col="white",
         bty="n",ncol=1, xpd = TRUE)
}
plot(0:100, 0:100, type="n", axes=FALSE, xlab="", ylab="") 
legend(-30,100,legend=c("RSV","IAV: σ= 0.0","IAV: σ= 0.5",
                        "IAV: σ= 1.5","IAV: σ= 2.0",
                        "IAV: σ= 2.5"),lwd=0.5,cex=1,
       col=c("black",x),lty=c(2,1,1,1,1,1),bty="n",ncol=1, xpd = TRUE)



A <- c("G","H","I")

for(a in 1:3){
  x <- c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c")
  i=sigma1_plan[a]
  plot(tim,re3[which(re3$virus=="virus 1"&re3$sigma1==0&re3$sigma2==i),"cases"],ylim=c(0,280000),type="l",xlab="",ylab="Cases",col=x[1],lwd=1,xaxt = "n",
       cex.main=0.9,cex.lab=1,main=paste0("IAV:σ =",-i))
  mtext("Year", side = 1, line = 2,cex=0.7)
  lines(tim,re3[which(re3$virus=="virus 2"&re3$sigma1==0&re3$sigma2==i),"cases"],lty=2,lwd=1)
  lines(tim,re3[which(re3$virus=="virus 1"&re3$sigma1==-0.5&re3$sigma2==i),"cases"],col=x[2],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 1"&re3$sigma1==-1&re3$sigma2==i),"cases"],col=x[3],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 1"&re3$sigma1==-1.5&re3$sigma2==i),"cases"],col=x[4],lwd=1)
  lines(tim,re3[which(re3$virus=="virus 1"&re3$sigma1==-2&re3$sigma2==i),"cases"],col=x[5],lwd=1)
  
  xmaks <- c("11","12","13","14")
  axis.Date(1, at = c(0,366,365*2+1,365*3+1), format = xmaks ,cex.axis=1)
  legend(-400,346000,legend=A[a],lwd=0.5,cex=1.2,col="white",
         bty="n",ncol=1, xpd = TRUE)
  
  
  
}
plot(0:100, 0:100, type="n", axes=FALSE, xlab="", ylab="") 
legend(-30,100,legend=c("IAV","RSV: σ= 0.0","RSV: σ= 0.5",
                        "RSV: σ= 1.5","RSV: σ= 2.0",
                        "RSV: σ= 2.5"),lwd=0.5,cex=1,
       col=c("black",x),lty=c(2,1,1,1,1,1),bty="n",ncol=1, xpd = TRUE)
dev.off()

