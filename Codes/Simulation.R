 ####################################################################################################
################# Parameters                                                        ################
####################################################################################################
n       <- 1000;                       #Number of simulations
T       <- 120;                        #Data points per time serie
N       <- 1000;                       #Number of companies in a given portfolio
shift   <- 3;                          #Expresses how many time steps the business cycle is leading on a default event
Z       <- matrix(nrow=n,ncol=T+shift);#BASIS Business Cycle no structural break
Z_b     <- matrix(nrow=n,ncol=T+shift);#ACTUAL Business Cycle shares structural break
Z_      <- matrix(nrow=n,ncol=T+shift);#Similar Business Cycle but less correlated than the basis (no structural break)
d       <- 0.7;                        #Dependancy of the similar Business Cycle
R       <- array(dim=c(n,T,N));        #Asset returns
r       <- 0.02;                        #Asset correlation
a       <- 0.9;                        #Autocorrelation;
t_break <- 60;                          #structural break at t=t_break
ADF_s   <- matrix(nrow=n,ncol=T);      #ADF over time per simulation
PD_s    <- matrix(nrow=n,ncol=T);      #PD (Expected Default Frequency) over time per simulation
beta    <- 1:n;                        #Slope without shift
beta_   <- 1:n;                        #Slope without shift for bad BC
beta_s  <- 1:n;                        #Slope with shift
beta_s_ <- 1:n;                        #Slope with shift for bad BC
beta_b  <- 1:n;                        #Slope with with shift + structural break
beta_b_ <- 1:n;                        #Slope with with shift + structural break, calibration DOLS

rho     <- 1:n;                        #Rho (Asset Correlation) without shift
rho_    <- 1:n;                        #Rho (Asset Correlation) without shift for bad BC
rho_s   <- 1:n;                        #Rho (Asset Correlation) with shift
rho_s_  <- 1:n;                        #Rho (Asset Correlation) with shift for bad BC
rho_b   <- 1:n;                        #Rho (Asset Correlation) with shift + structural break
rho_b_  <- 1:n;                        #Rho (Asset Correlation) with shift + structural break, calibration DOLS

gamma    <- 1:n;                        #Slope without shift
gamma_   <- 1:n;                        #Slope without shift for bad BC
gamma_s  <- 1:n;                        #Slope with shift
gamma_s_ <- 1:n;                        #Slope with shift for bad BC
gamma_b  <- 1:n;                        #Slope with shift + structural break

alpha   <- 1:n;                        #Alpha (PiTness factor) without shift
alpha_  <- 1:n;                        #Alpha (PiTness factor) without shift for bad BC
alpha_s <- 1:n;                        #Alpha (PiTness factor) with shift
alpha_s_<- 1:n;                        #Alpha (PiTness factor) with shift for bad BC
alpha_b <- 1:n;                        #Alpha (PiTness factor) with shift+ structural break

PiT     <- 0.8;                        #PiTness factor of the expected default rate (PD)
R2      <- matrix(nrow=n,ncol=5);      #Coefficient of determination per model per simulation
B       <- 1:n;                        #Beste linear model for a given simulation
e       <- matrix(nrow=n,ncol=5);      #Check autocorrelation models
################# Risk segments                                                     ################
nb      <- 3;                                                    #Nb of stratum
ADF     <- c(0.012,0.056,0.1);                                   #AVG Default rate per stratum
CV       <- array(dim=c(nb,T));                                #Critical value for each stratum per simulation
#simulate time dependent CVs
cycleLength =matrix(data= 60,nrow=1,ncol=nb) #initial cycle lenghts of CVs
for (t in 1:T){
  if (t==1){
    CV[,t] = qnorm(ADF) #initial
  }
  else if (t>=t_break){ #sintroduce structural break
    CV[, t] = qnorm(ADF+0.05) 
  }else{
    CV[,t] =CV[,(t-1)]
  }
}
I       <- c(round(N/nb/2),round(2*N/nb),N);                     #Index used to identify stratum
ADF_P   <- sum((I[2:nb]-I[1:(nb-1)])*ADF[2:nb]/N)+I[1]/N*ADF[1]; #Approximation of the Portfolio ADF
####################################################################################################
################# Simulation                                                        ################
####################################################################################################
Z [,1]  <- rnorm(n);                                                  #Starting value @t=0
Z_[,1]  <- sqrt(1-d) *rnorm(n) + sqrt(d) *Z[,1];                      #Starting value @t=0
for (t in 2:(T+shift)){
  Z [,t] <- sqrt(1-a) *rnorm(n) + sqrt(a) *Z[,(t-1)];                  #Autocorrelated noise
  Z_[,t] <- sqrt(1-d) *rnorm(n) + sqrt(d) *Z[,t];                      #Business Cycle that is dependent on the actual
  if (t>=t_break){
    Z_b[,t] = Z [,t]+(qnorm(ADF_P)-qnorm(ADF_P+0.05))
  }else{
    Z_b[,t] = Z [,t]
  }
}
Z        <- t(sapply(1:n,function(i) qnorm(rank(Z[i,])/(T+shift+1))));#Calculate the empirical CDF and transform to a standard normal
Z_b      <- t(sapply(1:n,function(i) qnorm(rank(Z_b[i,])/(T+shift+1))));#Calculate the empirical CDF and transform to a standard normal

for(i in 1:N) {
  R[,,i]<- sapply(1:T,         function(t)
    sqrt(1-r) *rnorm(n) + sqrt(r) *Z_b[,(t+shift)]);             #Bivariate normal per entity with a given business cycle
}
####################################################################################################
################# Estimate Asset Correlation and Pitness                            ################
####################################################################################################
ADF_s<-t(sapply(matrix(1:n,c(n,1)),  function(x) 
  sapply(matrix(1:T,c(1,T)),  function(t) sum(
    sapply(array(1:N,c(1,1,N)), function(i) 
      R[x,t,i]<=CV[which(I-i>=0)[1],t]))/N)));              #Calculate ADF per simulation and time observation
PD_s <-t(sapply(matrix(1:n,c(n,1)),  function(x)
  sapply(matrix(1:T,c(1,T)),  function(t) 
    pnorm((qnorm(ADF_P) -sqrt(r)*PiT *Z_b[x,(t+shift)])/sqrt(1-r*PiT**2)))));#PD with a PiTness factor  
ADF_n <- apply(ADF_s,c(1,2),qnorm);                                                   #Transform ADF to distance to default.
PD_n  <- apply(PD_s ,c(1,2),qnorm);                                                   #Transform PD to distance to default.
for (x in 1:n) { 
  v         <- which(ADF_s[x,]!=0);                                 #Exclude 0 observations for ADF
  y         <- diff(ADF_n[x,v])
  
  x1        <- diff(Z [x,v])
  beta[x]   <- sum(y*x1)/sum(x1**2) #cov(  ADF_n[x,v],        Z [x,v])/var(Z [x,v]);      #Calculate slope without drift
  e2        <- (y - beta[x]*x1)              # residuals
  R2[x,1]   <- sum((beta[x]*x1)**2)/sum(y**2) #1-var(e2)/var(y);     #Calculate R square
  rho[x]    <- beta[x]**2/(1+beta[x]**2);                           #Calculate asset correlation without shift
  e2        <- e2 -mean(e2)
  e[x,1]    <- cor(e2[1:(length(v)-2)]**2,e2[2:(length(v)-1)]**2)              #Calculate autocorrelation residuals
  
  x2        <- diff(Z_ [x,v])
  beta_[x]  <- sum(y*x2)/sum(x2**2) #cov(  ADF_n[x,v],         Z_[x,v])/var(Z_[x,v]);     #Calculate slope without drift for bad BC
  e2_       <- (y - beta_[x]*x2)              #Square residuals
  R2[x,2]   <- sum((beta_[x]*x2)**2)/sum(y**2) #1-var(e2_)/var(y);      #Calculate R square
  rho_[x]   <- beta_[x]**2/(1+beta_[x]**2);                         #Calculate asset correlation without shift for bad BC
  e2_       <- e2_ -mean(e2_)
  e[x,2]    <- cor(e2_[1:(length(v)-2)]**2,e2_[2:(length(v)-1)]**2)             #Calculate autocorrelation residuals
  
  x3        <- diff(Z [x,v+shift])
  beta_s[x] <- sum(y*x3)/sum(x3**2);#Calculate slope with shift (Note different from simulation)
  e2_s      <- (y - beta_s[x]*x3)            # residuals
  R2[x,3]   <- sum((beta_s[x]*x3)**2)/sum(y**2) #1-var(e2_s)/var(y);     #Calculate R square
  rho_s[x]  <- beta_s[x]**2/(1+beta_s[x]**2);                       #Calculate slope with shift
  e2_s      <- e2_s -mean(e2_s)
  e[x,3]    <- cor(e2_s[1:(length(v)-2)]**2,e2_s[2:(length(v)-1)]**2)   #Calculate autocorrelation residuals
  
  x4        <- diff(Z_ [x,v+shift])
  beta_s_[x]<- sum(y*x4)/sum(x4**2); #cov(ADF_n[x,v],Z_[x,(v+shift)])/var(Z_[x,(v+shift)]);#Calculate slope with shift (Note different from simulation) for bad BC
  e2_s_     <- (y - beta_s_[x]*x4)              # residuals
  R2[x,4]   <- sum((beta_s_[x]*x4)**2)/sum(y**2) #1-var(e2_s_)/var(y)     #Calculate R square
  rho_s_[x] <- beta_s_[x]**2/(1+beta_s_[x]**2);                     #Calculate slope with shift for bad BC
  e2_s_     <- e2_s_ -mean(e2_s_)
  e[x,4]    <- cor(e2_s_[1:(length(v)-2)]**2,e2_s_[2:(length(v)-1)]**2)  #Calculate autocorrelation residuals
  
  x5        <- diff(Z_b [x,v+shift])
  beta_b[x] <- sum(y*x5)/sum(x5**2) #cov(  ADF_n[x,v],        Z [x,v])/var(Z [x,v]);      #Calculate slope without drift
  e2_b      <- (y - beta_b[x]*x5)              #Square residuals
  R2[x,5]   <- sum((beta_b[x]*x5)**2)/sum(y**2) #1-var(e2_b)/var(y)     #Calculate R square
  rho_b[x]  <- beta_b[x]**2/(1+beta_b[x]**2);                         #Calculate asset correlation without shift for bad BC
  
  LL <- function(b1){
    re = y - diff(Z_b [x,v+shift]) *b1
    n  = length(re)
    re = lm(re[2:n]~re[1:(n-1)])$residuals
    return(-sum(log(dnorm(re,0,1))))
  }
  
  beta_b_[x]<-  nlm(LL,1)$estimate;
  rho_b_[x] <- beta_b_[x]**2/(1+beta_b_[x]**2);
  e2_b      <- e2_b -mean(e2_b)
  e[x,5]    <- cor(e2_b[1:(T-2)]**2,e2_b[2:(T-1)]**2)
  
  B[x]      <- which(R2[x,]==max(R2[x,]));                                  #Select best linear model
 
   ############# Calibration PiTness ##########################################################################################################################
  y_a        <- diff(PD_n[x,v])
  gamma[x]   <- sum(y_a*diff(Z [x,v]))/sum(diff(Z [x,v])**2)      #Calculate slope without drift
  alpha[x]   <- min(sqrt(gamma[x]**2/(1+gamma[x]**2)/rho[x]),1);    #Alpha (PiTness factor) without shift best BC
  
  gamma_[x]   <- sum(y_a*diff(Z_ [x,v]))/sum(diff(Z_ [x,v])**2)      #Calculate slope without drift
  alpha_[x]  <-  min(sqrt(gamma_[x]**2/(1+gamma_[x]**2)/rho_[x]),1);    #Alpha (PiTness factor) without shift bad BC
  
  gamma_s[x]   <- sum(y_a*diff(Z [x,v+shift]))/sum(diff(Z [x,v+shift])**2)      #Calculate slope without drift
  alpha_s[x] <-  min(sqrt(gamma_s[x]**2/(1+gamma_s[x]**2)/rho_s[x]),1);   #Alpha (PiTness factor) with shift best BC
  
  gamma_s_[x]   <-sum(y_a*diff(Z_ [x,v+shift]))/sum(diff(Z_ [x,v+shift])**2)      #Calculate slope without drift
  alpha_s_[x]<- min(sqrt(gamma_s_[x]**2/(1+gamma_s_[x]**2)/rho_s_[x]),1);    #Alpha (PiTness factor) with shift bad BC
  
  gamma_b[x]   <- sum(y_a*diff(Z_b [x,v+shift]))/sum(diff(Z_b [x,v+shift])**2)       #Calculate slope without drift
  alpha_b[x]   <- min(sqrt(gamma_b[x]**2/(1+gamma_b[x]**2)/rho_b[x]),1);    #Alpha (PiTness factor) with shift bad BC
  
  }
####################################################################################################
################# Plot Results                                                      ################
####################################################################################################
par(mar=c(4,3.9,0.5,0.4))
boxplot(ADF_s,xlab="Time", ylab="ADF", cex.lab = 0.8,cex.axis = 0.8); 
grid()

labs   = c(expression(atop(tau==0,Z[2]) ), 
           expression(atop(tau==0, Z[3])),
           expression(atop(tau==3, Z[2])), 
           expression(atop(tau==3, Z[3])),
           expression(atop(tau==3, Z[1])))
par(mar=c(2.5,3.9,0.4,0.4))
boxplot(cbind(rho,rho_,rho_s,rho_s_,rho_b), 
        xaxt="n",
        ylab=expression({rho}),
        cex.lab = 1.2,cex.axis = 0.8)
axis(side = 1, at =1:5, as.expression(labs), tick=TRUE, padj=0.3,  cex.axis=0.8)
grid()

par(mar=c(2.5,3.9,0.4,0.4))
boxplot(cbind(alpha,alpha_,alpha_s,alpha_s_,alpha_b), ylim=c(0,1.1),
        xaxt="n",
        ylab=expression({alpha}),
        cex.lab = 1.2,cex.axis = 0.8)
axis(side = 1, at =1:5, as.expression(labs), tick=TRUE, padj=0.3,  cex.axis=0.8)
grid()

plot(B)

boxplot(e, 
        xaxt="n",
        ylab="residual autocorrelation",
        cex.lab = 0.8,cex.axis = 0.8)
axis(side = 1, at =1:5, as.expression(labs), tick=TRUE, padj=0.3,  cex.axis=0.8)
grid()

boxplot(R2, 
        xaxt="n",
        ylab="R squared",
        cex.lab = 0.8,cex.axis = 0.8)
axis(side = 1, at =1:5, as.expression(labs), tick=TRUE, padj=0.3,  cex.axis=0.8)
grid()

par(mar=c(2.6,4.1,0.4,2))
plot(ADF_s[1,], xlab="Time", ylab="ADF", cex.lab = 0.8,cex.axis = 0.8)
par(new=T)
plot(Z[1,4:123],type="l", xlab=NA, ylab=NA,xaxt="n",xaxt="n",yaxt="n",cex.lab = 0.8,cex.axis = 0.8)
axis(side = 4,cex.lab = 0.8,cex.axis = 0.8)
legend("topright",legend=c("Z"),lty=c(1),pch=c(NA),merge=FALSE, cex=0.8)
legend("topleft",legend=c("ADF"),lty=c(0),pch=c(1),merge=FALSE, cex=0.8)
grid()

par(mar=c(2.6,4.1,0.4,2))
plot(ADF_n[1,], xlab="Time", ylab=expression({N^-1}(ADF)),cex.lab = 0.8,cex.axis = 0.8)
par(new=T)
plot(Z[1,4:123],type="l", xlab=NA, ylab=NA,xaxt="n",xaxt="n",yaxt="n",cex.lab = 0.8,cex.axis = 0.8)
axis(side = 4,cex.lab = 0.8,cex.axis = 0.8)
legend("topright",legend=c("Z"),lty=c(1),pch=c(NA),merge=FALSE, cex=0.8)
legend("topleft",legend=expression({N^-1}(ADF)),lty=c(0),pch=c(1),merge=FALSE, cex=0.8)
grid()

bla<-sqrt(r)*Z_b[,(shift+1):(T+shift)]+sqrt(1-r)*ADF_n
boxplot(apply(bla,c(1,2),pnorm))

labs   = c(expression(atop(tau==3, Z[1]: OLS)),
           expression(atop(tau==3, Z[1]: MLE)))
par(mar=c(2.5,3.9,0.4,0.4))
boxplot(cbind(rho_b,rho_b_), 
        xaxt="n",
        ylab=expression({rho}),
        cex.lab = 1.3,cex.axis = 0.8)
axis(side = 1, at =1:2, as.expression(labs), tick=TRUE, padj=0.3,  cex.axis=0.8)
grid()
# 
# plot(diff(Z_b [x,v+shift]),beta_b[x]*diff(Z_b [x,v+shift]), type='l', title="break")
# points(diff(Z_b [x,v+shift]), y)
# 
# plot(diff(Z [x,v]), y)
# lines(diff(Z [x,v]),beta[x]*diff(Z [x,v]))
