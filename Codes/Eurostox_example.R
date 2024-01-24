####################################################################################################
################# Parameters                                                        ################
####################################################################################################
N       <- 5000;                       #Number of companies in a given portfolio
r       <- 0.02;                       #Asset correlation
set.seed(1234);   
################# Risk segments                                                     ################
nb      <- 3;                                                    #Nb of stratum
ADF     <- c(0.012,0.056,0.1);                                   #Default rate per stratum
CV      <- sapply(ADF, qnorm);                                   #Critical value for each stratum
I       <- c(round(N/nb/2),round(2*N/nb),N);                     #Index used to identify stratum
ADF_P   <- sum((I[2:nb]-I[1:(nb-1)])*ADF[2:nb]/N)+I[1]/N*ADF[1]; #Approximation of the Portfolio ADF
####################################################################################################
################# Simulation                                                        ################
####################################################################################################
#install.packages('qrmdata');
library(qrmdata)
data("EURSTOXX")
ZZ      <- EURSTOXX[xts:::endof(EURSTOXX, "months")][,1]
time    <- time(EURSTOXX[xts:::endof(EURSTOXX, "months")])
zn      <- dim(ZZ)[1];
Z       <- t(log(rev(as.numeric(ZZ[1:(zn-12),1]))/rev(as.numeric(ZZ[13:zn,1]))));
Z_      <- sapply(Z,function(x) (x-mean(Z))/sd(Z));
ret     <- Z;
Z_CDF   <- ecdf(Z);                                                                     #Calculate the empirical CDF
Z       <- sapply(Z, function(x) qnorm(Z_CDF(x)));                                      #Output standard normal transformation of the Business Cycle
R       <- sapply(1:length(Z), function(t) sqrt(1-r) *rnorm(N) + sqrt(r) *Z[t]);        #Asset returns
####################################################################################################
################# Estimate Asset Correlation and Pitness                            ################
####################################################################################################
ADF_s<-sapply(matrix(1:length(Z),c(1,length(Z))),  function(t) sum(sapply(array(1:N,c(1,1,N)), function(i) R[i,t]<=CV[which(I-i>=0)[1]]))/N);#Calculate ADF per simulation and time observation
ADF_CDF <- ecdf(c(ADF_s, max(ADF_s)*(1+0.001)));

ADF_n    <- sapply(ADF_s,qnorm);                                                                                                             #Transform ADF to distance to default.
zero     <- which(ADF_s==0)-1;
ADF_n_df <- (ADF_n[2:length(Z)]-ADF_n[1:(length(Z)-1)])[-c(zero,zero+1)];
Z_df     <- (Z[2:length(Z)]-Z[1:(length(Z)-1)])[-c(zero,zero+1)];
Z_df_    <- (Z_[2:length(Z)]-Z_[1:(length(Z)-1)])[-c(zero,zero+1)];
b1       <- cov(ADF_n_df,Z_df)/var(Z_df);
b1_      <- cov(ADF_n_df,Z_df_)/var(Z_df_);
rho      <- b1**2/(1-b1**2);
rho_     <- b1_**2/(1-b1_**2);
####################################################################################################
################# Plot Results                                                      ################
####################################################################################################
####Plotting colors
yellow = '#FFE600'
gray   = '#80808088'
black  = 'black'

xts::plot.xts(ZZ)
par(mar=c(4,4.1,0.4,3.75))
cex_ = 1.3
#plot(1:length(Z),ret,xaxt="n",ylim=c(-0.6,0.6), xlab='Time',ylab="YoY monthly log changes",type="l",cex.lab = 0.8,cex.axis = 0.8, lwd=2) 
#+ lines(ADF_s,type="h", lwd=3, col=gray)
plot(1:length(Z),ADF_s,mgp=c(2,1,0), xaxt="n", col=gray, ylim=c(-0.25,0.25), xlab='Time',ylab="ADF",type="h",cex.lab = cex_,cex.axis = cex_, lwd=3) 
par(new=T)
plot(1:length(Z),ret,xaxt="n",yaxt="n",ylim=c(-4,4), xlab="",ylab="",type="l",cex.lab = 2,cex.axis = 2, lwd=2) + lines(Z,type = "l",lty=2, lwd=1) + lines(Z_,lwd=1)

axis(side = 4, cex.axis=cex_)
mtext(side = 4, line = 2.5, 'Eurostoxx and Z', cex = cex_)
legend("topleft",legend=bquote("ADF with" ~rho == .(r)),bg ="transparent", fill=c(gray),border=c("black"), lty=c(NA),pch=c(NA),lwd=c(NA),merge=TRUE, cex=cex_)

legend("bottomright",bg ="transparent", lty=c(1,2,1),lwd=c(2,1,1), legend=c("Eurostoxx returns - YoY Log changes",as.expression(bquote("Z (CDF) with callibrated "~rho == .(round(rho,digits=3)))),as.expression(bquote("Z (Standard score) with callibrated"~rho == .(round(rho_,digits=3))))),pt.cex=cex_)
grid()

tmp = time(ZZ)
table_out <- data.frame(
  #T = tmp[1:(zn-12)],
  T = 1:length(Z),
  ADF = ADF_s,
  ret = c(ret),
  Z = Z,
  Zscore = Z_
)
write.table(table_out, "clipboard")
time[seq(1,340,12)]
# 
# plot(Z_df,ADF_n_df)
# lines(Z_df, b1*Z_df)
# 
# plot(Z_df_,ADF_n_df)
# lines(Z_df_, b1_*Z_df_)
# 
# I = is.finite(ADF_n_df)
# Z_model = lm(ADF_n_df[I] ~0 +Z_df[I])
# summary(Z_model)
# Zscore_model = lm(ADF_n_df[I] ~0 +Z_df_[I])
# summary(Zscore_model)
# 

#Rubtsov
E <- function(x) mean(x)
V <- function(x) E(x^2)-E(x)^2
S <- function(x) E(x^3)-3*E(x)*V(x)-E(x)^3
ICDF <- function(x) qnorm(x, mean = 0, sd = 1) #, lower.tail= TRUE)

I = is.finite(ADF_n)
rho_petrov = 1/(1+1/V(ADF_n[I]))  #see (5.4) Rubtsov
B_petrov = E(ADF_n[I])*sqrt(1-rho_petrov)  #see (5.4) Rubtsov
ZZ = (-B_petrov+ADF_n[I]*sqrt(1-rho_petrov))/sqrt(rho_petrov)    #see (5.3) Rubtsov
Y = (S(ADF_n[I])/S(ZZ))^(-2/3)  #see (10) addendum to Rubtsov
rho_rub = 1/(1+Y)
