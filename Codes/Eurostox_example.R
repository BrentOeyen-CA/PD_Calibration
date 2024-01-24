####################################################################################################
################# Parameters                                                        ################
####################################################################################################
N       <- 1000;                       #Number of companies in a given portfolio
r       <- 0.08;                       #Asset correlation
set.seed(1234);   
################# Risk segments     ################
nb      <- 3;                                                    #Nb of stratum
ADF     <- c(0.012,0.056,0.1);                                   #Default rate per stratum
CV      <- sapply(ADF, qnorm);                                   #Critical value for each stratum
I       <- c(round(N/nb/2),round(2*N/nb),N);                     #Index used to identify stratum
ADF_P   <- sum((I[2:nb]-I[1:(nb-1)])*ADF[2:nb]/N)+I[1]/N*ADF[1]; #Approximation of the Portfolio ADF
####################################################################################################
################# Simulation                                                        ################
####################################################################################################

list.of.packages <- c("qrmdata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)
lapply(list.of.packages,library,character.only=TRUE)

data("EURSTOXX")
ZZ      <- EURSTOXX[xts:::endof(EURSTOXX, "months")][,1]
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

####EY plotting colors
EYyellow = '#FFE600'
EYgray   = '#808080'
EYblack = 'black'



xts::plot.xts(ZZ)

text  <- paste("ADF with rho: ",r,sep="");
text1 <- paste("Z (CDF) with calibrated rho: ",round(rho,digits=3),sep="");
text2 <- paste("Z (Standard score) with calibrated rho: ",round(rho_,digits=3),sep="");
par(mar = c(5,5,2,5))
plot(1:length(Z),ret,xlab='Time',ylab="YoY Log changes", type="l",col=EYgray, lwd="1", cex.lab = 0.8,cex.axis = 0.8)
par(new=T)
plot(1:length(Z),Z,type="l",lty=1, col=EYblack,lwd="1.5", xlab="",ylab="",xaxt="n",yaxt="n",cex.lab = 0.8,cex.axis = 0.8)
axis(4, cex.lab = 0.8,cex.axis = 0.8)
legend("bottomleft",legend=c("Eurostoxx returns - YoY Log changes"),col=EYgray,lty=c(1),pch=c(NA),lwd=c(1),merge=FALSE,cex=0.8)
legend("topright",legend=c("Z"),col=EYblack,lty=c(1),pch=c(NA),lwd=c(1),merge=FALSE,cex=0.8)
grid()
axis(side = 4)
mtext(side = 4, line = 3, 'Normal variates')
legend("bottomleft",pch=c(NA,16,10,2),lty=c(1,0,0,0),legend=c("Eurostoxx returns - YoY Log changes",text,text1,text2))

#write output
OUTPUT = cbind(c(ret), c(Z), c(ADF_s))

write.excel <- function(x,row.names=FALSE, col.names=TRUE,...){
  write.table(x,"clipboard", sep="\t", row.names=row.names, col.names=col.names,...)
}
write.excel(OUTPUT)

cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
  write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
}
cb(OUTPUT)
