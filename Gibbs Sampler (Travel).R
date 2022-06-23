##########################################################################
#Use Data to Populate y.t-vector and R.E.-matrix
##########################################################################

library(mvtnorm)

y.t=fitdata$Surveyed.Travel.Cost
R.t=matrix(0,nrow=length(fitdata$stategroup),ncol=max(fitdata$stategroup))

for(i in 1:length(fitdata$stategroup))
{
  R.t[i,fitdata$stategroup[i]]=1  
}

##########################################################################
#Create Storage Vectors, Initialize
##########################################################################

I=10000
n=length(y.t)
m=ncol(R.t)

alpha.t=matrix(NA,nrow=I,ncol=m)
mu.t=matrix(NA,nrow=I,ncol=1)
phi.t=matrix(NA,nrow=I,ncol=m)
psi.t=matrix(NA,nrow=I,ncol=1)
cost.t=matrix(NA,nrow=nrow(extrap),ncol=I)

phi.t[1,]=1
psi.t[1]=1
mu.t[1]=900
alpha.t[1,]=0
cost.t[,1]=extrap[,1]

##########################################################################
#Gibbs Sampler
##########################################################################

for(i in 2:I)
{
  v=2
  list = fitdata$stategroup
  a = (v+1)/2
  b = .5*phi.t[i-1,list]*(y.t-mu.t[i-1]-alpha.t[i-1,list])^2+v/2
  g = rgamma(n,shape=a,rate=b)
  
  #Create Covariance Matrix
  
  S = g*phi.t[i-1,list] 
  
  #Updating mu.t
  
  V=1/sum(S)
  E=V*crossprod(S,y.t-alpha.t[i-1,list])
  
  temp1=-1
  while(temp1<0)
  {
    temp1=rnorm(1,mean=E,sd=sqrt(V))
  }
  mu.t[i]=temp1
  
  #Updating alpha.t.t
  
  V=solve(t(R.t)%*%diag(S)%*%R.t+diag(psi.t[i-1],m))
  E=V%*%(t(R.t)%*%diag(S)%*%(y.t-mu.t[i]))
  
  c=-1
  while(c<0)
  {
    temp2=t(rmvnorm(1,mean=E,sigma=V))
    c=min(mu.t[i]+temp2)
  }
  alpha.t[i,]=temp2
  
  #Updating phi.t
  
  for(j in 1:m)
  {
    a=length(y.t[fitdata$stategroup==j])/2-1
    b=.5*t(g[list==j]*(y.t[list==j]-mu.t[i]-alpha.t[i,j]))%*%(y.t[list==j]-mu.t[i]-alpha.t[i,j])
    phi.t[i,j]=rgamma(1,shape=a,rate=b)  
  }
  
  #Updating psi.t
  
  a=m/2-0.5
  b=.5*t(alpha.t[i,])%*%(alpha.t[i,])
  psi.t[i]=rgamma(1,shape=a,rate=b)
  
  temp1 = as.data.frame(cbind(seq(1,7,by=1),mu.t[i]+c(alpha.t[i,],0)))
  colnames(temp1)[1] = "stategroup"
  temp2 = join(extrap,temp1,by="stategroup",type="left")
  cost.t[,i] = temp2[,4]
  
  if(i%%1000==0)
  {
    print(i)
  }
  
}

burn=4000

# par(mfrow=c(3,3))
# plot(mu.t[burn:I],cex=.25,type="l",ylab="mu.t",main="mu.t")
# plot(sqrt(1/phi.t[burn:I,1]),cex=.25,type="l",ylab="phi.t 1",main="phi.t 1")
# plot(sqrt(1/phi.t[burn:I,2]),cex=.25,type="l",ylab="phi.t 2",main="phi.t 2")
# plot(sqrt(1/phi.t[burn:I,3]),cex=.25,type="l",ylab="phi.t 3",main="phi.t 3")
# plot(sqrt(1/phi.t[burn:I,4]),cex=.25,type="l",ylab="phi.t 4",main="phi.t 4")
# plot(sqrt(1/phi.t[burn:I,5]),cex=.25,type="l",ylab="phi.t 5",main="phi.t 5")
# plot(sqrt(1/phi.t[burn:I,6]),cex=.25,type="l",ylab="phi.t 6",main="phi.t 6")
# plot(sqrt(1/psi.t[burn:I]),cex=.25,type="l",ylab="psi.t",main="psi.t")
# plot(alpha.t[burn:I,1],cex=.25,type="l",ylab="alpha.t 1",main="alpha.t 1")
# plot(alpha.t[burn:I,2],cex=.25,type="l",ylab="alpha.t 2",main="alpha.t 2")
# plot(alpha.t[burn:I,3],cex=.25,type="l",ylab="alpha.t 3",main="alpha.t 3")
# plot(alpha.t[burn:I,4],cex=.25,type="l",ylab="alpha.t 4",main="alpha.t 4")
# plot(alpha.t[burn:I,5],cex=.25,type="l",ylab="alpha.t 5",main="alpha.t 5")
# plot(alpha.t[burn:I,6],cex=.25,type="l",ylab="alpha.t 6",main="alpha.t 6")
# 
# state1=mu.t[burn:I]+alpha.t[burn:I,1]
# state2=mu.t[burn:I]+alpha.t[burn:I,2]
# state3=mu.t[burn:I]+alpha.t[burn:I,3]
# state4=mu.t[burn:I]+alpha.t[burn:I,4]
# state5=mu.t[burn:I]+alpha.t[burn:I,5]
# state6=mu.t[burn:I]+alpha.t[burn:I,6]
# 
# mean(state1);quantile(state1,c(.025,.975))
# mean(state2);quantile(state2,c(.025,.975))
# mean(state3);quantile(state3,c(.025,.975))
# mean(state4);quantile(state4,c(.025,.975))
# mean(state5);quantile(state5,c(.025,.975))
# mean(state6);quantile(state6,c(.025,.975))
# mean(mu.t[burn:I]);quantile(mu.t[burn:I],c(.025,.975))
