##########################################################################
#Use Data to Populate y.o-vector and R.E.-matrix
##########################################################################
library(mvtnorm)

y.o=fitdata$Surveyed.Other.Cost
R.o=matrix(0,nrow=length(fitdata$fieldgroup),ncol=max(fitdata$fieldgroup))

for(i in 1:length(fitdata$fieldgroup))
{
  R.o[i,fitdata$fieldgroup[i]]=1  
}

##########################################################################
#Create Storage Vectors, Initialize
##########################################################################

I=10000
n=length(y.o)
m=ncol(R.o)

alpha.o=matrix(NA,nrow=I,ncol=m)
mu.o=matrix(NA,nrow=I,ncol=1)
phi.o=matrix(NA,nrow=I,ncol=m)
psi.o=matrix(NA,nrow=I,ncol=1)
cost.o=matrix(NA,nrow=nrow(extrap),ncol=I)

phi.o[1,]=1
psi.o[1]=1
mu.o[1]=900
alpha.o[1,]=0
cost.o[,1]=extrap[,1]

##########################################################################
#Gibbs Sampler
##########################################################################

for(i in 2:I)
{
  v=1
  list = fitdata$fieldgroup
  a = (v+1)/2
  b = .5*phi.o[i-1,list]*(y.o-mu.o[i-1]-alpha.o[i-1,list])^2+v/2
  g = rgamma(n,shape=a,rate=b)
  
  #Create Covariance Matrix
  
  S = g*phi.o[i-1,list] 
    
  #Updating mu.o
  
  V=1/sum(S)
  E=V*crossprod(S,y.o-alpha.o[i-1,list])
  
  temp1=-1
  while(temp1<0)
  {
    temp1=rnorm(1,mean=E,sd=sqrt(V))
  }
  mu.o[i]=temp1
  
  #Updating alpha.o
  
  V=solve(t(R.o)%*%diag(S)%*%R.o+diag(psi.o[i-1],m))
  E=V%*%(t(R.o)%*%diag(S)%*%(y.o-mu.o[i]))
  
  c=-1
  while(c<0)
  {
    temp2=t(rmvnorm(1,mean=E,sigma=V))
    c=min(mu.o[i]+temp2)
  }
  alpha.o[i,]=temp2
  
  #Updating phi.o
  
  for(j in 1:m)
  {
    a=length(y.o[fitdata$fieldgroup==j])/2-1
    b=.5*t(g[list==j]*(y.o[list==j]-mu.o[i]-alpha.o[i,j]))%*%(y.o[list==j]-mu.o[i]-alpha.o[i,j])
    phi.o[i,j]=rgamma(1,shape=a,rate=b)  
  }
  
  #Updating psi.o
  
  a=m/2-0.5
  b=.5*t(alpha.o[i,])%*%(alpha.o[i,])
  psi.o[i]=rgamma(1,shape=a,rate=b)
  
  temp1 = as.data.frame(cbind(seq(1,5,by=1),mu.o[i]+alpha.o[i,]))
  colnames(temp1)[1] = "fieldgroup"
  temp2 = join(extrap,temp1,by="fieldgroup",type="left")
  cost.o[,i] = temp2[,4]
  
  if(i%%1000==0)
  {
  print(i)
  }
  
}

burn=4000

# par(mfrow=c(2,2))
# plot(mu.o[burn:I],cex=.25,type="l",ylab="mu.o",main="mu.o")
# plot(sqrt(1/phi.o[burn:I,1]),cex=.25,type="l",ylab="phi.o 1",main="phi.o 1")
# plot(sqrt(1/phi.o[burn:I,2]),cex=.25,type="l",ylab="phi.o 2",main="phi.o 2")
# plot(sqrt(1/phi.o[burn:I,3]),cex=.25,type="l",ylab="phi.o 3",main="phi.o 3")
# plot(sqrt(1/phi.o[burn:I,4]),cex=.25,type="l",ylab="phi.o 4",main="phi.o 4")
# plot(sqrt(1/phi.o[burn:I,5]),cex=.25,type="l",ylab="phi.o 5",main="phi.o 5")
# plot(sqrt(1/psi.o[burn:I]),cex=.25,type="l",ylab="psi.o",main="psi.o")
# plot(alpha.o[burn:I,1],cex=.25,type="l",ylab="alpha.o 1",main="alpha.o 1")
# plot(alpha.o[burn:I,2],cex=.25,type="l",ylab="alpha.o 2",main="alpha.o 2")
# plot(alpha.o[burn:I,3],cex=.25,type="l",ylab="alpha.o 3",main="alpha.o 3")
# plot(alpha.o[burn:I,4],cex=.25,type="l",ylab="alpha.o 4",main="alpha.o 4")
# plot(alpha.o[burn:I,5],cex=.25,type="l",ylab="alpha.o 5",main="alpha.o 5")
# 
# field1=mu.o[burn:I]+alpha.o[burn:I,1]
# field2=mu.o[burn:I]+alpha.o[burn:I,2]
# field3=mu.o[burn:I]+alpha.o[burn:I,3]
# field4=mu.o[burn:I]+alpha.o[burn:I,4]
# field5=mu.o[burn:I]+alpha.o[burn:I,5]
# 
# mean(field1);quantile(field1,c(.025,.975))
# mean(field2);quantile(field2,c(.025,.975))
# mean(field3);quantile(field3,c(.025,.975))
# mean(field4);quantile(field4,c(.025,.975))
# mean(field5);quantile(field5,c(.025,.975))
# 
# for(i in 1:5)
# {
#   temp = mean(mu.o[burn:I]+alpha.o[burn:I,i])
#   temp2 = quantile(mu.o[burn:I]+alpha.o[burn:I,i],c(.025,.975))
#   print(temp)
#   print(temp2)
# }


