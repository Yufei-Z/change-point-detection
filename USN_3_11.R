getAf = function(S){
    l2 = length(S)
    l1 = 1:l2
    cumsum(S^2) + (1+l1)*(2*l1+1)/(6*l1)*S^2 - 2/l1*S*cumsum((1:l2)*S)
}
getAb = function(S){
    l2 = length(S)
    l1b = 1:l2
    out0 = cumsum(S^2) + (1+l1b)*(2*l1b+1)/(6*l1b)*S^2 - 2/l1b*S*cumsum((1:l2)*S)
    out = rep(0,l2)
    out[1:(l2-1)] = out0[(l2-1):1]
    out
}
getTf = function(x,epsilon = 0.1){
    n = length(x)
    S = cumsum(x)
    xbar = S/(1:n)

    Tf = rep(NA,n)
    D0 = L0 = R0 = matrix(NA,nr=n,nc=n)
    for(l2 in 3:(n-1)){
        D = (xbar[1:l2]-xbar[l2])*(1:l2)/sqrt(l2)
        A = (getAf(S[1:l2]) + getAb(cumsum(x[l2:1])))/l2^2
        Tf[l2] = max((D^2/A)[2:(l2-1)])
        D0[2:(l2-1),l2] = D[2:(l2-1)]
        L0[2:(l2-1),l2] = (getAf(S[1:l2])/l2^2)[2:(l2-1)]
        R0[2:(l2-1),l2] = (getAb(cumsum(x[l2:1]))/l2^2)[2:(l2-1)]
    }
    max(Tf[(floor(epsilon*n)+1):ceiling((1-epsilon)*n)]) #list(D=D0,L=L0,R=R0)
}
getT = function(x,e = 0.1){
    getTf(x,epsilon=e) + getTf(rev(x),epsilon=e)
}
m_1 = rep(NA,10000)
m_2 = rep(NA,10000)
m_3 = rep(NA,10000)
for( i in 1:10000){
  x = arma11(1000)
 m_1[i]=getT(x,0.05)
 m_2[i]=getT(x)
 m_3[i]=getT(x,0.15)
}
q_1 = quantile(m_1,.95)
q_2 = quantile(m_2,.95)
q_3 = quantile(m_3,.95)


#find quantile
n1 = 200
out3 = rep(NA,n1)
out4 = rep(NA,n1)
out5 = rep(NA,n1)
out6 = rep(NA,n1)
for(i in 1:n1){
  set.seed(i)
  y = arma11(1000)  
  #y = rnorm(500)
  z = y+rep(c(1,1.5),each = 250)
  
  out4[i] = getT(z,.05)
  out5[i] = getT(z,.1)
  out3[i] = getT(z,.15)
  out6[i] = getT(z,0)
  if(i%%50 == 0) {print(i)}
}
q0=quantile(out6,.95)
q1=quantile(out4,.95)
q2=quantile(out5,.95)
q3=quantile(out3,.95)
#
n1 = 200
out3 = rep(NA,n1)
out4 = rep(NA,n1)
out5 = rep(NA,n1)
out6 = rep(NA,n1)
for(i in 1:n1){
    set.seed(i)
    y = arma11(1000)  
    #y = rnorm(500)
    z = y+rep(c(0,1),each = 250)
    
    out4[i] = getT(z,.05)>q1
    out5[i] = getT(z,.1)>q2
    out3[i] = getT(z,.15)>q3
    #out6[i] = getT(z,0)>q4
    if(i%%50 == 0) {print(i)}
}
(p_0.5 = mean(out4))
(p_0.1 = mean(out5))
(p_0.15 = mean(out3))
matplot(power,type="b",ylim=c(0,1),xlab=expression(delta),lwd=2,lty=1)
abline(h=c(0,0.05,1),col=c("black","red","black"),lty=2)


###-----------------------------------------------------------------------
### Sample Experiments
###-----------------------------------------------------------------------
delta = (1:4)/2
newdelta=cbind(delta,(5:8)/2)
delta2=cbind(delta,2*delta)
n = 600
nRep = 999
nCase = length(delta)
test=rep(0,99)
out = array(NA,dim=c(nRep,nCase,3))
for(iRep in 1:nRep){
  X = arma11(n)
  test[iRep]=SN(X)
	for(iCase in 1:nCase){
  Y = X + rep(c(0,delta[iCase]),each=n/3)
  
  out[iRep,iCase,1] = getT(Y,0.05)>q1
  out[iRep,iCase,2] = getT(Y,0.1)>q2
  out[iRep,iCase,3] = getT(Y,0.15)>q3
  }
}










