# SAMPLE SIZE FOR ASSESSING NON-INFERIORITY OF TWO TESTS
#' Sample size for testing the non-inferiority when the accuracy is sensitivity
#' @param Se1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Se2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=1|T2=1)
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_Se<-function(Se1,Se2,p,Delta,alpha,beta,R){
  V0<-Se1+Se2-2*Se2*p
  VA<-V0-(Se1-Se2)^2
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Se1-Se2-Delta)^2)
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}

#' Sample size for testing the non-inferiority when the accuracy is specificity
#' @param Sp1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Sp2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=1|T2=1)
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_Sp<-function(Sp1,Sp2,p,Delta,alpha,beta,R){
  V0<-Sp1+Sp2-2*Sp2*p
  VA<-V0-(Sp1-Sp2)^2
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Sp1-Sp2-Delta)^2)
  N<-ceiling(n*(1+1/R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}

#' Sample size for testing the non-inferiority when the accuracy is AUC
#' under the binormal assumption
#' @param a1 the binormal parameter a for test 1 under the alternative hypothesis
#' @param b1 the binormal parameter b for test 1 under the alternative hypothesis
#' @param a2 the binormal parameter a for test 2 under the alternative hypothesis
#' @param b2 the binormal parameter b for test 2 under the alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_AUC_binormal<-function(a1,b1,a2,b2,rD,rN,Delta,alpha,beta,R){
  VA1<-0.0099*exp(-a1^2/2)*(5*a1^2+8+(a1^2+8)/R)
  VA2<-0.0099*exp(-a2^2/2)*(5*a2^2+8+(a2^2+8)/R)
  CA<-exp(-(a1^2+a2^2)/4)/12.5664*(rD+rN/R+rD^2*a1*a2/2)+exp(-(a1^2+a2^2)/4)/50.2655*(a1*a2*(rN^2+R*rD^2))/(2*R)-exp(-(a1^2+a2^2)/4)/25.1327*rD^2*a1*a2
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*(VA1+VA2-2*CA)/(pnorm(a1/sqrt(1+b1^2))-pnorm(a2/sqrt(1+b2^2))-Delta)^2)
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing the non-inferiority when the accuracy is AUC
#' under any distribution
#' @param theta1 AUC of test 1
#' @param theta2 AUC of test 2
#' @param r the correlation between the tests because of the paired design
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_AUC_any<-function(theta1,theta2,r,Delta,alpha,beta,R){
  VA<-theta1*(1-theta1)+theta2*(1-theta2)-2*r*sqrt(theta1*(1-theta1)*theta2*(1-theta2))
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(theta1-theta2-Delta)^2)
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing the non-inferiority when the accuracy is sensitivity at fixed FPR
#' We compare the z_transformed sensitivity
#' @param a1 binormal parameter for test 1 under alternative hypothesis
#' @param b1 binormal parameter for test 1 under alternative hypothesis
#' @param a2 binormal parameter for test 2 under alternative hypothesis
#' @param b2 binormal parameter for test 2 under alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param e the fixed FPR rate
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_fixedFPR<-function(a1,b1,a2,b2,rD,rN,e,Delta,alpha,beta,R){
  VA1<-1+b1^2/R+a1^2/2+(qnorm(e))^2*(b1^2*(1+R))/(2*R)
  VA2<-1+b2^2/R+a2^2/2+(qnorm(e))^2*(b2^2*(1+R))/(2*R)
  CA<-rD+(rN*b1*b2)/R+(rD^2*a1*a2)/2+(qnorm(e))^2*b1*b2*(rN^2+R*rD^2)/(2*R)+qnorm(e)*rD^2*(a1*b2+a2*b1)/2
  VA<-VA1+VA2-2*CA
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(a1+b1*qnorm(e)-a2-b2*qnorm(e)-Delta)^2)
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}



#' Sample size for testing the non-inferiority when the accuracy is partial AUC
#' @param a1 binormal parameter for test 1 under alternative hypothesis
#' @param b1 binormal parameter for test 1 under alternative hypothesis
#' @param a2 binormal parameter for test 2 under alternative hypothesis
#' @param b2 binormal parameter for test 2 under alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param e1 lower bound for partial area
#' @param e2 upper bound for partial area
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
non_inferiority_partialAUC<-function(a1,b1,a2,b2,rD,rN,e1,e2,Delta,alpha,beta,R){
  a<-c(a1,a2)
  b<-c(b1,b2)
  r1<-c(0,0)
  r2<-c(0,0)
  r3<-c(0,0)
  r4<-c(0,0)
  f<-c(0,0)
  g<-c(0,0)
  V<-c(0,0)
  for(i in 1:2){
    r1[i]<-exp(-a[i]^2/(2*(1+b[i]^2)))
    r2[i]<-1+b[i]^2
    r3[i]<-pnorm((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))-pnorm((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))
    r4[i]<-exp(-((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)-exp(-((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)
    f[i]<-r1[i]/sqrt(2*pi*r2[i])*r3[i]
    g[i]<-r1[i]*r4[i]/(2*pi*r2[i])-a[i]*b[i]*r1[i]*r3[i]/sqrt(2*pi*r2[i]^3)
    V[i]<-f[i]^2*(1+b[i]^2/R+a[i]^2/2)+g[i]^2*(b[i]^2*(1+R)/(2*R))
  }
  CA<-f[1]*f[2]*(rD+rN*b[1]*b[2]/R+rD^2*a[1]*a[2]/2)+g[1]*g[2]*b[1]*b[2]*(rN^2+R*rD^2)/(2*R)+f[1]*g[2]*rD^2*a[1]*b[2]/2+f[2]*g[1]*rD^2*a[2]*b[1]/2
  VA<-V[1]+V[2]-2*CA
  tmp1<-function(x){
    return(pnorm(a[1]+b[1]*qnorm(x)))
  }
  tmp2<-function(x){
    return(pnorm(a[2]+b[2]*qnorm(x)))
  }
  n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(integrate(tmp1,e1,e2)$value-integrate(tmp2,e1,e2)$value-Delta)^2)
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


# SAMPLE SIZE FOR ASSESSING EQUIVALENCY OF TWO TESTS
#' Sample size for testing equivalency when the accuracy is sensitivity
#' @param Se1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Se2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=1|T2=1)
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_Se<-function(Se1,Se2,p,Delta,alpha,beta,R){
  V0<-Se1+Se2-2*Se2*p
  VA<-V0-(Se1-Se2)^2
  if(Se1-Se2>0){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Se1-Se2-Delta)^2)
  }
  else if(Se1-Se2<0){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Se1-Se2+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*VA/Delta^2)
  }
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}

#' Sample size for testing equivalency when the accuracy is specificity
#' @param Sp1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Sp2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=1|T2=1)
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_Sp<-function(Sp1,Sp2,p,Delta,alpha,beta,R){
  V0<-Sp1+Sp2-2*Sp2*p
  VA<-V0-(Sp1-Sp2)^2
  if(Sp1-Sp2>0){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Sp1-Sp2-Delta)^2)
  }
  else if(Sp1-Sp2<0){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(Sp1-Sp2+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*VA/Delta^2)
  }
  N<-ceiling(n*(1+1/R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing equivalency when the accuracy is AUC
#' under the binormal assumption
#' @param a1 the binormal parameter a for test 1 under the alternative hypothesis
#' @param b1 the binormal parameter b for test 1 under the alternative hypothesis
#' @param a2 the binormal parameter a for test 2 under the alternative hypothesis
#' @param b2 the binormal parameter b for test 2 under the alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_AUC_binormal<-function(a1,b1,a2,b2,rD,rN,Delta,alpha,beta,R){
  VA1<-0.0099*exp(-a1^2/2)*(5*a1^2+8+(a1^2+8)/R)
  VA2<-0.0099*exp(-a2^2/2)*(5*a2^2+8+(a2^2+8)/R)
  CA<-exp(-(a1^2+a2^2)/4)/12.5664*(rD+rN/R+rD^2*a1*a2/2)+exp(-(a1^2+a2^2)/4)/50.2655*(a1*a2*(rN^2+R*rD^2))/(2*R)-exp(-(a1^2+a2^2)/4)/25.1327*rD^2*a1*a2
  if(a1/sqrt(1+b1^2)>a2/sqrt(1+b2^2)){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*(VA1+VA2-2*CA)/(pnorm(a1/sqrt(1+b1^2))-pnorm(a2/sqrt(1+b2^2))-Delta)^2)
  }
  else if(a1/sqrt(1+b1^2)<a2/sqrt(1+b2^2)){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*(VA1+VA2-2*CA)/(pnorm(a1/sqrt(1+b1^2))-pnorm(a2/sqrt(1+b2^2))+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*(VA1+VA2-2*CA)/Delta^2)
  }
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}

#' Sample size for testing equivalency when the accuracy is AUC
#' under any distribution
#' @param theta1 AUC of test 1
#' @param theta2 AUC of test 2
#' @param r the correlation between the tests because of the paired design
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_AUC_any<-function(theta1,theta2,r,Delta,alpha,beta,R){
  VA<-theta1*(1-theta1)+theta2*(1-theta2)-2*r*sqrt(theta1*(1-theta1)*theta2*(1-theta2))
  if(theta1>theta2){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(theta1-theta2-Delta)^2)
  }
  else if(theta1<theta2){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(theta1-theta2+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*VA/Delta^2)
  }
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing equivalency when the accuracy is sensitivity at fixed FPR
#' We compare the z_transformed sensitivity
#' @param a1 binormal parameter for test 1 under alternative hypothesis
#' @param b1 binormal parameter for test 1 under alternative hypothesis
#' @param a2 binormal parameter for test 2 under alternative hypothesis
#' @param b2 binormal parameter for test 2 under alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param e the fixed FPR rate
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_fixedFPR<-function(a1,b1,a2,b2,rD,rN,e,Delta,alpha,beta,R){
  VA1<-1+b1^2/R+a1^2/2+(qnorm(e))^2*(b1^2*(1+R))/(2*R)
  VA2<-1+b2^2/R+a2^2/2+(qnorm(e))^2*(b2^2*(1+R))/(2*R)
  CA<-rD+(rN*b1*b2)/R+(rD^2*a1*a2)/2+(qnorm(e))^2*b1*b2*(rN^2+R*rD^2)/(2*R)+qnorm(e)*rD^2*(a1*b2+a2*b1)/2
  VA<-VA1+VA2-2*CA
  if(a1+b1*qnorm(e)>a2+b2*qnorm(e)){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(a1+b1*qnorm(e)-a2-b2*qnorm(e)-Delta)^2)
  }
  else if(a1+b1*qnorm(e)<a2+b2*qnorm(e)){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(a1+b1*qnorm(e)-a2-b2*qnorm(e)+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*VA/Delta^2)
  }
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing equivalency when the accuracy is partial AUC
#' @param a1 binormal parameter for test 1 under alternative hypothesis
#' @param b1 binormal parameter for test 1 under alternative hypothesis
#' @param a2 binormal parameter for test 2 under alternative hypothesis
#' @param b2 binormal parameter for test 2 under alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param e1 lower bound for partial area
#' @param e2 upper bound for partial area
#' @param Delta the smallest difference in accuracy which is not acceptable
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size(paired design)
#' @export
equivalency_partialAUC<-function(a1,b1,a2,b2,rD,rN,e1,e2,Delta,alpha,beta,R){
  a<-c(a1,a2)
  b<-c(b1,b2)
  r1<-c(0,0)
  r2<-c(0,0)
  r3<-c(0,0)
  r4<-c(0,0)
  f<-c(0,0)
  g<-c(0,0)
  V<-c(0,0)
  for(i in 1:2){
    r1[i]<-exp(-a[i]^2/(2*(1+b[i]^2)))
    r2[i]<-1+b[i]^2
    r3[i]<-pnorm((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))-pnorm((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))
    r4[i]<-exp(-((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)-exp(-((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)
    f[i]<-r1[i]/sqrt(2*pi*r2[i])*r3[i]
    g[i]<-r1[i]*r4[i]/(2*pi*r2[i])-a[i]*b[i]*r1[i]*r3[i]/sqrt(2*pi*r2[i]^3)
    V[i]<-f[i]^2*(1+b[i]^2/R+a[i]^2/2)+g[i]^2*(b[i]^2*(1+R)/(2*R))
  }
  CA<-f[1]*f[2]*(rD+rN*b[1]*b[2]/R+rD^2*a[1]*a[2]/2)+g[1]*g[2]*b[1]*b[2]*(rN^2+R*rD^2)/(2*R)+f[1]*g[2]*rD^2*a[1]*b[2]/2+f[2]*g[1]*rD^2*a[2]*b[1]/2
  VA<-V[1]+V[2]-2*CA
  tmp1<-function(x){
    return(pnorm(a[1]+b[1]*qnorm(x)))
  }
  tmp2<-function(x){
    return(pnorm(a[2]+b[2]*qnorm(x)))
  }
  if(integrate(tmp1,e1,e2)$value>integrate(tmp2,e1,e2)$value){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(integrate(tmp1,e1,e2)$value-integrate(tmp2,e1,e2)$value-Delta)^2)
  }
  else if(integrate(tmp1,e1,e2)$value<integrate(tmp2,e1,e2)$value){
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*VA/(integrate(tmp1,e1,e2)$value-integrate(tmp2,e1,e2)$value+Delta)^2)
  }
  else{
    n<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta/2,lower.tail=FALSE))^2*VA/Delta^2)
  }
  N<-ceiling(n*(1+R))
  result<-list(n=NULL,N=NULL)
  result$n<-n
  result$N<-N
  return(result)
}


#' Sample size for testing the non-inferiority when the accuracies are sensitivity and specificity
#' calculate the sample size needed for testing the rTPR parameter by constructing joint confidence intervals for rTPR and rFPR
#' proposed by Alonzo and Pepe(2002),for paired design
#' @param TPR2 the sensitivity for test 2
#' @param TPPR the probability that both tests 1 and 2 are positive for patients with the condition
#' @param gamma a specific value of interest for rTPR under the alternative hypothesis
#' @param delta1 The null hypothesis is rTPR<delta1
#' @param alpha 1-alpha is the level of the joint confidence interval for rTPR and rFPR
#' @param beta 1-beta is the power of the test
#' @return n_pro total sample size for prospective study
#'         n_retro total sample size for retrospective study
#' @export
non_inferiority_rTPR<-function(TPR2,TPPR,gamma,delta1,alpha,beta){
  n_pro<-ceiling((qnorm(beta,lower.tail=FALSE)+qnorm(sqrt(1-alpha)))^2/log(gamma/delta1)^2*((gamma+1)*TPR2-2*TPPR)/(gamma*TPR2))
  n_retro<-ceiling((qnorm(sqrt(1-beta))+qnorm(sqrt(1-alpha)))^2/log(gamma/delta1)^2*((gamma+1)*TPR2-2*TPPR)/(gamma*TPR2))
  result<-list(n_pro=NULL,n_retro=NULL)
  result$n_pro<-n_pro
  result$n_retro<-n_retro
  return(result)
}


#' Sample size for testing the non-inferiority when the accuracies are sensitivity and specificity
#' calculate the sample size needed for testing the rFPR parameter by constructing joint confidence intervals for rTPR and rFPR
#' proposed by Alonzo and Pepe(2002),for paired design
#' @param FPR2 the false positive rate for test 2
#' @param FPPR the probability that both tests 1 and 2 are positive for patients without the condition
#' @param gamma a specific value of interest for rFPR under the alternative hypothesis
#' @param delta2 The null hypothesis is rFPR>delta2
#' @param alpha 1-alpha is the level of the joint confidence interval for rTPR and rFPR
#' @param beta 1-beta is the power of the test
#' @return n_pro total sample size for prospective study
#'         n_retro total sample size for retrospective study
#' @export
non_inferiority_rFPR<-function(FPR2,FPPR,gamma,delta2,alpha,beta){
  n_pro<-ceiling((qnorm(beta,lower.tail=FALSE)+qnorm(sqrt(1-alpha)))^2/log(gamma/delta2)^2*((gamma+1)*FPR2-2*FPPR)/(gamma*FPR2))
  n_retro<-ceiling((qnorm(sqrt(1-beta))+qnorm(sqrt(1-alpha)))^2/log(gamma/delta2)^2*((gamma+1)*FPR2-2*FPPR)/(gamma*FPR2))
  result<-list(n_pro=NULL,n_retro=NULL)
  result$n_pro<-n_pro
  result$n_retro<-n_retro
  return(result)
}


#' Sample size for testing non-inferiority when the accuracy is PPV
#' using the method proposed by Moskowitz and Pepe(2006)
#' @param PPV2 the conjectured value of PPV of test 2
#' @param p a vector with eight dimensions
#'        For patients without the condition：
#'        p[1] the proportion of patients testing positive on both tests
#'        p[2] the proportion testing positive on test 1 and negative on test 2
#'        p[3] the proportion testing positive on test 2 and negative on test 1
#'        p[4] the proportion testing negative on both tests.
#'        For patients with the condition:
#'        p[5] the proportion of patients testing positive on both tests
#'        p[6] the proportion testing positive on test 1 and negative on test 2
#'        p[7] the proportion testing positive on test 2 and negative on test 1
#'        p[8] the proportion testing negative on both tests
#' @param gamma a specific value of interest for rPPV under the alternative hypothesis
#' @param delta The null hypothesis is rPPV≤delta.
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @return n_one total sample size in a one-sided hypothesis testing
#'         n_two total sample sizein a two-sided hypothesis testing
#' @export
non_inferiority_PPV<-function(PPV2,p,gamma,delta,alpha,beta){
  tmp1_one<-((qnorm(beta,lower.tail=FALSE)+qnorm(alpha,lower.tail=FALSE))/log(gamma/delta))^2
  tmp2<-1/((p[5]+p[6])*(p[5]+p[7]))
  tmp3<-2*(p[7]+p[3])*gamma*PPV2^2+PPV2*(p[5]*(1-gamma)-p[6])+p[6]+p[7]*(1-3*gamma*PPV2)
  tmp1_two<-((qnorm(beta,lower.tail=FALSE)+qnorm(alpha/2,lower.tail=FALSE))/log(gamma/delta))^2
  n_one<-ceiling(tmp1_one*tmp2*tmp3)
  n_two<-ceiling(tmp1_two*tmp2*tmp3)
  result<-list(n_one=NULL,n_two=NULL)
  result$n_one<-n_one
  result$n_two<-n_two
  return(result)
}

#' Sample size for testing non-inferiority when the accuracy is NPV
#' using the method proposed by Moskowitz and Pepe(2006)
#' @param NPV2 the conjectured value of NPV of test 2
#' @param p a vector with eight dimensions
#'        For patients without the condition：
#'        p[1] the proportion of patients testing positive on both tests
#'        p[2] the proportion testing positive on test 1 and negative on test 2
#'        p[3] the proportion testing positive on test 2 and negative on test 1
#'        p[4] the proportion testing negative on both tests.
#'        For patients with the condition：
#'        p[5] the proportion of patients testing positive on both tests
#'        p[6] the proportion testing positive on test 1 and negative on test 2
#'        p[7] the proportion testing positive on test 2 and negative on test 1
#'        p[8] the proportion testing negative on both tests
#' @param gamma a specific value of interest for rNPV under the alternative hypothesis
#' @param delta The null hypothesis is rNPV≤delta.
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @return n_one total sample size in a one-sided hypothesis testing
#'         n_two total sample size in a two-sided hypothesis testing
#' @export
non_inferiority_NPV<-function(NPV2,p,gamma,delta,alpha,beta){
  tmp1_one<-((qnorm(beta,lower.tail=FALSE)+qnorm(alpha,lower.tail=FALSE))/log(gamma/delta))^2
  tmp2<-1/((p[2]+p[4])*(p[3]+p[4]))
  tmp3<--2*(p[4]+p[8])*gamma*NPV2^2+NPV2*(gamma*(p[4]-p[2])+p[4]-p[3])+p[2]+p[3]
  tmp1_two<-((qnorm(beta,lower.tail=FALSE)+qnorm(alpha/2,lower.tail=FALSE))/log(gamma/delta))^2
  n_one<-ceiling(tmp1_one*tmp2*tmp3)
  n_two<-ceiling(tmp1_two*tmp2*tmp3)
  result<-list(n_one=NULL,n_two=NULL)
  result$n_one<-n_one
  result$n_two<-n_two
  return(result)
}

