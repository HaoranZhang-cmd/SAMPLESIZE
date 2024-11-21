#' Sample Size for Comparing Tests Partial Area Under the ROC Curve
#' @param a1_null binormal parameter for test 1 under null hypothesis
#' @param b1_null binormal parameter for test 1 under null hypothesis
#' @param a2_null binormal parameter for test 2 under null hypothesis
#' @param b2_null binormal parameter for test 2 under null hypothesis
#' @param a1_alter binormal parameter for test 1 under alternative hypothesis
#' @param b1_alter binormal parameter for test 1 under alternative hypothesis
#' @param a2_alter binormal parameter for test 2 under alternative hypothesis
#' @param b2_alter binormal parameter for test 2 under alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param e1 lower bound for partial area
#' @param e2 upper bound for partial area
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing(paired design)
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing(paired design)
#' @export
samplepartialAUC<-function(a1_null,b1_null,a2_null,b2_null,a1_alter,b1_alter,a2_alter,b2_alter,rD,rN,e1,e2,alpha,beta,R){
  a<-c(a1_null,a2_null,a1_alter,a2_alter)
  b<-c(b1_null,b2_null,b1_alter,b2_alter)
  r1<-rep(0,length=4)
  r2<-rep(0,length=4)
  r3<-rep(0,length=4)
  r4<-rep(0,length=4)
  f<-rep(0,length=4)
  g<-rep(0,length=4)
  V<-rep(0,length=4)
  for(i in 1:4){
    r1[i]<-exp(-a[i]^2/(2*(1+b[i]^2)))
    r2[i]<-1+b[i]^2
    r3[i]<-pnorm((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))-pnorm((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))
    r4[i]<-exp(-((qnorm(e1)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)-exp(-((qnorm(e2)+a[i]*b[i]/(1+b[i]^2))*sqrt(1+b[i]^2))^2/2)
    f[i]<-r1[i]/sqrt(2*pi*r2[i])*r3[i]
    g[i]<-r1[i]*r4[i]/(2*pi*r2[i])-a[i]*b[i]*r1[i]*r3[i]/sqrt(2*pi*r2[i]^3)
    V[i]<-f[i]^2*(1+b[i]^2/R+a[i]^2/2)+g[i]^2*(b[i]^2*(1+R)/(2*R))
  }
  C0<-f[1]*f[2]*(rD+rN*b[1]*b[2]/R+rD^2*a[1]*a[2]/2)+g[1]*g[2]*b[1]*b[2]*(rN^2+R*rD^2)/(2*R)+f[1]*g[2]*rD^2*a[1]*b[2]/2+f[2]*g[1]*rD^2*a[2]*b[1]/2
  CA<-f[3]*f[4]*(rD+rN*b[3]*b[4]/R+rD^2*a[3]*a[4]/2)+g[3]*g[4]*b[3]*b[4]*(rN^2+R*rD^2)/(2*R)+f[3]*g[4]*rD^2*a[3]*b[4]/2+f[4]*g[3]*rD^2*a[4]*b[3]/2
  V0<-V[1]+V[2]-2*C0
  VA<-V[3]+V[4]-2*CA
  tmp1<-function(x){
    return(pnorm(a[3]+b[3]*qnorm(x)))
  }
  tmp2<-function(x){
    return(pnorm(a[4]+b[4]*qnorm(x)))
  }
  Delta<-abs(integrate(tmp1,e1,e2)$value-integrate(tmp2,e1,e2)$value)
  n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/Delta^2)
  n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/Delta^2)
  N_one<-ceiling(n_one*(1+R))
  N_two<-ceiling(n_two*(1+R))
  result<-list(n_one=NULL,n_two=NULL,N_one=NULL,N_two=NULL)
  result$n_one<-n_one
  result$n_two<-n_two
  result$N_one<-N_one
  result$N_two<-N_two
  return(result)
}

