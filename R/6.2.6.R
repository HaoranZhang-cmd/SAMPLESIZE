#' Sample Size for Comparing Tests Sensitivity at Fixed FPR
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
#' @param e the fixed FPR rate
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing(paired design)
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing(paired design)
#' @export
samplefixedFPR<-function(a1_null,b1_null,a2_null,b2_null,a1_alter,b1_alter,a2_alter,b2_alter,rD,rN,e,alpha,beta,R){
  z1_alter<-a1_alter+b1_alter*qnorm(e)
  z2_alter<-a2_alter+b2_alter*qnorm(e)
  V01<-1+b1_null^2/R+a1_null^2/2+(qnorm(e))^2*(b1_null^2*(1+R))/(2*R)
  V02<-1+b2_null^2/R+a2_null^2/2+(qnorm(e))^2*(b2_null^2*(1+R))/(2*R)
  VA1<-1+b1_alter^2/R+a1_alter^2/2+(qnorm(e))^2*(b1_alter^2*(1+R))/(2*R)
  VA2<-1+b2_alter^2/R+a2_alter^2/2+(qnorm(e))^2*(b2_alter^2*(1+R))/(2*R)
  Delta<-abs(z1_alter-z2_alter)
  C0<-rD+(rN*b1_null*b2_null)/R+(rD^2*a1_null*a2_null)/2+(qnorm(e))^2*b1_null*b2_null*(rN^2+R*rD^2)/(2*R)+qnorm(e)*rD^2*(a1_null*b2_null+a2_null*b1_null)/2
  CA<-rD+(rN*b1_alter*b2_alter)/R+(rD^2*a1_alter*a2_alter)/2+(qnorm(e))^2*b1_alter*b2_alter*(rN^2+R*rD^2)/(2*R)+qnorm(e)*rD^2*(a1_alter*b2_alter+a2_alter*b1_alter)/2
  V0<-V01+V02-2*C0
  VA<-VA1+VA2-2*CA
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

