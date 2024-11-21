#' Sample Size for Comparing Tests Area Under the ROC Curve
#' under the binormal assumption
#' @param a1_null the binormal parameter a for test 1 under the null hypothesis
#' @param a2_null the binormal parameter a for test 2 under the null hypothesis
#' @param a1_alter the binormal parameter a for test 1 under the alternative hypothesis
#' @param a2_alter the binormal parameter a for test 2 under the alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing(paired design)
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing(paired design)
#' @export
samplecompareAUC_binormal<-function(a1_null,a2_null,a1_alter,a2_alter,rD,rN,alpha,beta,R){
  V01<-0.0099*exp(-a1_null^2/2)*(5*a1_null^2+8+(a1_null^2+8)/R)
  V02<-0.0099*exp(-a2_null^2/2)*(5*a2_null^2+8+(a2_null^2+8)/R)
  VA1<-0.0099*exp(-a1_alter^2/2)*(5*a1_alter^2+8+(a1_alter^2+8)/R)
  VA2<-0.0099*exp(-a2_alter^2/2)*(5*a2_alter^2+8+(a2_alter^2+8)/R)
  C0<-exp(-(a1_null^2+a2_null^2)/4)/12.5664*(rD+rN/R+rD^2*a1_null*a2_null/2)+exp(-(a1_null^2+a2_null^2)/4)/50.2655*(a1_null*a2_null*(rN^2+R*rD^2))/(2*R)-exp(-(a1_null^2+a2_null^2)/4)/25.1327*rD^2*a1_null*a2_null
  CA<-exp(-(a1_alter^2+a2_alter^2)/4)/12.5664*(rD+rN/R+rD^2*a1_alter*a2_alter/2)+exp(-(a1_alter^2+a2_alter^2)/4)/50.2655*(a1_alter*a2_alter*(rN^2+R*rD^2))/(2*R)-exp(-(a1_alter^2+a2_alter^2)/4)/25.1327*rD^2*a1_alter*a2_alter
  Delta<-abs(pnorm(a1_alter/sqrt(2))-pnorm(a2_alter/sqrt(2)))
  n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)*sqrt(V01+V02-2*C0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA1+VA2-2*CA))^2/Delta^2)
  n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V01+V02-2*C0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA1+VA2-2*CA))^2/Delta^2)
  N_one<-ceiling(n_one*(1+R))
  N_two<-ceiling(n_two*(1+R))
  result<-list(n_one=NULL,n_two=NULL,N_one=NULL,N_two=NULL)
  result$n_one<-n_one
  result$n_two<-n_two
  result$N_one<-N_one
  result$N_two<-N_two
  return(result)
}

#' Sample Size for Comparing Tests Area Under the ROC Curve
#' under any distribution
#' @param theta AUC of both tests under the null hypothesis
#' @param theta1 AUC of test 1 under the alternative hypothesis
#' @param theta2 AUC of test 2 under the alternative hypothesis
#' @param r the correlation between the tests because of the paired design
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing(paired design)
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing(paired design)
#' @export
samplecompareAUC_any<-function(theta,theta1,theta2,r,alpha,beta,R){
  V0<-2*theta*(1-theta)-2*r*theta*(1-theta)
  VA<-theta1*(1-theta1)+theta2*(1-theta2)-2*r*sqrt(theta1*(1-theta1)*theta2*(1-theta2))
  Delta<-abs(theta1-theta2)
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

