#' Sample Size for Estimating Sensitivity at Fixed FPR
#' @param a0 the parameters of the assumed underlying binormal distribution
#' @param b0 the parameters of the assumed underlying binormal distribution
#' @param e the fixed FPR rate
#' @param L half of the length of CI for the sensitivity
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing
#' @export
fixedFPRsample<-function(a0,b0,e,L,alpha,beta,R){
  Se0<-pnorm(a0+b0*qnorm(e))
  a_lower<-qnorm(Se0-L)-b0*qnorm(e)
  a_upper<-qnorm(Se0+L)-b0*qnorm(e)
  print(c(a_lower,a_upper))
  z_lower<-a_lower+b0*qnorm(e)
  z_upper<-a_upper+b0*qnorm(e)
  length<-(z_upper-z_lower)/2
  print(length)
  V<-1+b0^2/R+a0^2/2+qnorm(e)^2*b0^2*(1+R)/(2*R)
  print(V)
  n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*V/length^2)
  n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*V/length^2)
  N_one<-ceiling(n_one*(1+R))
  N_two<-ceiling(n_two*(1+R))
  result<-list(n_one=NULL,N_one=NULL,n_two=NULL,N_two=NULL)
  result$n_one<-n_one
  result$N_one<-N_one
  result$n_two<-n_two
  result$N_two<-N_two
  return(result)
}

