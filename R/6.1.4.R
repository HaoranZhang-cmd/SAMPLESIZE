#' Testing the Hypothesis that the ROC Area is Equal to a Particular Value
#' @param null area under the null hypothesis
#' @param alter area under the alternative hypothesis
#' @param alpha significance level
#' @param beta 1-beta is power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing
#' @export
hypoarea<-function(null,alter,alpha,beta,R){
  a_null<-1.414*qnorm(null)
  a_alter<-1.414*qnorm(alter)
  V_null<-0.0099*exp(-a_null^2/2)*(5*a_null^2+8+(a_null^2+8)/R)
  print(V_null)
  V_alter<-0.0099*exp(-a_alter^2/2)*(5*a_alter^2+8+(a_alter^2+8)/R)
  n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)*sqrt(V_null)+qnorm(beta,lower.tail=FALSE)*sqrt(V_alter))^2/(null-alter)^2)
  N_one<-ceiling(n_one*(1+R))
  n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V_null)+qnorm(beta,lower.tail=FALSE)*sqrt(V_alter))^2/(null-alter)^2)
  N_two<-ceiling(n_two*(1+R))
  result<-list(n_one=NULL,N_one=NULL,n_two=NULL,N_two=NULL)
  result$n_one<-n_one
  result$N_one<-N_one
  result$n_two<-n_two
  result$N_two<-N_two
  return(result)
}

