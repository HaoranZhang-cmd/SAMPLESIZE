#' Sample Size for Comparing Tests Sensitivity and/or Specificity
#' @param Se1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Se2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=1|T2=1)
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing
#' @export
samplecompareSe<-function(Se1,Se2,p,alpha,beta,R){
    V0<-Se1+Se2-2*Se2*p
    VA<-V0-(Se1-Se2)^2
    n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/(Se1-Se2)^2)
    n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/(Se1-Se2)^2)
    result<-list(n_one=NULL,n_two=NULL,N_one=NULL,N_two=NULL)
    result$n_one<-n_one
    result$n_two<-n_two
    result$N_one<-ceiling(n_one*(1+R))
    result$N_two<-ceiling(n_two*(1+R))
    return(result)
}


#' Sample Size for Comparing Tests Sensitivity and/or Specificity
#' @param Sp1 the conjectured value of sensitivity under the alternative hypothesis of test 1
#' @param Sp2 the conjectured value of sensitivity under the alternative hypothesis of test 2
#' @param p P(T1=0|T2=0)
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased subjects and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing
#' @export
samplecompareSp<-function(Sp1,Sp2,p,alpha,beta,R){
  V0<-Sp1+Sp2-2*Sp2*p
  VA<-V0-(Sp1-Sp2)^2
  n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/(Sp1-Sp2)^2)
  n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/(Sp1-Sp2)^2)
  result<-list(n_one=NULL,n_two=NULL,N_one=NULL,N_two=NULL)
  result$n_one<-n_one
  result$n_two<-n_two
  result$N_one<-ceiling(n_one*(1+1/R))
  result$N_two<-ceiling(n_two*(1+1/R))
  return(result)
}


