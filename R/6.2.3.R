#' Sample Size for Comparing Tests Positive Predictive Values
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
#' @param delta The null hypothesis is rPPV=delta.In most cases we let delta=1.
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @return n_one total sample size in a one-sided hypothesis testing
#'         n_two total sample size in a two-sided hypothesis testing
#' @export
samplecomparePPV<-function(PPV2,p,gamma,delta,alpha,beta){
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

#' Sample Size for Comparing Tests Negative Predictive Values
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
#' @param gamma a specific value of interest for rPPV under the alternative hypothesis
#' @param delta The null hypothesis is rPPV≤delta.In most cases we let delta=1.
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @return n_one total sample size in a one-sided hypothesis testing
#'         n_two total sample size in a two-sided hypothesis testing
#' @export
samplecompareNPV<-function(NPV2,p,gamma,delta,alpha,beta){
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

