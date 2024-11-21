# Fixed- Reader MRMC Design
#' this function calculate the sample size needed for the fixed-reader MRMC design
#' using the method proposed by Blume(2009)
#' @param theta AUC of both tests under the null hypothesis
#' @param theta1 AUC of test 1 under the alternative hypothesis
#' @param theta2 AUC of test 2 under the alternative hypothesis
#' @param r the correlation between the tests because of the paired design
#' @param J number of readers
#' @param rho the correlation between different readers on the difference in ROC areas
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size needed
#' @export
sample_fixedreader_AUC_Blume<-function(theta,theta1,theta2,r,J,rho,alpha,beta,R){
    V0<-(2*theta*(1-theta)-2*r*theta*(1-theta))*(1/J+(J-1)*rho/J)
    VA<-(theta1*(1-theta1)+theta2*(1-theta2)-2*r*sqrt(theta1*(1-theta1)*theta2*(1-theta2)))*(1/J+(J-1)*rho/J)
    print(c(V0,VA))
    n<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/(theta1-theta2)^2)
    N<-ceiling(n*(1+R))
    result<-list(n=NULL,N=NULL)
    result$n<-n
    result$N<-N
    return(result)
}

#' this function calculate the sample size needed for the fixed-reader MRMC design
#' using the method proposed by Obuchowski(1995)
#' @param a1_null the binormal parameter a for test 1 under the null hypothesis
#' @param a2_null the binormal parameter a for test 2 under the null hypothesis
#' @param a1_alter the binormal parameter a for test 1 under the alternative hypothesis
#' @param a2_alter the binormal parameter a for test 2 under the alternative hypothesis
#' @param rD the correlation of the underlying bivariate binormal distribution for patients with the condition
#' @param rN the correlation of the underlying bivariate binormal distribution for patients without the condition
#' @param J number of readers
#' @param rho The difference between estimated correlation when the same patients are evaluated by different readers using the same test
#'            and The estimated correlation when the same patients are evaluated by different readers using different tests
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param R the ratio of undiseased and diseased subjects
#' @return n sample size of diseased subjects
#'         N total sample size needed
#' @export
sample_fixedreader_AUC_Obuchowski<-function(a1_null,a2_null,a1_alter,a2_alter,rD,rN,J,rho,alpha,beta,R){
  V01<-0.0099*exp(-a1_null^2/2)*(5*a1_null^2+8+(a1_null^2+8)/R)
  V02<-0.0099*exp(-a2_null^2/2)*(5*a2_null^2+8+(a2_null^2+8)/R)
  VA1<-0.0099*exp(-a1_alter^2/2)*(5*a1_alter^2+8+(a1_alter^2+8)/R)
  VA2<-0.0099*exp(-a2_alter^2/2)*(5*a2_alter^2+8+(a2_alter^2+8)/R)
  C0<-exp(-(a1_null^2+a2_null^2)/4)/12.5664*(rD+rN/R+rD^2*a1_null*a2_null/2)+exp(-(a1_null^2+a2_null^2)/4)/50.2655*(a1_null*a2_null*(rN^2+R*rD^2))/(2*R)-exp(-(a1_null^2+a2_null^2)/4)/25.1327*rD^2*a1_null*a2_null
  CA<-exp(-(a1_alter^2+a2_alter^2)/4)/12.5664*(rD+rN/R+rD^2*a1_alter*a2_alter/2)+exp(-(a1_alter^2+a2_alter^2)/4)/50.2655*(a1_alter*a2_alter*(rN^2+R*rD^2))/(2*R)-exp(-(a1_alter^2+a2_alter^2)/4)/25.1327*rD^2*a1_alter*a2_alter
  Delta<-abs(pnorm(a1_alter/sqrt(2))-pnorm(a2_alter/sqrt(2)))
  V0<-(V01+V02-2*C0)*(1/J+(J-1)*rho/J)
  VA<-(VA1+VA2-2*CA)*(1/J+(J-1)*rho/J)
  result<-list(n=NULL,N=NULL)
  result$n<-ceiling((qnorm(alpha/2,lower.tail=FALSE)*sqrt(V0)+qnorm(beta,lower.tail=FALSE)*sqrt(VA))^2/Delta^2)
  result$N<-ceiling(result$n*(1+R))
  return(result)
}

