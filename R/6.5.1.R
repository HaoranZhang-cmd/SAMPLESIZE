# SAMPLE SIZE DETERMINATION FOR MULTI-READER STUDIES
#' calculate the power of the test with J readers and n individuals with the condition when pilot data is not available
#' the accuracy is sensitivity
#' @param Se the average sensitivity of readers under the null hypothesis
#' @param delta the difference of the mean of sensitivity under the alternative hypothesis
#' @param var_b The estimated variability from different readers interpreting the results of the same patients using the same test
#' @param var_w The estimated variability from the same reader interpreting the results of the same patients using the same test on different occasions
#' @param rho_1 The estimated correlation when the same patients are evaluated by the same reader using different tests
#' @param rho_2 The estimated correlation when the same patients are evaluated by different readers using the same test
#' @param rho_3 The estimated correlation when the same patients are evaluated by different readers using different tests
#' @param rho_b The estimated correlation when the same readers evaluate patients using different tests
#' @param J number of readers for each diagnostic test
#' @param n number of individuals with the condition
#'          For the paired-patient study designs,n is the total number of individuals with the condition needed.
#'          For unpaired-patient designs,n is the number of individuals with the condition needed per diagnostic test.
#'          For the paired-patient per reader design,n individuals with the condition are needed for each of the J readers.
#' @param Q number of times each reader interprets the test results of each patient using the same test (often Q=1)
#' @param alpha significance level
#' @return the power of the test with J readers and n individuals with the condition when pilot data is not available
#' @export
MRMC_nopilot_Se<-function(Se,delta,var_b,var_w,rho_1,rho_2,rho_3,rho_b,J,n,Q,alpha){
  var_c<-Se*(1-Se)/n
  lambda<-J*delta^2/(2*(var_b*(1-rho_b)+var_w/Q+var_c*(1-rho_1+(J-1)*(rho_2-rho_3))))
  return((1-pf(qf(1-alpha,1,J-1),1,J-1,ncp=lambda)))
}


#' calculate the power of the test with J readers and n individuals with the condition when pilot data is not available
#' the accuracy is AUC
#' @param A the average AUC of readers under the null hypothesis
#' @param delta the difference of the mean of AUC under the alternative hypothesis
#' @param var_b The estimated variability from different readers interpreting the results of the same patients using the same test
#' @param var_w The estimated variability from the same reader interpreting the results of the same patients using the same test on different occasions
#' @param rho_1 The estimated correlation when the same patients are evaluated by the same reader using different tests
#' @param rho_2 The estimated correlation when the same patients are evaluated by different readers using the same test
#' @param rho_3 The estimated correlation when the same patients are evaluated by different readers using different tests
#' @param rho_b The estimated correlation when the same readers evaluate patients using different tests
#' @param J number of readers for each diagnostic test
#' @param n number of individuals with the condition
#'          For the paired-patient study designs,n is the total number of individuals with the condition needed.
#'          For unpaired-patient designs,n is the number of individuals with the condition needed per diagnostic test.
#'          For the paired-patient per reader design,n individuals with the condition are needed for each of the J readers.
#' @param Q number of times each reader interprets the test results of each patient using the same test (often Q=1)
#' @param R the ratio of number of individuals without the condition and the number of individuals with the condition
#' @param alpha significance level
#' @return the power of the test with J readers and n individuals with the condition when pilot data is not available
#' @export
MRMC_nopilot_AUC<-function(A,delta,var_b,var_w,rho_1,rho_2,rho_3,rho_b,J,n,Q,R,alpha){
  a<-qnorm(A)*1.414
  var_c<-0.0099*exp(-a^2/2)*(5*a^2+8+(a^2+8)/R)/n
  lambda<-J*delta^2/(2*(var_b*(1-rho_b)+var_w/Q+var_c*(1-rho_1+(J-1)*(rho_2-rho_3))))
  return((1-pf(qf(1-alpha,1,J-1),1,J-1,ncp=lambda)))
}



#' calculate the power of the test with J readers and n individuals with the condition when pilot data is not available
#' the accuracy is sensitivity at fixed FPR
#' @param a binormal parameter a under the null hypothesis
#' @param b binormal parameter b under the null hypothesis
#' @param e the fixed FPR rate
#' @param delta the difference of the mean of z-transformed sensitivity at fixed FPR under the alternative hypothesis
#' @param var_b The estimated variability from different readers interpreting the results of the same patients using the same test
#' @param var_w The estimated variability from the same reader interpreting the results of the same patients using the same test on different occasions
#' @param rho_1 The estimated correlation when the same patients are evaluated by the same reader using different tests
#' @param rho_2 The estimated correlation when the same patients are evaluated by different readers using the same test
#' @param rho_3 The estimated correlation when the same patients are evaluated by different readers using different tests
#' @param rho_b The estimated correlation when the same readers evaluate patients using different tests
#' @param J number of readers for each diagnostic test
#' @param n number of individuals with the condition
#'          For the paired-patient study designs,n is the total number of individuals with the condition needed.
#'          For unpaired-patient designs,n is the number of individuals with the condition needed per diagnostic test.
#'          For the paired-patient per reader design,n individuals with the condition are needed for each of the J readers.
#' @param Q number of times each reader interprets the test results of each patient using the same test (often Q=1)
#' @param R the ratio of number of individuals without the condition and the number of individuals with the condition
#' @param alpha significance level
#' @return the power of the test with J readers and n individuals with the condition when pilot data is not available
#' @export
MRMC_nopilot_fixedFPR<-function(a,b,e,delta,var_b,var_w,rho_1,rho_2,rho_3,rho_b,J,n,Q,R,alpha){
  var_c<-((1+b^2/R+a^2/2)+qnorm(e)^2*(b^2*(1+R)/(2*R)))/n
  lambda<-J*delta^2/(2*(var_b*(1-rho_b)+var_w/Q+var_c*(1-rho_1+(J-1)*(rho_2-rho_3))))
  return((1-pf(qf(1-alpha,1,J-1),1,J-1,ncp=lambda)))
}



#' calculate the power of the test with J readers and n individuals with the condition when pilot data is not available
#' the accuracy is partial AUC
#' @param a binormal parameter a under the null hypothesis
#' @param b binormal parameter b under the null hypothesis
#' @param e1 lower bound for FPR of partial area
#' @param e2 upper bound for FPR of partial area
#' @param delta the difference of the mean of partial area index under the alternative hypothesis
#' @param var_b The estimated variability from different readers interpreting the results of the same patients using the same test
#' @param var_w The estimated variability from the same reader interpreting the results of the same patients using the same test on different occasions
#' @param rho_1 The estimated correlation when the same patients are evaluated by the same reader using different tests
#' @param rho_2 The estimated correlation when the same patients are evaluated by different readers using the same test
#' @param rho_3 The estimated correlation when the same patients are evaluated by different readers using different tests
#' @param rho_b The estimated correlation when the same readers evaluate patients using different tests
#' @param J number of readers for each diagnostic test
#' @param n number of individuals with the condition
#'          For the paired-patient study designs,n is the total number of individuals with the condition needed.
#'          For unpaired-patient designs,n is the number of individuals with the condition needed per diagnostic test.
#'          For the paired-patient per reader design,n individuals with the condition are needed for each of the J readers.
#' @param Q number of times each reader interprets the test results of each patient using the same test (often Q=1)
#' @param R the ratio of number of individuals without the condition and the number of individuals with the condition
#' @param alpha significance level
#' @return the power of the test with J readers and n individuals with the condition when pilot data is not available
#' @export
MRMC_nopilot_partialAUC<-function(a,b,e1,e2,delta,var_b,var_w,rho_1,rho_2,rho_3,rho_b,J,n,Q,R,alpha){
  r1<-exp(-a^2/(2*(1+b^2)))
  r2<-1+b^2
  r3<-pnorm((qnorm(e2)+a*b/(1+b^2))*sqrt(1+b^2))-pnorm((qnorm(e1)+a*b/(1+b^2))*sqrt(1+b^2))
  r4<-exp(-((qnorm(e1)+a*b/(1+b^2))*sqrt(1+b^2))^2/2)-exp(-((qnorm(e2)+a*b/(1+b^2))*sqrt(1+b^2))^2/2)
  f<-r1/sqrt(2*pi*r2)*r3
  g<-r1*r4/(2*pi*r2)-a*b*r1*r3/sqrt(2*pi*r2^3)
  var_c<-(f^2*(1+b^2/R+a^2/2)+g^2*(b^2*(1+R)/(2*R)))/n
  lambda<-J*delta^2/(2*(var_b*(1-rho_b)+var_w/Q+var_c*(1-rho_1+(J-1)*(rho_2-rho_3))))
  return((1-pf(qf(1-alpha,1,J-1),1,J-1,ncp=lambda)))
}


