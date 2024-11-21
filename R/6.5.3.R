# MRMC Sample Size Calculations with Pilot Data
#' calculate the power of the test with J readers and N total patients when pilot data is available
#' @param delta the difference of accuracy of test 1 under the alternative hypothesis
#' @param J_pilot number of readers in the pilot study
#' @param N_pilot number of total patients in the pilot study
#' @param J number of readers in the current study
#' @param N number of total patients in the current study
#' @param MSTRP the modality-by- reader-by-patient interaction mean squares
#' @param MSTR the modality-by-reader interaction mean squares
#' @param MSTP the modality-by-patient interaction mean squares
#' @param alpha significance level
#' @return the power of the test with J readers and N total patients when pilot data is available
#' @export
MRMC_pilot<-function(delta,J_pilot,N_pilot,J,N,MSTRP,MSTR,MSTP,alpha){
  sigma<-max(MSTRP,0)
  sigma_TR<-max((MSTR-MSTRP)/N_pilot,0)
  sigma_TP<-max((MSTP-MSTRP)/J_pilot,0)
  df2<-J-1
  if(sigma_TP>0){
    tmp1<-(N*sigma_TR+J*sigma_TP+sigma)^2
    tmp2<-(N*sigma_TR+sigma)^2/(J-1)
    tmp3<-(J*sigma_TP+sigma)^2/(N-1)
    tmp4<-sigma^2/((J-1)*(N-1))
    df2<-tmp1/(tmp2+tmp3+tmp4)
  }
  lambda<-(delta)^2*J*N/(2*(N*sigma_TR+J*sigma_TP+sigma))
  return((1-pf(qf(1-alpha,1,df2),1,df2,ncp=lambda)))
}

