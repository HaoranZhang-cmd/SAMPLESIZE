#' Sample Size for Estimating the Area Under the ROC Curve
#' @param A the conjectured AUC of the test
#' @param R the ratio of the number of patients without the condition to patients with the condition in the study sample
#' @param alpha 1-alpha is the size of the desired confidence interval
#' @param beta 1-beta is the desired power
#' @param L the desired width of one-half of the CI
#' @return n.normal is a more precise estimation of the number of patients with condition
#'         n.rough is a rough estimation of the number of patients with condition
#'         N.normal is a more precise estimation of the total number of patients
#'         N.rough is a rough estimation of the total number of patients
#' @export
samplearea<-function(A,R,alpha,beta,L){
  V.A<-c(0,0,0,0)
  n.A.normal<-c(0,0,0,0)
  n.A.rough<-c(0,0,0,0)
  #for specific distributions such as exponential distribution
  V.A[1]<-A/((2-A)*R)+2*A^2/(1+A)-A^2*(1/R+1)
  a<-1.414*qnorm(A)
  V.A[2]<-(0.0099)*exp(-a^2/2)*((5*a^2+8)+(a^2+8)/R)
  #for any distribution
  V.A[3]<-A*(1-A)
  #for binormal distribution
  V.A[4]<-(0.0099)*exp(-a^2/2)*((5*a^2+8)+(a^2+8)/R)-0.0398*a^2*exp(-a^2/2)
  n.A.normal<-ceiling(((qnorm(1-alpha/2)+qnorm(1-beta))*sqrt(V.A))^2/L^2)
  n.A.rough<-ceiling((qnorm(1-alpha/2)*sqrt(V.A))^2/L^2)
  N.A.normal<-ceiling(n.A.normal*(1+R))
  N.A.rough<-ceiling(n.A.rough*(1+R))
  output<-list(n.A.normal=NULL,n.A.rough=NULL,N.A.normal=NULL,N.A.rough=NULL)
  output$n.A.normal<-n.A.normal
  output$n.A.rough<-n.A.rough
  output$N.A.normal<-N.A.normal
  output$N.A.rough<-N.A.rough
  return(output)
}

