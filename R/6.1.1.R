#' Sample Size Calculations for Estimating Sensitivity and Specificity
#' @param Se the conjectured sensitivity of the test
#' @param Sp the conjectured specificity of the test
#' @param alpha 1-alpha is the size of the desired confidence interval
#' @param beta  1-beta is the desired power
#' @param L the desired width of one-half of the CI
#' @param p the prevalence of the condition in the population
#' @return n.normal is a more precise estimation of the number of patients with condition
#'         n.rough is a rough estimation of the number of patients with condition
#'         N.normal is a more precise estimation of the total number of patients
#'         N.rough is a rough estimation of the total number of patients
#' @export
sampleSeSp<-function(Se,Sp,alpha,beta,L,p){
  V.Se<-Se*(1-Se)
  V.Sp<-Sp*(1-Sp)

  n.Se.normal<-ceiling(((qnorm(1-alpha/2)+qnorm(1-beta))*sqrt(V.Se))^2/L^2)
  n.Sp.normal<-ceiling(((qnorm(1-alpha/2)+qnorm(1-beta))*sqrt(V.Sp))^2/L^2)
  n.Se.rough<-ceiling((qnorm(1-alpha/2)*sqrt(V.Se))^2/L^2)
  n.Sp.rough<-ceiling((qnorm(1-alpha/2)*sqrt(V.Sp))^2/L^2)

  a.Se.normal<-p^2
  b.Se.normal<-(-2)*n.Se.normal*p-(qnorm(1-beta))^2*p*(1-p)
  c.Se.normal<-n.Se.normal^2
  N.Se.normal<-ceiling((-b.Se.normal+sqrt(b.Se.normal^2-4*a.Se.normal*c.Se.normal))/(2*a.Se.normal))
  a.Sp.normal<-(1-p)^2
  b.Sp.normal<-(-2)*n.Sp.normal*(1-p)-(qnorm(1-beta))^2*p*(1-p)
  c.Sp.normal<-n.Sp.normal^2
  N.Sp.normal<-ceiling((-b.Sp.normal+sqrt(b.Sp.normal^2-4*a.Sp.normal*c.Sp.normal))/(2*a.Sp.normal))

  a.Se.rough<-p^2
  b.Se.rough<-(-2)*n.Se.rough*p-(qnorm(1-beta))^2*p*(1-p)
  c.Se.rough<-n.Se.rough^2
  N.Se.rough<-ceiling((-b.Se.rough+sqrt(b.Se.rough^2-4*a.Se.rough*c.Se.rough))/(2*a.Se.rough))
  a.Sp.rough<-(1-p)^2
  b.Sp.rough<-(-2)*n.Sp.rough*(1-p)-(qnorm(1-beta))^2*p*(1-p)
  c.Sp.rough<-n.Sp.rough^2
  N.Sp.rough<-ceiling((-b.Sp.rough+sqrt(b.Sp.rough^2-4*a.Sp.rough*c.Sp.rough))/(2*a.Sp.rough))

  output<-list(n.Se.normal=NULL,n.Sp.normal=NULL,n.Se.rough=NULL,n.Sp.rough=NULL,N.Se.normal=NULL,N.Sp.normal=NULL,N.Se.rough=NULL,N.Sp.rough=NULL)
  output$n.Se.normal<-n.Se.normal
  output$n.Sp.normal<-n.Sp.normal
  output$n.Se.rough<-n.Se.rough
  output$n.Sp.rough<-n.Sp.rough
  output$N.Se.normal<-N.Se.normal
  output$N.Sp.normal<-N.Sp.normal
  output$N.Se.rough<-N.Se.rough
  output$N.Sp.rough<-N.Sp.rough
  return(output)
}


