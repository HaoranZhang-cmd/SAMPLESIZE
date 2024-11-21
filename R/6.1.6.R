#' Sample Size for Estimating the Partial Area Under the ROC Curve
#' @param a binormal parameters
#' @param b binormal parameters
#' @param e1 lower bound for partial area
#' @param e2 upper bound for partial area
#' @param alpha significance level
#' @param beta 1-beta is the power of the test
#' @param L the desired width of one-half of the CI
#' @param R ratio of undiseased subjects and diseased subjects
#' @return n_one sample size of diseased subjects in a one-sided hypothesis testing
#'         N_one total sample size in a one-sided hypothesis testing
#'         n_two sample size of diseased subjects in a two-sided hypothesis testing
#'         N_two total sample size in a two-sided hypothesis testing
#' @export
samplepartialarea<-function(a,b,e1,e2,alpha,beta,L,R){
   r1<-exp(-a^2/(2*(1+b^2)))
   r2<-1+b^2
   r3<-pnorm((qnorm(e2)+a*b/(1+b^2))*sqrt(1+b^2))-pnorm((qnorm(e1)+a*b/(1+b^2))*sqrt(1+b^2))
   r4<-exp(-((qnorm(e1)+a*b/(1+b^2))*sqrt(1+b^2))^2/2)-exp(-((qnorm(e2)+a*b/(1+b^2))*sqrt(1+b^2))^2/2)
   f<-r1/sqrt(2*pi*r2)*r3
   g<-r1*r4/(2*pi*r2)-a*b*r1*r3/sqrt(2*pi*r2^3)
   V<-f^2*(1+b^2/R+a^2/2)+g^2*(b^2*(1+R)/(2*R))
   n_one<-ceiling((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*V/L^2)
   n_two<-ceiling((qnorm(alpha/2,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2*V/L^2)
   N_one<-ceiling(n_one*(1+R))
   N_two<-ceiling(n_two*(1+R))
   result<-list(n_one=NULL,N_one=NULL,n_two=NULL,N_two=NULL)
   result$n_one<-n_one
   result$N_one<-N_one
   result$n_two<-n_two
   result$N_two<-N_two
   return(result)
}

