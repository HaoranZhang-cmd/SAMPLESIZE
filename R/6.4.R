#' SAMPLE SIZE FOR DETERMINING A SUITABLE CUTOFF VALUE
#' @param SP the minimum specificity required at the cutoff
#' @param SE_con the conjectured sensitivity at the cutoff
#' @param SE_min the sensitivity at which the cutoff value for the test would no longer be useful
#' @param b the binormal parameter b for the test result
#' @param alpha The choice of alpha determines an overall confidence probability of 1-alpha for the statement "the specificity at the chosen cutoff is at least SP"
#' @param beta The choice of 1-beta determines the probability of finding such a cutoff, when such a cutoff truly exists
#' @return N the rough total sample size when the cutoff point is not determined based on a prespecified specificity and sensitivity
#'         N_pre the precise total sample size when the cutoff point is determined based on a prespecified specificity and sensitivity
#' @export
samplecutoff<-function(SP,SE_con,SE_min,b,alpha,beta){
  lambda<-qnorm(SE_con)-qnorm(SE_min)
  vx<-b*sqrt(1+qnorm(SP)^2/2)
  vy<-sqrt(1+qnorm(SE_min)^2/2)
  calculateN<-function(a){
    if(a<=0||a>=1){
      return(Inf)
    }
    else{
      return((qnorm(sqrt(1-alpha))*(vx/sqrt(a)+vy/sqrt(1-a))+qnorm(1-beta)*sqrt(vx^2/a+vy^2/(1-a)))^2/lambda^2)
    }
  }
  N<-(ceiling(optim(0.5,calculateN,method="BFGS")$value))
  calculateN_pre<-function(a){
    if(a<=0||a>=1){
      return(Inf)
    }
    else{
      return((-sqrt(2)*qnorm(sqrt(alpha/2))*(vx/sqrt(a)+vy/sqrt(1-a))+qnorm(1-beta)*sqrt(vx^2/a+vy^2/(1-a)))^2/lambda^2)
    }
  }
  N_pre<-(ceiling(optim(0.5,calculateN_pre,method="BFGS")$value))
  result<-list(N=NULL,N_pre=NULL)
  result$N<-N
  result$N_pre<-N_pre
  return(result)
}

