MGLMpredInt <- function(object){
  dist <- object$dist
  beta <- as.vector(object$coeff)
  H <- object$Hessian
  vcov <- chol2inv(chol(-H)) 
  Design_mat <- object$data$X
  old_count <- object$data$Y
  
  #1. beta ~ MVN(Beta, vcov(beta))
  x <- rmvnorm(n=1,mean=beta, sigma=vcov)
  sampling_beta <- matrix(x,ncol=ncol(object$coeff), nrow=nrow(object$coeff))# 
  #2 Calc pi or alpha.
  if(dist=="Multinomial"){
    sampling_fitted <-  exp(Design_mat%*%sampling_beta)/(rowSums(exp(Design_mat%*%sampling_beta))+1)
    sampling_fitted <- cbind(sampling_fitted, (1-rowSums(sampling_fitted)))
  #3. sample from Multinom(n;pi)
    new_count <- rmultinomial(rowSums(old_count), sampling_fitted) # Sample from multinomial  
  } else if(dist=="Dirichlet Multinomial"){
    sampling_fitted <- exp(Design_mat%*%sampling_beta)
  #3. sample from DirMultinom(n;alpha)
    new_count <-rdirm(rowSums(old_count), sampling_fitted)
    
  }
  #4. Calc proportions of sampled counts.
  sampling_val  <- new_count/rowSums(old_count)
  return(sampling_val) 
}
