###########################################################
##### several functions
###########################################################
h_alpha = function(alpha1, alpha2, eta){
  if (eta > 0){
    if (alpha1 > 0){
      res = (exp(alpha1 * eta) - 1) / alpha1;
    }		
    else if (alpha1 == 0){
      res = eta;
    }
    else{ ## alpha1 < 0
      res = -log(1 - alpha1 * eta) / alpha1;
    }
  }
  else{ 
    if (alpha2 > 0){
      res = -(exp(-alpha2 * eta) - 1) / alpha2;
    }		
    else if (alpha2 == 0){
      res = eta;
    }		
    else{ 
      res = log(1 + alpha2 * eta) / alpha2;
    }
  }
  return(res)
}

###########################################################
##### generate response Y
Generate_Y_2PGlogit <- function(a,b,alpha1,alpha2,theta){
  n_student  <- length(theta)
  n_item<- length(a)
  Y <- matrix(NA,n_student,n_item)
  for(i in 1:n_student){
    for(j in 1:n_item){
      eta_ij= a[j]*(theta[i]-b[j])
      h_ij=h_alpha(alpha1[j],alpha2[j],eta_ij)
      prob <- exp(h_ij)/(1+exp(h_ij))
      if(is.nan(prob)){
        Y[i,j]=1
      }else{
        Y[i,j] <- rbinom(1,1,prob)
      }
    }
  }
  return(Y)
}

###########################################################
##### generate reponse time logT
Generate_logT <- function(tao,lambda,sigma){
  n_student <- length(tao)
  n_item <- length(lambda)
  logT <- matrix(NA,n_student,n_item)
  for(i in 1:n_student){
    for(j in 1:n_item){
      logT[i,j] <- rnorm(1,lambda[j]-tao[i],sigma)
    }
  }
  return(logT)
}

