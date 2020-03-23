#' Simmple Contingency Table Analysis 
#' 
#' 
#' This function returns simple measures of indepence for I x J contingency tables. Namely, 
#' it gives Pearsons X^2 and the likelihood ratio test statistic, G^2. Both are asymptotically distributed 
#' chi-squared with (I-1) x (J-1) degrees of freedom. It aslo returns pearson standardized residuals for the 
#' contingency table. 
#' 
#' 
#' @author Ben Brennan (brennben at umich dot edu)
#' 
#' @param obs This is a contingency table, or it should be. That is it. If it not - then idk, you shouldn't be using this function?
#' 
#' @return 


contingency_table = function(obs){
  
  ni_plus = rowSums(obs)
  
  n_plusj = colSums(obs)
  
  n = sum(as.vector(obs))
  
  chi = 0
  
  g = 0
  
  rows = NULL
  
  cols = NULL
  
  stand_residuals = matrix(rep(0,length(as.vector(obs))),nrow(obs),ncol(obs))
  
  for(j in 1:ncol(obs)){
    
    for(i in 1:nrow(obs)){
      
      mu = ni_plus[i] * n_plusj[j] / n
      
      chi = chi + (obs[i,j] - mu)^2 / mu
      
      g = g + obs[i,j]*(log(obs[i,j]/mu))
      
      stand_residuals[i,j] = (obs[i,j] - mu) / (sqrt(mu*(1 - ni_plus[i]/n)*(1 - n_plusj[j]/n)))
      
      rows = c(rows,rep(i,obs[i,j]))
      cols = c(cols,rep(j,obs[i,j]))
    }
  }
  g = 2*g
  
  p_chi = 1 - pchisq(chi, (nrow(obs) - 1 ) * (ncol(obs) - 1))
  
  p_g = 1 - pchisq(g, (nrow(obs) - 1 ) * (ncol(obs) - 1))
  
  m = (n - 1)*cor(rows,cols)^2
  
  return(list("chi" = chi, "g2" = g, "p_chi" = p_chi, "p_g2" = p_g, 
              "pearson_residuals" = stand_residuals, 
              "m" = m))
}
