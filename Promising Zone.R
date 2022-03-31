## Promising Zone (Two arms, two-stage)
## Reference: Mehta & Pocock 2011

library("ggplot2")

#Eq. 6 (To calculate conditional power at interim look)
CP <- function(z1, n2_tilde, n1, alpha = 0.05, sides = 2){
  n2 = n1 + n2_tilde
  if(sides == 2){alpha = alpha / 2} #alpha means alpha/2
  left = (qnorm(1 - alpha)*sqrt(n2) - z1*sqrt(n1)) / sqrt(n2_tilde)
  right = (z1*sqrt(n2_tilde)) / sqrt(n1)
  result = 1 - pnorm(left-right)
  if(result < 1 & result > 0){
    return(result) 
  }else{cat("ERROR, CP = ", result)}
}
# CP(1.275874, n2_tilde = 74-38, 38)

#Eq.6 (To calculate critical value at stage 1 given promising zone)
Z1 <- function(cp_min, n2_tilde, n1, alpha = 0.05, sides = 2){
  n2 = n1 + n2_tilde
  if(sides == 2){alpha = alpha / 2} #alpha means alpha/2
  up = (qnorm(1 - alpha)*sqrt(n2)/sqrt(n2_tilde)) - qnorm(1 - cp_min)
  down = (sqrt(n1)/sqrt(n2_tilde)) + (sqrt(n2_tilde)/sqrt(n1))
  result = up / down
  return(result)
}
# Z1(0.05, n2_tilde = 74-38, 38)
#Eq. 7 ~ Eq. 9 (To estimate new sample size for stage 2)
N2_STAR <- function(z1, n2_tilde, n1, n_max, power, alpha = 0.05, sides = 2){
  n2 = n1 + n2_tilde
  if(sides == 2){alpha = alpha / 2} #alpha means alpha/2
  n2_tilde_prime = (n1/(z1^2))*((qnorm(1-alpha)*sqrt(n2)-z1*sqrt(n1))/sqrt(n2-n1) + qnorm(power))^2
  n2_prime = n2_tilde_prime + n1
  #n2_prime = ceiling(n2_prime)
  if(sides == 1){
    temp = CP(z1, n2_prime, n1, alpha, sides)
  }else if(sides == 2){
    temp = CP(z1, n2_prime, n1, alpha*2, sides)
  } #need to recover the value of alpha
  result = max(min(n2_prime, n_max), n2) ####NEED A CONSTRAINT
  #cat("The conditional power is:", temp)
  return(result)
}
#Eq. 11 (To estimate the new critical value)
B <- function(z1, n1, n2_tilde, n2_star, alpha = 0.05, sides = 2){
  n2 = n1 + n2_tilde
  n2_tilde_star = n2_star - n1
  if(sides == 2){alpha = alpha / 2} #alpha means alpha/2
  result = ((n2_star)^(-0.5))*
    (sqrt(n2_tilde_star/n2_tilde)*(qnorm(1-alpha)*sqrt(n2)-z1*sqrt(n1))+z1*sqrt(n1))
  return(result)
}


#Decision boundary
#' @CP_min  Lower bound of promising zone
#' @Z1      The critical value needed to achieve the conditional power at interim look (Given Z1, n1 & n2)
#' @N2_STAR The adjusted total sample size (Given Z1)
#' @B       Modified critical value for final analysis
#' @DIFF    B - Z_alpha
#' 
decision <- function(n1, n2_tilde, n_max, power = 0.9, alpha=0.05, sides=2){
  num = power*10000 - 1
  sol = matrix(NA, nrow = num, ncol = 5)
  colnames(sol) = c("CP_min","Z1","N2_STAR","B","DIFF")
  sol[,1] = rep(1:num)/10000
  for(i in 1:num){
    sol[i,2] = Z1(cp_min=sol[i,1], n2_tilde, n1, alpha, sides)
    sol[i,3] = N2_STAR(sol[i,2], n2_tilde, n1, n_max, power, alpha, sides)
    sol[i,4] = B(sol[i,2], n1, n2_tilde, sol[i,3], alpha, sides)
  }
  sol = as.data.frame(sol)
  #plot(sol[,4]~sol[,1])
  if(sides == 2){alpha = alpha / 2}
  sol[,5] = round(abs(sol[,4] - qnorm(1-alpha)),4)
  lower = min(which(sol[,5] == min(sol[,5])))/10000
  cat("CP min is ",lower)
  # return(sol)
}

# decision(38,74-38,150,alpha = 0.02442549,sides = 1)
