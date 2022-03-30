if (!require("gsDesign")) {
  install.packages("gsDesign")
  library("gsDesign")
}


# parameter description ---------------------------------------------------

#' @param p1 treatment group response rate
#' @param p2 control group response rate
#' @param alpha type I error (one sided)
#' @param n_max maximum acceptable sample size
#' @param timing information time of interim analysis
#' @param pz_low lower bound of promising zone
#' @param pz_high upper bound of promising zone
#' @param nsim number of simulations
#' @param k number of looks
#' @param sfu alpha spending function , check package `gsDesign`
#' @param sfupar alpha spending function parameter, check package `gsDesign`
#' @param cp_fu   futility boundary based on conditional power
#' @param scale estimator: "difference", "or", "rr" 
#' @param overrun overrun during interim analysis
#' @return operating characteristics of adaptive design
#' @note planed sample size calculation is based on the best scenario 
# which is that p1 = 0.6  p2 = 0.25 power = 0.86 under two-stage sequential design


pz.sim.binary <- function(p1,
                          p2,
                          alpha,
                          n_max,
                          timing,
                          pz_low,
                          pz_high,
                          nsim,
                          cp_fu,
                          sfu = NULL,
                          sfupar = NULL,
                          k = 2,
                          scale = "difference",
                          seed = 1,
                          overrun = 15,
                          n1 = NULL,
                          n2 = NULL,
                          futility_b = NULL
){
  set.seed(seed)
  scale = match.arg(tolower(scale), c("difference", "rr", "or"))
  
  # size = nBinomial(p1=p1,p2=p2,
  #                  alpha=alpha,beta = beta, sided=1,
  #                  scale = scale,
  #                  n=NULL)
  
  design = gsDesign(k = k,test.type = 1,  alpha = alpha,
                    n.fix= 100, # sample size does not affect boundary calculation
                    timing = c(timing),
                    sfu = sfu,
                    sfupar = 1,
                    endpoint = 'Binomial')
  
  
  #stage 1 boundary
  z1_bound = design$upper$bound[1]
  #stage 2 boundary
  z2_bound = design$upper$bound[2]
  
  #stage 1 sample size
  # if(is.null(n1)){n1 = ceiling(design$n.I)[1]}
  #planned sample size
  # if(is.null(n2)){n2 = ceiling(design$n.I)[2]}
  #maintain equal allocation
  is.even <- function(x){ x %% 2 == 0 }
  if(is.even(n1) == F){n1 = n1 + 1; n2 = n2 + 1}
  if(is.even(n2) == F){n2 = n2 + 1}
  #stage 2 sample size
  n2_tilde = n2 - n1
  #simulate trials
  results = as.data.frame(matrix(NA, nrow=nsim, ncol=18)) 
  colnames(results) <- c('Stage 1 sample size', 'Stage 2 sample size (planned)',
                         'Total sample size (planned)',
                         'Stage 2 sample size (re-estimated)','Total sample size (re-estimated)',
                         'Stage 1 TRT success','Stage 1 TRT failure',
                         'Stage 1 CTRL success','Stage 1 CTRL failure',
                         'Chisq','Wald Z  w./ Cont','CP','Relative risk reduction','Risk diff',
                         'Wald_adaptive', "Success_adaptive",'Wald_sequential',"Success_sequential")
  results[,1] = n1
  results[,2] = n2_tilde
  results[,3] = n2
  #data
  results[,6] = rbinom(nsim, n1/2,p1)#response in treatment group
  results[,7] = n1/2 - results[,6]
  results[,8] = rbinom(nsim, n1/2,p2)#response in control group
  results[,9] = n1/2 - results[,8]
  # results[,13] = (results[,6]/(n1/2) - results[,8]/(n1/2)) / (results[,8]/(n1/2))
  results[,14] = (results[,6]/(n1/2)) - (results[,8]/(n1/2))
  #test statistic in stage 1: z1
  results[,11] = apply(results[,6:9],1,function(x) 
    testBinomial(x[1],x[3],x[1]+x[2],x[3]+x[4],adj=1,scale = scale))
  #CP
  cp <- function(n1,n2,z1){
    n2_tilde = n2-n1
    1 - pnorm(((z2_bound*sqrt(n2) - z1*sqrt(n1)) / sqrt(n2_tilde)) -
                ((z1*sqrt(n2_tilde)) / sqrt(n1)))
  }
  results[,12] = apply(results[,1:11],1,function(x)
    cp(x[1], x[3], x[11]))
  
  
  #re-estimate  sample size
  re.est <- function(pz_low,pz_high,n_max,n1,n2,z1,cp){
    if(cp >= pz_low & cp < pz_high & z1 < z1_bound){
      n2_tilde_prime = (n1/(z1^2))*((z2_bound*sqrt(n2)-z1*sqrt(n1))/
                                      sqrt(n2-n1) + qnorm(pz_high))^2
      n2_prime = n2_tilde_prime + n1
      result = ceiling(max(min(n2_prime, n_max), n2))
      is.even <- function(x){ x %% 2 == 0 }
      if(is.even(result) == F){result = result + 1}
    }else{
      if(cp < cp_fu | z1 >= z1_bound){
        result = n1 + overrun
      }else{
        result = n2
      }
    }
    return(result)
  }
  #sample size(re-estimated)
  results[,5] = apply(results[,1:12],1,function(x)
    re.est(pz_low,pz_high,n_max,x[1], x[3], x[11],x[12]))
  #stage 2 sample size(re-estimated)
  results[,4] = results[,5] - results[,1]
  # final test  in stage 2
  final.test <- function(z1,x1,x2,n2_new,n2,p1,p2,
                         cp, n2_plan = NULL, n_plan = NULL){
    if (is.null(n2_plan)) {
      if (cp < cp_fu & z1 < z1_bound) {
        result = -Inf #if stop for futiltiy 
      }else if(cp >= cp_fu & z1 >= z1_bound){
        result = Inf #if stop for efficacy
      }else{
        trt = rbinom(1, n2_new/2,p1)
        ctrl = rbinom(1, n2_new/2,p2)
        result = testBinomial(x1+trt,x2+ctrl,n2/2,n2/2,
                              adj=1,scale = scale)
      } #cont. to stage 2
      
    }else{
      trt = rbinom(1, n2_plan/2,p1)
      ctrl = rbinom(1, n2_plan/2,p2)
      result = testBinomial(x1+trt,x2+ctrl,n_plan/2,n_plan/2,
                            adj=1,scale = scale)
    } #no change made to sample size
    return(result) 
  }
  # test statistic in stage 2: z2_adaptive
  results[,15] = apply(results[,1:13],1,function(x)
    final.test(x[11], x[6], x[8], x[4], x[5], p1, p2, cp = x[12]))
  results[,16] = as.numeric(results[,15] >= z2_bound) 
  
  wald = rep(0,nsim)
  wald =  apply(results[,1:13],1,function(x){    
    final.test(x[11],x[6], x[8], x[4],x[5],p1,p2,
               cp = x[12],n2_plan = x[2] ,n_plan = x[3])
  })
  
  results[ ,17] = results[ ,15]
  # GSD stage 2 test statistic
  results[results[ ,12] >= pz_low & results[ ,12] < pz_high, 17] = wald[results[ ,12] >= pz_low & results[ ,12] < pz_high]
  # test statistic in stage 2: z2_sequential 
  results[ ,18] = as.numeric(results[,17] >= z2_bound)
  
  #temp = results[results[,12] >= pz_low & results[,12] < pz_high,]
  #library("tidyverse")
  #temp = as.data.frame(temp)
  #temp1 = temp %>%
  #  count(`Total sample size (re-estimated)`,`Stage 1 TRT success`,
  #        `Stage 1 TRT failure`,`Stage 1 CTRL success`,`Stage 1 CTRL failure`,
  #        `CP`)
  #temp2 = temp %>%
  #  count(`Total sample size (re-estimated)`, `CP`)
  #temp2$total = temp2$`Total sample size (re-estimated)`*temp2$n
  #sum(temp2$total)/sum(temp2$n)
  
  cat("\nProb. futility", sum(results[,12] < cp_fu)/nsim)
  cat("\nProb. unfavorable", sum(results[,12] < pz_low)/nsim)
  cat("\nProb. promising", sum(results[,12] >= pz_low & results[,12] < pz_high)/nsim)
  cat("\nProb. favorable", sum(results[,12] >= pz_high)/nsim)
  cat("\n")
  cat("\nAverage sample size (unfavorable)", mean(results[results[,12] < pz_low,5]))
  cat("\nAverage sample size (promising)", mean(results[results[,12] >= pz_low & results[,12] < pz_high,5]))
  cat("\nAverage sample size (favorable)", mean(results[results[,12] >= pz_high,5]))
  cat("\n")
  cat("\n##Sequential Design##")
  cat("\nPower conditional on interim outcome", sum(results[results[,12] < pz_low, 18])/sum(results[,12] < pz_low))
  cat("\nPower conditional on interim outcome", sum(results[results[,12] >= pz_low & results[,12] < pz_high, 18])/sum(results[,12] >= pz_low & results[,12] < pz_high))
  cat("\nPower conditional on interim outcome", sum(results[results[,12] >= pz_high, 18])/sum(results[,12] >= pz_high))
  cat("\nOverall Power ", sum(results[, 18])/nsim)
  cat("\n")
  cat("\n##Adaptive Design##")
  cat("\nPower conditional on interim outcome", sum(results[results[,12] < pz_low, 16])/sum(results[,12] < pz_low))
  cat("\nPower conditional on interim outcome", sum(results[results[,12] >= pz_low & results[,12] < pz_high, 16])/sum(results[,12] >= pz_low & results[,12] < pz_high))
  cat("\nPower conditional on interim outcome", sum(results[results[,12] >= pz_high, 16])/sum(results[,12] >= pz_high))
  cat("\nOverall Power ", sum(results[, 16])/nsim)
  
  
  # return(results)
}

# HS006 igG4-RD -----------------------------------------------------
# 
# 0.60 vs 0.25 n1 = 38 n2 = 74, nmax = 148,CPmin = 0.3556 
# Scenario 1
pz.sim.binary(p1 = 0.60,
              p2 = 0.25,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 1,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 2
pz.sim.binary(p1 = 0.55,
              p2 = 0.25,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 2,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 3
pz.sim.binary(p1 = 0.50,
              p2 = 0.25,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 3,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 4
pz.sim.binary(p1 = 0.45,
              p2 = 0.25,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 4,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 5
pz.sim.binary(p1 = 0.60,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              cp_fu = 0.15,
              nsim = 10000,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 5,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 6
pz.sim.binary(p1 = 0.55,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 6,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 7
pz.sim.binary(p1 = 0.50,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 7,
              overrun = 0,
              n1 = 38,
              n2 = 74)
# Scenario 8
pz.sim.binary(p1 = 0.45,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 8,
              overrun = 0,
              n1 = 38,
              n2 = 74)


# HS006 (RA) TRIAL --------------------------------------------------------
# 0.56 vs 0.36 n1 = 110 n2 = 220 nmax = 330 CP min is  0.4061 
# Scenario 1
pz.sim.binary(p1 = 0.56,
              p2 = 0.36,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 12,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 2
pz.sim.binary(p1 = 0.60,
              p2 = 0.18,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 22,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 3
pz.sim.binary(p1 = 0.55,
              p2 = 0.18,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 32,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 4
pz.sim.binary(p1 = 0.50,
              p2 = 0.18,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 42,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 5
pz.sim.binary(p1 = 0.45,
              p2 = 0.18,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 52,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 6
pz.sim.binary(p1 = 0.60,
              p2 = 0.24,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 62,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 7
pz.sim.binary(p1 = 0.55,
              p2 = 0.24,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 72,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 8
pz.sim.binary(p1 = 0.50,
              p2 = 0.24,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 82,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 9
pz.sim.binary(p1 = 0.45,
              p2 = 0.24,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 92,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 10
pz.sim.binary(p1 = 0.60,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 102,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 11
pz.sim.binary(p1 = 0.55,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 112,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 12
pz.sim.binary(p1 = 0.50,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 122,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 13
pz.sim.binary(p1 = 0.45,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 132,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 14
pz.sim.binary(p1 = 0.60,
              p2 = 0.36,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 142,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 15
pz.sim.binary(p1 = 0.50,
              p2 = 0.36,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 152,
              overrun = 0,
              n1 = 110,
              n2 = 220)
# Scenario 16
pz.sim.binary(p1 = 0.45,
              p2 = 0.36,
              alpha = 0.025,
              n_max = 330,
              timing = 0.5,
              pz_low = 0.4061,
              pz_high = 0.8,
              nsim = 10000,
              cp_fu = 0.15,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 162,
              overrun = 0,
              n1 = 110,
              n2 = 220)



#CP to risk diff
cp.to.diff <- function(x1=NULL, x2=NULL, n1, n2, cp_fu=0.15){
  #
  if(is.null(x1) == F & is.null(x2) == F){
    n_trt = n1/2
    n_ctrl = n1/2
    z = testBinomial(x1 = x1, x2 = x2,
                     n1 = n_trt, n2 = n_ctrl, adj=1, scale = "Difference")
    cp_temp = cp(n1, n2, z)
    cat("\nCP", cp_temp)
  }else{
    #get 2nd stage boundary
    design = gsDesign(k = 2,test.type = 2,  alpha = 0.025,  
                      beta = 0.2, 
                      n.fix= n2, timing = c(0.5), 
                      sfu = sfLDOF,
                      sfupar = 1,
                      endpoint = 'Binomial')
    z2_bound = design$upper$bound[2]
    alpha = 2*(1-pnorm(z2_bound))
    #CP
    cp <- function(n1,n2,z1){
      n2_tilde = n2-n1
      1 - pnorm(((z2_bound*sqrt(n2) - z1*sqrt(n1)) / sqrt(n2_tilde)) -
                  ((z1*sqrt(n2_tilde)) / sqrt(n1)))
    }
    #search for risk diff
    x1_temp = x2
    x2_temp = x1
    n_trt = n1/2
    n_ctrl = n1/2
    cp_temp = -Inf
    while(cp_temp <= cp_fu){
      if(is.null(x1) == F & is.null(x2) == T){
        x2_temp = x2_temp - 1
        z = testBinomial(x1 = x1, x2 = x2_temp,
                         n1 = n_trt, n2 = n_ctrl, adj=1, scale = "Difference")
        cp_temp = cp(n1, n2, z)
      }else if(is.null(x1) == T & is.null(x2) == F){
        x1_temp = x1_temp + 1
        z = testBinomial(x1 = x1_temp, x2 = x2,
                         n1 = n_trt, n2 = n_ctrl, adj=1, scale = "Difference")
        cp_temp = cp(n1, n2, z)
      }
    }
    if(is.null(x1) == F & is.null(x2) == T){
      cat("\nno. success trt", x1, "(", x1/n_trt*100,"%)")
      cat("\nno. success", x2_temp, "(", x2_temp/n_ctrl*100,"%)")
      cat("\nRisk diff", x1-x2_temp, "(", ((x1/n_trt)-(x2_temp/n_ctrl))*100,"%)")
      cat("\nCP", cp_temp)
    }else if(is.null(x1) == T & is.null(x2) == F){
      cat("\nno. success trt", x1_temp, "(", x1_temp/n_trt*100,"%)")
      cat("\nno. success", x2, "(", x2/n_ctrl*100,"%)")
      cat("\nRisk diff", x1_temp-x2, "(", ((x1_temp/n_trt)-(x2/n_ctrl))*100,"%)")
      cat("\nCP", cp_temp)
    }
  }
}
cp.to.diff(x1=11, x2=NULL, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=10, x2=NULL, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=9, x2=NULL, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=8, x2=NULL, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=11, x2=9, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=10, x2=8, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=9, x2=5, n1=38, n2=74, cp_fu=0.15)
cp.to.diff(x1=8, x2=6, n1=38, n2=74, cp_fu=0.15)



cp.to.diff(x1=33, x2=NULL, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=30, x2=NULL, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=28, x2=NULL, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=25, x2=NULL, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=NULL, x2=20, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=NULL, x2=17, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=NULL, x2=13, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=NULL, x2=10, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=33, x2=28, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=33, x2=29, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=32, x2=28, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=21, x2=17, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=18, x2=14, n1=110, n2=220, cp_fu=0.15)
cp.to.diff(x1=17, x2=13, n1=110, n2=220, cp_fu=0.15)


