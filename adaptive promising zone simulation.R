# simulating trials with binomial endpoint to evaluate operating characteristics
<<<<<<< HEAD
# of an adaptive design based on promising zone with futility and efficacy stop
=======
# of a adaptive design based on promising zone with futility and efficacy stop
>>>>>>> 06e01525d6f5953854f215188515d922d59871f2

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
#' @param adjust   Miettinen and Nurminen test statistic This includes a factor of n / (n - 1) where n is the total sample size. If adj is not 1, this factor is not applied. The default is adj=0 since nominal Type I error is generally conservative with adj=1 
#' @param scale estimator: "difference", "or", "rr" 
#' @param overrun overrun during interim analysis
#' @return operating characteristics of adaptive design
#' @note  sample size calculation is based on fixed design

pz.sim.binary <- function(p1,
                          p2,
                          alpha,
                          n_max,
                          timing,
                          pz_low,
                          pz_high,
                          nsim,
                          cp_fu,
                          adjust = 0,
                          sfu = NULL,
                          sfupar = NULL,
                          k = 2,
                          scale = "difference",
                          seed = 1,
                          overrun = 15,
                          n1 = NULL,
                          n2 = NULL
){
  set.seed(seed)
  scale = match.arg(tolower(scale), c("difference", "rr", "or"))
  # efficacy stop boundary
  design = gsDesign(k = k,test.type = 1,  alpha = alpha,
                    n.fix= 100, # sample size does not affect boundary calculation
                    timing = c(timing),
                    sfu = sfu,
                    sfupar = sfupar,
                    endpoint = 'Binomial')
  
  
  #stage 1 boundary
  z1_bound = design$upper$bound[1]
  #stage 2 boundary
  z2_bound = design$upper$bound[2]
  
  is.even <- function(x){ x %% 2 == 0 }
  if(!is.even(n1)){n1 = n1 + 1; n2 = n2 + 1}
  if(!is.even(n2)){n2 = n2 + 1}
  #stage 2 sample size
  n2_tilde = n2 - n1
  #simulate trials
  results = as.data.frame(matrix(NA, nrow=nsim, ncol=18)) 
  colnames(results) <- c('Stage 1 sample size', # 1
                         'Stage 2 sample size (planned)',# 2
                         'Total sample size (planned)',# 3
                         'Stage 2 sample size (re-estimated)',# 4
                         'Total sample size (re-estimated)',# 5
                         'Stage 1 TRT success',# 6
                         'Stage 1 TRT failure',# 7
                         'Stage 1 CTRL success',# 8
                         'Stage 1 CTRL failure',# 9
                         'Chisq',# 10
                         'Wald Z  w./ Cont',# 11
                         'CP',# 12
                         'Relative risk reduction',# 13
                         'Risk diff',# 14
                         'Wald_adaptive', # 15
                         "Success_adaptive",# 16
                         'Wald_sequential',# 17
                         "Success_sequential")# 18
  results[ ,1] = n1
  results[ ,2] = n2_tilde
  results[ ,3] = n2
  #data
  results[ ,6] = rbinom(nsim, n1/2, p1)# response in treatment group
  results[ ,7] = n1/2 - results[ ,6]
  results[ ,8] = rbinom(nsim, n1/2, p2)# response in control group
  results[ ,9] = n1/2 - results[ ,8]
  results[ ,14] = (results[ ,6]/(n1/2)) - (results[ ,8]/(n1/2))
  #test statistic in stage 1: z1
  results[,11] = apply(results[,6:9], 1, function(x) 
    testBinomial(x[1], x[3], x[1]+x[2], x[3]+x[4], adj = adjust, scale = scale))
  #CP
  cp <- function(n1, n2, z1){
    n2_tilde = n2-n1
    1 - pnorm(((z2_bound*sqrt(n2) - z1*sqrt(n1)) / sqrt(n2_tilde)) -
                ((z1*sqrt(n2_tilde)) / sqrt(n1)))
  }
  results[,12] = apply(results[ ,1:11], 1, function(x)
    cp(x[1], x[3], x[11]))
  
  
  #re-estimate  sample size
  re.est <- function(pz_low, pz_high, n_max, n1, n2, z1, cp){
    if(cp >= pz_low & cp < pz_high & z1 < z1_bound){
      n2_tilde_prime = (n1/(z1^2))*((z2_bound*sqrt(n2)-z1*sqrt(n1))/
                                      sqrt(n2-n1) + qnorm(pz_high))^2
      n2_prime = n2_tilde_prime + n1
      result = ceiling(max(min(n2_prime, n_max), n2))

      if(!is.even(result)){result = result + 1}
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
  results[ ,5] = apply(results[ ,1:12],1,function(x)
    re.est(pz_low, pz_high, n_max, x[1], x[3], x[11], x[12]))
  #stage 2 sample size(re-estimated)
  results[ ,4] = results[ ,5] - results[ ,1]
  # final test  in stage 2
  final.test <- function(z1, x1, x2, n2_new, n2, p1, p2,
                         cp, n2_plan = NULL, n_plan = NULL){
    if (is.null(n2_plan)) {
      if (cp < cp_fu & z1 < z1_bound) {
        result = -Inf #if stop for futility 
      }else if(cp >= cp_fu & z1 >= z1_bound){
        result = Inf #if stop for efficacy
      }else{
        trt = rbinom(1, n2_new/2, p1)
        ctrl = rbinom(1, n2_new/2, p2)
        result = testBinomial(x1+trt, x2+ctrl, n2/2, n2/2,
                              adj = adjust, scale = scale)
      } #cont. to stage 2
      
    }else{
      trt = rbinom(1, n2_plan/2, p1)
      ctrl = rbinom(1, n2_plan/2, p2)
      result = testBinomial(x1+trt, x2+ctrl, n_plan/2, n_plan/2,
                            adj = adjust, scale = scale)
    } #no change made to sample size
    return(result) 
  }
  # test statistic in stage 2: z2_adaptive
  results[ ,15] = apply(results[,1:13], 1, function(x)
    final.test(x[11], x[6], x[8], x[4], x[5], p1, p2, cp = x[12]))
  results[ ,16] = as.numeric(results[,15] >= z2_bound) 
  
  wald = rep(0, nsim)
  wald =  apply(results[ ,1:13], 1, function(x){    
    final.test(x[11], x[6], x[8], x[4], x[5], p1, p2,
               cp = x[12], n2_plan = x[2] ,n_plan = x[3])
  })
  
  results[ ,17] = results[ ,15]
  # GSD stage 2 test statistic
  results[results[ ,12] >= pz_low & results[ ,12] < pz_high, 17] = wald[results[ ,12] >= pz_low & results[ ,12] < pz_high]
  # test statistic in stage 2: z2_sequential 
  results[ ,18] = as.numeric(results[ ,17] >= z2_bound)
  
  
  cat("\nProb. futility", sum(results[ ,12] < cp_fu)/nsim)
  cat("\nProb. unfavorable", sum(results[ ,12] < pz_low)/nsim)
  cat("\nProb. promising", sum(results[ ,12] >= pz_low & results[ ,12] < pz_high)/nsim)
  cat("\nProb. favorable", sum(results[ ,12] >= pz_high)/nsim)
  cat("\nProb. efficacy", sum(results[ ,11] >= z1_bound)/nsim)
  cat("\n")
  cat("\nAverage sample size (unfavorable)", mean(results[results[ ,12] < pz_low, 5]))
  cat("\nAverage sample size (promising)", mean(results[results[ ,12] >= pz_low & results[ ,12] < pz_high, 5]))
  cat("\nAverage sample size (favorable)", mean(results[results[ ,12] >= pz_high, 5]))
  cat("\n")
  cat("\n##Sequential Design##")
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] < pz_low, 18])/sum(results[ ,12] < pz_low))
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] >= pz_low & results[ ,12] < pz_high, 18])/sum(results[ ,12] >= pz_low & results[ ,12] < pz_high))
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] >= pz_high, 18])/sum(results[ ,12] >= pz_high))
  cat("\nOverall Power ", sum(results[, 18])/nsim)
  cat("\n")
  cat("\n##Adaptive Design##")
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] < pz_low, 16])/sum(results[ ,12] < pz_low))
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] >= pz_low & results[ ,12] < pz_high, 16])/sum(results[ ,12] >= pz_low & results[ ,12] < pz_high))
  cat("\nPower conditional on interim outcome", sum(results[results[ ,12] >= pz_high, 16])/sum(results[ ,12] >= pz_high))
  cat("\nOverall Power ", sum(results[ ,16])/nsim)
  
  
  # return(results)
}

pz.sim.binary(p1 = 0.60,
              p2 = 0.30,
              alpha = 0.025,
              n_max = 148,
              timing = 0.5,
              pz_low = 0.3556,
              pz_high = 0.8,
              cp_fu = 0.15,
              adjust = 0,
              nsim = 10000,
              k = 2,
              sfu = sfLDOF,
              sfupar = 1,
              scale = "diff",
              seed = 5,
              overrun = 0,
              n1 = 38,
              n2 = 74)

#CP to risk diff
cp.to.diff <- function(x1 = NULL, x2 = NULL, n1, n2, adjust = 0, cp_fu=0.15){
  #
  if(!is.null(x1) & !is.null(x2)){
    n_trt = n1/2
    n_ctrl = n1/2
    z = testBinomial(x1 = x1, x2 = x2,
                     n1 = n_trt, n2 = n_ctrl, adj = adjust, scale = "Difference")
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
    cp <- function(n1, n2, z1){
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
                         n1 = n_trt, n2 = n_ctrl, adj = adjust, scale = "Difference")
        cp_temp = cp(n1, n2, z)
      }else if(is.null(x1) == T & is.null(x2) == F){
        x1_temp = x1_temp + 1
        z = testBinomial(x1 = x1_temp, x2 = x2,
                         n1 = n_trt, n2 = n_ctrl, adj = adjust, scale = "Difference")
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
cp.to.diff(x1=11, x2=NULL, n1=38, n2=74, adjust = 0, cp_fu=0.15)

