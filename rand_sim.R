#############################################################
##                                           
##  Program name     : 
##  Project          : 
##  Written by       : 
##  Date of creation : 
##  Description      : simulating to evaluate treatment imbalance
##  input parameters : 
#' @param nsim, number of simulation
#' @param pat_num, number of patients enrolled
#' @param prob1, enrollment capacity between sites
#' @param prob2, disease status
#############################################################

if (!require(tidyverse)) {
  install.packages(tidyverse, dependencies = TRUE)
  library(tidyverse)
}

set.seed(1)
# loading Dummy rand list
library(readxl)
rand_list <- read_excel("dummy list.xlsx")
view(rand_list)

D_list <- rand_list %>% 
  group_by(STRATA) %>% 
  mutate(id = 1:n()) 

head(D_list)

# XX sites & X disease status
site_num <- unlist(unique(D_list[,"RANDR1"]))
diease_num <- unlist(unique(D_list[,"RANDR2"]))


rand_simulation <- function(nsim, 
                            pat_num, 
                            prob1, 
                            prob2) {
  # save difference
  count <- as.data.frame(matrix(NA, nrow = nsim, ncol = 5))
  colnames(count) <- c("Trt", "One-sided", "Two-sided" ,"None", "Sex")
  # save  subjects enrolled
  df <-as.data.frame(matrix(NA, nrow = pat_num, ncol = 4)) 
  colnames(df) <- c("Pat.id", "RANDR1", "RANDR2","Sex")
  # simulation
  for (i in 1:nsim) {
    df[ ,1] <- 1:pat_num
    df[ ,2] <- sample(site_num,pat_num, replace = TRUE, prob1)
    df[ ,3] <- sample(diease_num, pat_num, replace = TRUE, prob2)
    df[ ,4] <- sample(1:2, pat_num, replace = TRUE)
    suppressMessages(
      sumry <- df %>% group_by(RANDR1, RANDR2) %>% 
        arrange(Pat.id) %>% 
        mutate(id = 1:n()) %>% 
        inner_join(D_list) %>% 
        ungroup())
    
    df1 <- sumry %>% count(RANDARM)
    df2 <- sumry %>% count(RANDARM, RANDR2)
    df3 <- sumry %>% count(RANDARM, Sex)
    # number diff between two arm
    diff1 <- abs(df1[1,2]-df1[2,2])
    # disease status
    # one-sided
    diff2 <- abs(df2[1,3]-df2[4,3])
    # two-sided
    diff3 <- abs(df2[2,3]-df2[5,3])
    # none
    diff4 <- abs(df2[3,3]-df2[6,3])
    # hypothetical  variable XX
    diff5 <- abs(df3[1,3]-df3[3,3])
    
    count[i, ] <-unname(unlist(c(diff1, diff2, diff3, diff4, diff5)))
  }
  apply(count, 2, function(x)quantile(x,c(0.99, 0.95, 0.90, 0.80, 0.70, 0.60)))
}


rand_simulation(nsim = 5000,
                pat_num = 212,
                prob1 = c(60/212, rep(30/212, 3), rep(6/212,5), rep(4/212,8)),
                prob2 = 1:3/6)





