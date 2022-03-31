# simulating to evaluate operating characteristic of i3+3


CreData <- function(nDose, cDosenames = paste("dose", 1:nDose, sep = " ")){
  data <- data.frame(dose = 1:nDose, 
                     npt = rep(0, nDose), 
                     ndlt = rep(0, nDose), 
                     dosenames = cDosenames)
  return(data)
}
Updateda <- function(data, nPt, nDlt, nLastdose){
  if (!(nLastdose %in% data$dose)) {
    stop("dose do not treat")
  }
  idx <- which(data$dose == nLastdose)
  data[idx, "npt"] <- data[idx, "npt"] + nPt
  data[idx, "ndlt"] <- data[idx, "ndlt"] + nDlt
  return(data) 
}

Pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1) {
    return(x)
  }
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) {
      break
    }
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

Decision <- function(data = data, lastdose, dTarget,
                     dEiu, dEil, nPt, nCohort, nStop) 
{
  pt <- sum(data$npt)
  idx <- which(data$dose == lastdose)
  curpt <- data$npt[idx]
  ndlt <- data$ndlt[idx]
  npt <- data$npt[idx]
  EI <- c(dTarget - dEil, dTarget + dEiu)
  
  if (ndlt/npt < EI[1]) {
    nextdose <-  lastdose + 1
  }
  else if (ndlt/npt >= EI[1] && ndlt/npt <= EI[2]) {
    nextdose <- lastdose
  }
  else if (ndlt/npt > EI[2]) {
    if ((ndlt - 1)/npt < EI[1]) {
      nextdose <- lastdose
    } 
    else if ((ndlt - 1)/npt > EI[1]) {
      nextdose <- lastdose - 1
    }
  }
  if (!(nextdose %in% data$dose) | pt >= nPt * nCohort | curpt >= nStop) {
    nextdose <- NA
  }
  return(nextdose)
}


Si3p3 <- function(vTrue, dTarget, dEiu, 
                  dEil, nPt, nCohort, 
                  nStop)
{
  nextdose <- 1
  ndose <- length(vTrue)
  data <- CreData(ndose)
  
  while (nextdose %in% data$dose) {
    lastdose <- nextdose
    ndlt <- sum(runif(nPt) < vTrue[lastdose])
    data <- Updateda(data, nPt, ndlt, lastdose)
    nextdose <- Decision(data, lastdose, dTarget,
                         dEiu, dEil, nPt, nCohort, nStop)
  }
  y <- data$ndlt 
  n <- data$npt 
  phat = (y + 0.05)/(n + 0.1)
  p.var = (y + 0.05) * (n - y + 0.05)/((n + 0.1)^2 * (n + 0.1 + 1))
  phat = Pava(phat, wt = 1/p.var)
  selectd = sort(abs(phat - dTarget), index.return = T)$ix[1]
  
  selectdose = data$dosename[selectd]
  
  list(data = data, mtd = selectdose, lastdose = lastdose)
}


#' @vTrue   a vector of assumed true toxicity rate of each dose
#' @dTarget target toxicity rate 
#' @dEiu    equivalence interval upper bound
#' @dEil    equivalence interval upper bound
#' @nPt     number of patients treated in each cohort 
#' @nCohort number of cohorts
#' @nStop   stop criteria: total number of patients treated under one dose 
#' @seed    seed 
#' @sim     number of simulation

SSi3p3 <- function(vTrue, dTarget = 0.3, 
                   dEiu = 0.05, dEil = 0.05, 
                   nPt = 3, nCohort = 10, 
                   nStop = 12, seed = 1, 
                   sim = 1000){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  pat <- dlt <- matrix(0,sim,length(vTrue))
  MTD <- numeric(sim)
  for (l in 1:sim) {
    onesimu <- Si3p3(vTrue, dTarget, 
                     dEiu, dEil, 
                     nPt, nCohort, 
                     nStop)
    pat[l, ] <- onesimu$data$npt
    dlt[l, ] <- onesimu$data$ndlt
    MTD[l] <- onesimu$mtd
  }
  
  d <- factor(MTD, levels = paste("dose", 1:length(vTrue), sep = " "))
  
  list(doselevel = 1:length(vTrue),
       truetox = vTrue,
       choose = unname(table(d)/sim),
       patients = colMeans(pat),
       tox = colMeans(dlt))
}

SSi3p3(vTrue =  c(0.1,0.25,0.3,0.4,0.5), seed = 2022)
