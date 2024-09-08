# Functions used for Weibull and Gumbel inferences and simulations : ####
# proposed in "Practical and Optimal Likelihood Intervals and Regions
# for Weibull and Gumbel Distributions". 
# Tech23_068 ####
# Updated by E. Díaz-Francés on September 7, 2024
#----------------------- --

# Overall arguments and variables used: ####
# Y is the vector of Weibull data
# X is vector of logdata, X=ln(Y).
# T is sum of logdata, T=sum(X)
# N is sample size, N=length(Y)
# P is the probability of the quantile of interest
# AMLE is the mle of location parameter A
#           (equal to the logarithm of the Weibull scale parameter)
# BETAMLE is the mle of Weibull shape parameter BETA
# LVMAX is the maximum value of the loglikelihood at the mles
# INFOA is the individual Fisher's observed information of A
# INFOBETA is the individual Fisher's observed information of BETA
# ______________________________________####

##  Functions to calculate Likelihoods and Mles of Weibull parameters ####

#  This function evaluates the log likelihood of (a,beta) at 
#   the point VEC2=c(A,BETA).
lvabetaWei = function(VEC2, X, T, N) {
  A = VEC2[1]
  BETA = VEC2[2]
  if (BETA > 0) {
    k1 = sum(exp(BETA * X))
    sali = N * log(BETA) + (BETA * T) - (N * A * BETA) -
      exp(-A * BETA) * k1
  } else {
    sali = -999999999999
  }
  return(sali)
}

#------------------------------------- --
# Relative Likelihood R(a,beta) evaluated at point (A,BETA):
Rabeta = function(VEC2, X, T, N, LVMAX) {
  sali = lvabetaWei(VEC2, X, T, N) - LVMAX
  return(exp(sali))
}

#------------------------------------- --
#   Computing the first derivative of the log profile likelihood of beta
#   equated to zero. Its root provides the mle of beta.
#   Both mles are calculated in function Weibullmles() below.
Weibetamle = function(BETA, X, T, N) {
  if (BETA > 0) {
    k1 = sum(exp(BETA * X))
    k2 = sum(X * exp(BETA * X))
    sali = (N / BETA) + T - (N * k2 / k1)
  } else {
    sali = 5555
  }
  return(sali)
}

#------------------------ --
# Mles of Weibull parameters: a,beta, and sigma=exp(a).
Weibullmles = function(Y) {
  N = length(Y)
  X = log(Y)
  T = sum(X)
  # The Gumbel method of moments estimate of Weibull shape beta is:
  betamom = pi * sqrt(N) / sqrt(6 * sum((X - (T / N))**2))
  # The mle of beta is calculated here:
  BETAMLE = uniroot(Weibetamle, c(betamom / 10, (betamom + 20)), X, T, N,
    lower = betamom / 10, upper = (betamom + 20), tol = 0.000001
  )$root
  #  Calculating the mle of a Gumbel:
  k1 = sum(exp(BETAMLE * X))
  AMLE = log(k1 / N) / BETAMLE
  Weimles = list(AMLE, BETAMLE, exp(AMLE))
  names(Weimles) = c("AGumbel", "BetaShape", "SigmaScale")
  return(Weimles)
}

#----------------------------- --

# Weibull distribution of parameters (a,beta) is evaluated at point YY:

pWeibullmin = function(Y, vec2) {
  X = log(Y)
  agum = vec2[1] # Gumbel location of logdata
  betawei = vec2[2] # Weibull shape parameter
  if (betawei > 0) {
    Fwei = 1 - exp(-exp((X - agum) * betawei))
  } else {
    Fwei = 999
  }
  return(Fwei)
}


#---------------------------------- --
## For Weibull model validation the following 2 functions are needed:   ####
#----------------------------- --
# 1. Anderson Darling Modified statistic by Stephens (1977) for testing
# the Weibull distribution

AndersonDarlingModif = function(Y, AMLE, BETAMLE) {
  N = length(Y)
  ZZ = sort(pWeibullmin(Y, c(AMLE, BETAMLE)))
  ns = seq(1, N, by = 1)
  aux1 = (2 * ns - 1)
  aux2 = (2 * N + 1 - (2 * ns))
  AD = -N - (1 / N) * sum((aux1 * log(ZZ)) + aux2 * log(1 - ZZ))
  AD = AD * (1 + 0.2 / sqrt(N))
  return(AD)
}

#----------------------------- --
# 2. Checking if the Weibull distribution is reasonable for Data:
Weibullvalidation = function(Data) {
  Y=Data
  ysort = sort(Y)
  if (ysort[1] <= 0) {
    print("Error: All data must be strictly positive for")
    print("a Weibull distribution to be reasonable.")
  } else {
    X = log(ysort)
    N = length(X)
    H = Weibullmles(Y)
    AMLE = H$AGumbel
    BETAMLE = H$BetaShape
    ModAD = AndersonDarlingModif(Y, AMLE, BETAMLE)
    if (ModAD <= 0.637) {
      print("According to Stephens (1977) Modified Anderson Darling Statistic,")
      print("the Weibull is a reasonable model for the data.")
    } else if ((ModAD > 0.637) && (ModAD < 0.757)) {
      print("By Stephens (1977) Modified Anderson Darling Statistic,")
      print("the Weibull is NOT a reasonable model for the data,")
      print("at a significance level B, such that 0.05 < B < 0.10.")
    } else {
      print("By Stephens (1977) Modified Anderson Darling Statistic,")
      print("the Weibull is NOT a reasonable model for the data,")
      print("at a significance level smaller than 0.05. ")
    }
    #-------------------- --
    # Weibull-Gumbel probability plot:
    # Transformed ordered data with Weibull distribution
    uwei = 1 - exp(-exp(BETAMLE * (X - AMLE)))
    enes = seq(1, N, by = 1)
    aenes = N + 1 - enes
    # Theoretical standard uniform quantiles:
    Ebetas = enes / (N + 1)
    ## 95% approximate confidence band for the PP plot:
    tau = 1 - (.05 / N)
    tau1 = (1 - tau) / 2
    tau2 = (1 + tau) / 2
    Ban1 = qbeta(tau1, enes, aenes)
    Ban2 = qbeta(tau2, enes, aenes)
    #----------------- --
    plot(Ebetas, uwei,
      pch = 19, cex = .5, ylab = "F(x) Transform",
      xlab = "Uniform (0,1) Quantiles", xlim = c(0, 1), ylim = c(0, 1),
      main = "Weibull Probability Plot", pty = "m", col = 1
    )
    lines(c(0, 1), c(0, 1), lty = 2, col = 1)
    lines(Ebetas, uwei, lty = 2, col = 1)
    lines(Ebetas, Ban1, lty = 2, col = 1)
    lines(Ebetas, Ban2, lty = 2, col = 1)
    #------------------ --
    Weivalid = list(AMLE, BETAMLE, ModAD)
    names(Weivalid) = c("Amle", "Betamle", "ModAD")
    return(Weivalid)
  }
}
# ______________________________________####

## Fisher's Observed Information Matrix and individual informations ####
# of single parameters are computed here:  
InformationABeta = function(YY, AMLE, BETAMLE) {
  N = length(YY)
  X = log(YY)
  k1 = sum(exp(BETAMLE * X))
  k2 = sum(X * exp(BETAMLE * X))
  k3 = sum((X**2) * exp(BETAMLE * X))
  aux = (2 * AMLE * k2) - k3 - ((AMLE**2) * k1)
  # negative of second derivatives of l(a,beta) evaluated at mles:
  D2a = (BETAMLE**2) * exp(-AMLE * BETAMLE) * k1
  D2beta = (N / (BETAMLE**2)) - exp(-AMLE * BETAMLE) * aux
  D2abeta = N - exp(-AMLE * BETAMLE) * (k1 * (1 - AMLE * BETAMLE) +
    BETAMLE * k2)
  Infomat = matrix(0, nrow = 2, ncol = 2)
  Infomat[1, 1] = D2a
  Infomat[2, 2] = D2beta
  Infomat[1, 2] = D2abeta
  Infomat[2, 1] = D2abeta
  #  Determinant of Fisher's Information Matrix of (a,beta):
  DETERM = Infomat[1, 1] * Infomat[2, 2] - (Infomat[1, 2] * Infomat[2, 1])
  #  Individual Observed Fisher's Information of a and beta:
  Infoa = DETERM / D2beta
  Infobeta = DETERM / D2a
  Infolist = list(Infomat, Infoa, Infobeta)
  names(Infolist) = c("Infabeta", "Inf_a", "Inf_beta")
  return(Infolist)
}

#-------------------------------- --
## Likelihood levels for intervals and regions are computed here ####
# as a function of sample size N, for confidence levels 90,95, and 99% :
clikelevels = function(N) {
  # For Single parameter profile likelihood interval:
  cs90 = 0.2585 - (1 / (0.775 + (1.641 * N)))
  cs95 = 0.1465 - (1 / (2.757 + (2.082 * N)))
  cs99 = 0.0362 - (1 / (18.294 + (5.229 * N)))
  ccessig = c(cs90, cs95, cs99)
  names(ccessig) = c("90%", "95%", "99%")
  # Two Parameters for likelihood region (beta,sigma)
  ctous90 = 0.10 - (1 / (1.11 + (4.655 * N)))
  ctous95 = 0.05 - (1 / (1.61 + (7.27 * N)))
  ctous99 = 0.01 - (1 / (19.51 + (23.59 * N)))
  ctous = c(ctous90, ctous95, ctous99)
  names(ctous) = c("90%", "95%", "99%")
  cesn = list("csig" = ccessig, "cbetasig" = ctous)
  return(cesn)
}

# ______________________________________####
## Computing Profile likelihood functions and PL intervals ####
#  ----------------------- --
#  Computing at BETA the first partial derivative of loglikelihood
#  l(a,beta) with respect to beta for a fixed value of a, denoted here as AA.
Der1BetLabeta = function(BETA, X, T, N, AA) {
  if (BETA > 0) {
    k1 = sum(exp(BETA * X))
    k2 = sum(X * exp(BETA * X))
    sali = (N / BETA) + T - (N * AA) + exp(-AA * BETA) * ((AA * k1) - k2)
  } else {
    sali = 999999999999
  }
  return(sali)
}

#--------------------------------------------- --
# Computing R(a), the relative profile likelihood of A,
#  numerically at a single value of A, denoted here as AA :
RpA = function(AA, X, T, N, LVMAX, BETAMLE) {
  BETARMLE = uniroot(Der1BetLabeta, c(BETAMLE / 10, (BETAMLE + 20)), X,
    T, N, AA,
    lower = BETAMLE / 10,
    upper = (BETAMLE + 20), tol = 0.000001
  )$root
  RPA = exp(lvabetaWei(c(AA, BETARMLE), X, T, N) - LVMAX)
  sali = list("Rpas" = RPA, "Betarmle" = BETARMLE)
  return(sali)
}

#------------------------------ --
# Computing R*(a),the proposed approximating function to R(a),
#  at single value AA :
RpStarA = function(AA, X, N, AMLE, BETAMLE, INFOA) {
  # auxilliary constants:
  H1n = (0.025 + (0.7 / (N**0.25))) * BETAMLE
  H2n = (0.18 - (1.6 / (N**0.54))) * BETAMLE
  aux1 = -INFOA / 2
  if (N < 10) {
    print("N must be larger than 10")
    rpa = -14
  } else if ((N < 101) && (AA <= AMLE)) {
    rpaleft = aux1 * ((1 - exp(H1n * (AA - AMLE)))**2) / (H1n**2)
    rpa = rpaleft
  } else if ((N < 101) && (AA > AMLE)) {
    rparight = aux1 * ((1 - exp(H2n * (AA - AMLE)))**2) / (H2n**2)
    rpa = rparight
  } else { # N>100
    rpa = aux1 * ((AMLE - AA)**2)
  }
  return(exp(rpa))
}

#------------------------------- --
## PL Intervals for A and SIGMA of confidence levels 90,95, and 99%
# These are the proposed profile likelihood and confidence intervals.
# The likelihood interval of level c=0.001 helps to determine plotting range.
IntervalsA = function(N, AMLE, BETAMLE, INFOA) {
  ces = clikelevels(N)$csig
  cess = c(ces, 0.001)
  #  Proposed profile likelihood intervals for a are calculated here:
  H1n = (0.0247 + (0.7002 / (N**0.2534))) * BETAMLE
  H2n = (0.1767 - (1.6333 / (N**0.5374))) * BETAMLE
  auxA = sqrt(-2 * log(cess) / INFOA)
  #----------------- --
  if (N < 10) {
    cat("WARNING: Sample size n must be at least 10")
  } else if ((N >= 10) && (N <= 100)) {
    InterA = cbind(AMLE + log(1 - H1n * auxA) / H1n, AMLE +
      log(1 + H2n * auxA) / H2n)
    AA1 = InterA[4, 1]
    AA2 = InterA[4, 2]
    InterSigma = exp(InterA)
  } else { # (only when n>100 the following is calculated:)
    InterA = cbind(AMLE - auxA, AMLE + auxA)
    AA1 = InterA[4, 1]
    AA2 = InterA[4, 2]
    InterSigma = exp(InterA)
  }
  intsali = list(
    "SigmaIntervals" = InterSigma[1:3, 1:2],
    "AIntervals" = InterA[1:3, 1:2],
    "LikeLevels" = ces, "EndpointsA" = c(AA1, AA2)
  )
  return(intsali)
}

#--------------------------------------- --
## Computing R(beta), the relative profile likelihood of beta,
#  at the single value BETA

RpbetaWei = function(BETA, X, T, N, BETAMLE) {
  if (BETA > 0) {
    k1 = sum(exp(BETA * X))
    k1MV = sum(exp(BETAMLE * X))
    sali = N * log(BETA / BETAMLE) + T * (BETA - BETAMLE) +
      N * log(k1MV / k1)
    sali = exp(sali)
  } else {
    sali = 0
  }
  return(sali)
}

#--------------------------------------- --
## Computing R*(beta), the proposed approximating function to R(beta),
# at the single value BETA:
RpbetaStWei = function(BETA, BETAMLE, INFOBETA) {
  aux = sqrt(BETA / BETAMLE)
  rstbet = -2 * INFOBETA * (BETAMLE^2) * ((1 - aux)^2)
  return(exp(rstbet))
}

#--------------------------------------------------- --
## Computing PL Intervals for BETA and B=1/BETA,
# of confidence levels 90, 95, and 99%
# as well as their likelihood interval of level 0.001 for plotting purposes:
IntervalsBeta = function(N, BETAMLE, INFOBETA) {
  # likelihood levels for 90,95, and 99% confidence levels are:
  ces = clikelevels(N)$csig
  cess = c(ces, 0.001)
  aux = sqrt(-log(cess) / (2 * INFOBETA)) / BETAMLE
  auxmat = cbind((1 - aux)**2, (1 + aux)**2)
  #   Profile likelihood intervals for beta:
  InterBeta = BETAMLE * auxmat
  # Plotting limits for Rp(beta):
  BBet1 = InterBeta[4, 1]
  BBet2 = InterBeta[4, 2]
  # Intervals for Weibull Shape Beta in first 2 columns of BetaMatrix
  BetaMatrix = InterBeta[1:3, ]
  #---------------------------- --
  #   Profile Likelihood intervals for Gumbel b:
  InterB = cbind(1 / InterBeta[1:3, 2], 1 / InterBeta[1:3, 1])
  # Intervals for Gumbel scale b:
  intsali = list(
    "BetaIntervals" = BetaMatrix, "BIntervals" = InterB,
    "LikeLevels" = cess[1:3], "EndpointsBeta" = c(BBet1, BBet2)
  )
  return(intsali)
}


# ______________________________________####
# Quantiles Mles, Weibull and Gumbel, of probability P   ####
#-------------------------------------------- --
# Weibull Quantiles mles:

WeiQpmle = function(P, AMLE, BETAMLE) {
  WP = log(-log(1 - P))
  Qweip = exp(AMLE + WP / BETAMLE)
  return(Qweip)
}

# Gumbel Quantiles mles:
GumQpmle = function(P, AMLE, BETAMLE) {
  WP = log(-log(1 - P))
  Qgump = (AMLE + WP / BETAMLE)
  return(Qgump)
}

#--------------------------------------------------- --
# Quantiles: Profile Likelihoods and  Intervals  ####
# ----------------------- --

# Computing the Weibull loglikelihood parametrized in terms of beta
# and the Gumbel quantile QGP of probability P.
lpQpGum = function(BETA, X, T, N, P, QGP) {
  HP = log(-log(1 - P))
  if (BETA > 0) {
    k1 = sum(exp(BETA * X))
    lv = (N * log(BETA)) + BETA * (T - (N * QGP)) -
      k1 * exp(-(BETA * QGP) + HP) + (N * HP)
  } else {
    lv = -999999999999
  }
  return(lv)
}
#-------------------------- --

# Computing the relative profile likelihood function of the Gumbel 
# quantile QGp at point QGPP. For this purpose, the log likelihood lpQpGum 
# is maximized over BETA at a fixed value QGPP.
RpQGp = function(QGPP, P, X, T, N, LVMAX, BETAMLE) {
  SALIS = optimize(lpQpGum, c(BETAMLE / 10, (BETAMLE + 20)), X, T, N,
    P, QGPP,
    lower = BETAMLE / 10,
    upper = (BETAMLE + 20), tol = 0.000001,
    maximum = TRUE
  )
  # The restricted mle of beta given QGPP is:
  BETARMLE = SALIS$maximum
  rpQGP = (SALIS$objective) - LVMAX
  resul = list("betarmle" = BETARMLE, "RpQ" = exp(rpQGP))
  return(resul)
}

#  ----------------------- --
#  Computing function [ rp(QGp)-lnCC] in order to find numerically
#  profile likelihood intervals of QGp having a likelihood level CC.
InterQp = function(QGPP, X, T, N, LVMAX, P, BETAMLE, QGpMLE, CC) {
  SALIS = optimize(lpQpGum, c(BETAMLE / 10, (BETAMLE + 20)), X, T, N,
    P, QGPP,
    lower = BETAMLE / 10,
    upper = (BETAMLE + 20), tol = 0.000001,
    maximum = TRUE
  )
  # BETARMLE=SALIS$maximum
  rpQGPP = (SALIS$objective) - LVMAX
  return(rpQGPP - log(CC))
}

#  ----------------------- --
# Computing profile likelihoo intervals of the Gumbel quantile QGP of
# probability P having confidence levels 90, 95, and 99%
IntervalQP = function(P, Y, N, AMLE, BETAMLE) {
  X = log(Y) # log data XX are Gumbel
  T = sum(X)
  # Gumbel quantile MLE:
  Wpp = log(-log(1 - P))
  QG = AMLE + Wpp / BETAMLE
  # Gumbel quantile of very smalle probability:
  Wpini = -35 # PP = 6.66 e-16
  QGini = AMLE + Wpini / BETAMLE
  # Gumbel quantile of large probability:
  Wpsup = log(-log(1 - 0.9999999))
  QGsup = AMLE + Wpsup / BETAMLE
  # The maximum loglikelihood at the mle is used to calculate them:
  LVMAX = lvabetaWei(c(AMLE, BETAMLE), X, T, N)
  # Likelihood levels of intervals:
  ces = clikelevels(N)$csig
  c90 = ces[1]
  c95 = ces[2]
  c99 = ces[3]
  csmall = 0.01
  #----------------- --
  QGp90inf = uniroot(InterQp, c(QGini, QG), X, T, N, LVMAX, P,
                     BETAMLE, QG, c90,
                     lower = QGini, upper = QG, tol = 0.000001
  )$root
  QGp90sup = uniroot(InterQp, c(QG, QGsup), X, T, N, LVMAX, P,
                     BETAMLE, QG, c90,
                     lower = QG, upper = QGsup, tol = 0.000001
  )$root
  #  ----------------- --
  QGp95inf = uniroot(InterQp, c(QGini, QG), X, T, N, LVMAX, P,
                     BETAMLE, QG, c95,
                     lower = QGini, upper = QG, tol = 0.000001
  )$root
  QGp95sup = uniroot(InterQp, c(QG, QGsup), X, T, N, LVMAX, P,
                     BETAMLE, QG, c95,
                     lower = QG, upper = QGsup, tol = 0.000001
  )$root
  #  --------------------  --
  QGp99inf = uniroot(InterQp, c(QGini, QG), X, T, N, LVMAX, P,
                     BETAMLE, QG, c99,
                     lower = QGini, upper = QG, tol = 0.000001
  )$root
  QGp99sup = uniroot(InterQp, c(QG, QGsup), X, T, N, LVMAX, P,
                     BETAMLE, QG, c99,
                     lower = QG, upper = QGsup, tol = 0.000001
  )$root
  #  --------------------  --
  QGp1 = uniroot(InterQp, c(QGini, QG), X, T, N, LVMAX, P,
                 BETAMLE, QG, csmall,
                 lower = QGini, upper = QG, tol = 0.000001
  )$root
  QGp2 = uniroot(InterQp, c(QG, QGsup), X, T, N, LVMAX, P,
                 BETAMLE, QG, csmall,
                 lower = QG, upper = QGsup, tol = 0.000001
  )$root
  # --------------------- --
  Intizq = c(QGp90inf, QGp95inf, QGp99inf)
  names(Intizq) = c("90%", "95%", "99%")
  Intder = c(QGp90sup, QGp95sup, QGp99sup)
  # Gumbel quantiles:
  IntervalsQGp = cbind(Intizq, Intder)
  # Weibull quantiles:
  IntervalsQWp = exp(IntervalsQGp)
  # -------------- --
  intsali = list(
    "QGpmle" = QG, "QWpmle" = exp(QG),
    "GumbelQPIntervals" = IntervalsQGp,
    "WeibullQPIntervals" = IntervalsQWp,
    "LikeLevels" = ces, "EndpointsQG" = c(QGp1, QGp2)
  )
  return(intsali)
}

# ______________________________________####
# Principal function for computing Weibull inferences and plots ####
# for parameters and the quantile of interest. 
# All data must be positive and the ModAD test statistic must be 0 < ModAD < 1.
# The function WeibullValidation() calculates ModAD and must be run previously.

Weibullinferences = function(Data, QpProb = 0.50) {
  # The probability of the quantile Qp of interest is specified in QpProb.
  # As a default, inferences for the median are calculated,
  # corresponding to QpPRob= 0.50.
  Y = Data
  N = length(Y)
  P = QpProb
  ysort = sort(Y)
   if ((ysort[1] < 0) || (P <= 0) || (P >= 1)) {
    print("Warning: Note that all data must be positive since otherwise")
    print("the Weibull distribution is not reasonable. ")
    print("Check too that a value within (0,1) was given as QpProb, the")
    print("probability of the quantile of interest .")
    ModAD = 99
  } else {
    # Mles of A and BETA and Stephens' statistic are computed here:
    X = log(ysort)
    T = sum(X)
    H = Weibullmles(Y)
    AMLE = H$AGumbel
    BETAMLE = H$BetaShape
    ModAD = AndersonDarlingModif(Y, AMLE, BETAMLE)
    ModADQuantiles = c(0.637, 0.757, 0.877, 1.038)
    names(ModADQuantiles) = c("Q90", "Q95", "Q975", "Q99")
  }

  if (ModAD <= 1) {
    # If ModAD is smaller than one, the Weibull distribution is reasonable
    # so inferences and plots of relative likelihoods are calculated here.
    # ----------------- --
    # Maximum value of the loglikelihood at the mles:
    LVmax = lvabetaWei(c(AMLE, BETAMLE), X, T, N)
    #-------------- --
    # Fisher's Observed Informations used to calculate intervals:
    B2 = InformationABeta(Y, AMLE, BETAMLE)
    # Information of a:
    Infoa = B2$Inf_a
    # Information of beta:
    Infobeta = B2$Inf_beta
    # Observed information matrix of (a,beta):
    Infoabeta = B2$Infabeta
    #--------------- --
    ## INFERENCES FOR BETA:
    # Computing PL intervals for beta, R(beta) and R*(beta)
    #  of likelihood levels c90, c95, c99 and c1=0.001
    B3 = IntervalsBeta(N, BETAMLE, Infobeta)
    Betaint = B3$BetaIntervals
    Betaend = B3$EndpointsBeta
    cs = B3$LikeLevels
    c90 = cs[1]
    c95 = cs[2]
    c99 = cs[3]
    Nbet = 150
    betis = seq(Betaend[1], Betaend[2], length.out = Nbet)
    # R*(beta) Approximating function fpr R(beta)
    RpbetaSt = RpbetaStWei(betis, BETAMLE, Infobeta)
    # Exact relative profile likelihood of beta, R(beta),
    # for comparison
    Rpbetas = rep(0, Nbet)
    for (i in (1:Nbet)) {
      BETA = betis[i]
      Rpbetas[i] = RpbetaWei(BETA, X, T, N, BETAMLE)
    }
    #-------------------- --
    ## INFERENCES FOR A AND SIGMA:
    B4 = IntervalsA(N, AMLE, BETAMLE, Infoa)
    Aend = B4$EndpointsA
    SigmaInt = B4$SigmaIntervals
    AInt = B4$AIntervals
    # Setting range for plot of Rp(a) and Rp*(a) that includes amle:
    nhalf = 100
    As1 = seq(Aend[1], AMLE, length.out = nhalf)
    difi = As1[2] - As1[1]
    As2 = seq((AMLE + difi), Aend[2], length.out = nhalf)
    Ais = c(As1, As2)
    Sigis = exp(Ais)
    nk = length(Ais)
    #--- Calculating Rp(a) and Rp*(a) over plot range here:
    RpaExact = rep(0, nk)
    BetarmleA = rep(0, nk)
    RpaSt = rep(0, nk)
    for (i in 1:nk) {
      AA = Ais[i]
      B4 = RpA(AA, X, T, N, LVmax, BETAMLE)
      RpaExact[i] = B4$Rpas
      BetarmleA[i] = B4$Betarmle
      RpaSt[i] = RpStarA(AA, X, N, AMLE, BETAMLE, Infoa)
    }
    ## ----------------------------- --
    # SIMULTANEOUS INFERENCES FOR (A,BETA):
    # Calculating likelihood regions:
    Ccontours = clikelevels(N)$cbetasig
    NC = 50 # number of points in grid
    labetamat = matrix(rep(0, NC * NC), nrow = NC, ncol = NC)
    # Values for axis limits for contour plot:
    AA1 = Aend[1]
    AA2 = Aend[2]
    AAS = seq(AA1, AA2, length = NC)
    BBet1 = Betaend[1]
    BBet2 = Betaend[2]
    BETAS = seq(BBet1, BBet2, length = NC)
    # Log relative likelihood r(a,beta) is evaluated over grid:
    for (i in 1:NC) {
      for (j in 1:NC) {
        bivec = c(AAS[i], BETAS[j])
        labetamat[i, j] = lvabetaWei(bivec, X, T, N)
      }
    }
    Rabetamat = exp(labetamat - LVmax) # R(a,beta) calculated over grid
    Labelvec = c("90%", "95%", "99%") # contour labels
    #------------------------ --
    # Inferences for the Weibull quantile QWp of probability P:
    IntQWP = IntervalQP(P, Y, N, AMLE, BETAMLE)
    #  Mles of Gumbel and Weibull quantiles of probability P:
    QGPmle = IntQWP$QGpmle
    QWPmle = IntQWP$QWpmle
    # Weibull PL intervals for QWp:
    IntsQWP = IntQWP$WeibullQPIntervals
    QpLikelevels = IntQWP$LikeLevels
    # Gumbel PL intervals for QGp:
    IntsQGP = log(IntsQWP)
    # endpoints for plotting R(Qp):
    QGPend = IntQWP$EndpointsQG
    # Calculating the quantile relative profile likelihood for plotting:
    NQp = 50
    qgums = seq(QGPend[1], QGPend[2], length.out = NQp)
    qweis = exp(qgums)
    RQweis = rep(0, NQp)
    RQgums = rep(0, NQp)
    for (i in (1:NQp)) {
      RelG = RpQGp(qgums[i], P, X, T, N, LVmax, BETAMLE)$RpQ
      RQgums[i] = RelG
    }
    #-------------------------- --
    # EXPORTING inferences to folder WeibullResults:
    dir.create("./WeibullResults", showWarnings = FALSE)
    sink(file = "WeibullResults/Weibull_Inferences_Output.txt")
    cat("Weibull inferences results for the provided data vector Y ")
    cat("\n")
    cat("are given here, and are equivalent to Gumbel inferences")
    cat("\n")
    cat("for the log data X=lnY: ")
    cat("\n")
    cat("____________________________________________")
    cat("\n")
    cat("\n")
    cat("The sample size was:")
    cat(N)
    cat("\n")
    cat("The mles of the Weibull shape and scale parameters, beta,sigma, were:")
    cat("\n")
    cat(BETAMLE)
    cat("\n")
    cat("\n")
    cat(exp(AMLE))
    cat("\n")
    cat("The mles of Gumbel location and scale parameters, a,b,")
    cat("\n")
    cat(" for log data X were:")
    cat("\n")
    cat(AMLE)
    cat("\n")
    cat(1 / BETAMLE)
    cat("\n")
    cat("The maximum value of the loglikelihood (a,beta) at the mle was:")
    cat("\n")
    cat(LVmax)
    cat("\n")
    cat("------------------------------------")
    cat("\n")
    cat("The observed Fisher's informations for a and beta were:")
    cat("\n")
    cat(Infoa)
    cat("\n")
    cat(Infobeta)
    cat("\n")
    cat("Fisher's observed information matrix for (a,beta) was:")
    cat("\n")
    print(Infoabeta)
    cat("\n")
    cat("------------------------------------")
    cat("\n")
    cat("The Modified Anderson Darling Statistic (ModAD) proposed by ")
    cat("\n")
    cat("Stephens (1977) for testing if the Weibull distribution ")
    cat("\n")
    cat("is reasonable for the data was:")
    cat("\n")
    cat(ModAD)
    cat("\n")
    cat("to be compared with the following quantiles of its distribution")
    cat("\n")
    cat("when the Weibull is reasonable:")
    cat("\n")
    print(ModADQuantiles)
    cat("\n")
    cat("WARNING: If ModAD is larger than these percentiles, then")
    cat("\n")
    cat("the Weibull/Gumbel models are not reasonable for the data/logdata")
    cat("\n")
    cat("at the corresponding significance level.")
    cat("\n")
    cat("-------------------------------------")
    cat("\n")
    cat("\n")
    cat("The PL intervals for the Weibull shape parameter beta were:")
    cat("\n")
    print(Betaint)
    cat("\n")
    cat("The PL intervals for the Weibull scale parameter sigma were:")
    cat("\n")
    cat("Note that sigma is also the Weibull quantile of probability 0.632:")
    cat("\n")
    print(SigmaInt)
    cat("\n")
    cat("The PL intervals for the Gumbel location parameter a=ln(sigma) were:")
    cat("\n")
    print(log(SigmaInt))
    cat("\n")
    cat("------------------------------------")
    cat("\n")
    cat("For the quantile QP of probability P:")
    cat("\n")
    cat(P)
    cat("\n")
    cat("the Mles of the Gumbel and Weibull Quantiles are:")
    cat("\n")
    cat(paste(round(QGPmle, 4), ",", round(QWPmle, 4)))
    cat("\n")
    cat("The confidence intervals for the Weibull quantile Qp are:")
    cat("\n")
    print(IntsQWP)
    cat("\n")
    cat("------------------------------------")
    cat("\n")
    cat("The likelihood levels of intervals of a single parameter ")
    cat("\n")
    cat("depend on the sample size n and :")
    cat("\n")
    cat("for the following confidence levels were:")
    cat("\n")
    print(QpLikelevels)
    cat("\n")
    cat("The likelihood levels of regions for simultaneous estimation")
    cat("\n")
    cat("of both parameters (a,beta) for this sample size were:")
    cat("\n")
    print(Ccontours)
    cat("\n")
    cat("------------------------------------")
    sink(file = NULL)
    #-------------------- --
    ## EXPORTING relevant PLOTS TO FOLDER WeibullPlots:
    dir.create("./WeibullPlots", showWarnings = FALSE)
    #------------------------------- ---
    # Plot of R*(beta) and R(beta) wiht PL intervals:
    # R(beta) in solid line, R*(beta) in dashes
    #------------------------------------ --
    png("WeibullPlots/Rbeta.png", width = 600, height = 350)
    plot(betis, RpbetaSt,
      type = "l", ylim = c(0, 1), lty = 2,
      main = expression(paste("R(", beta, ")")), lwd = 2,
      ylab = expression(paste("R(", beta, ")")),
      xlab = expression(paste("Weibull shape ", beta)),
      yaxs = "i", col = 2
    )
    # lines(c(1, 1), c(0, 1), lty = 2, col = 1, lwd = 2)
    lines(betis, Rpbetas, type = "l", lty = 1, col = 1, lwd = 2)
    # PL intervals :
    lines(Betaint[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
    lines(Betaint[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
    lines(Betaint[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)
    dev.off()
    #------------------------------------- --
    # Plot of R(sigma), R*(sigma), and PL intervals:
    # R(sigma) in solid line, R*(sigma) in dashes:
    png("WeibullPlots/Rsigma.png", width = 600, height = 350)
    plot(Sigis, RpaExact,
      type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
      xlim = exp(c(Aend[1], Aend[2])),
      main = expression(paste("R(", sigma, ")")),
      xlab = expression(paste("Weibull scale ", sigma, "=QW(0.63)")),
      ylab = expression(paste("R(", sigma, ")")),
      yaxs = "i"
    )
    # R*(sigma) in dashes:
    lines(Sigis, RpaSt, lty = 2, col = 1, lwd = 2)
    #----- PL Confidence intervals for sigma:
    lines(SigmaInt[1, ], rep(c90, 2), lty = 1, col = 2, lwd = 2)
    lines(SigmaInt[2, ], rep(c95, 2), lty = 1, col = 2, lwd = 2)
    lines(SigmaInt[3, ], rep(c99, 2), lty = 1, col = 2, lwd = 2)
    dev.off()
    #------------------------------------- --
    # Plotting R(a),R*(a), and PL intervals.
    # R(a) in solid line, R*(a) in dashes:
    png("WeibullPlots/RA.png", width = 600, height = 350)
    plot(Ais, RpaExact,
      type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
      xlim = c(Aend[1], Aend[2]),
      main = expression(paste("R(", a, ")")), xlab = expression(a),
      ylab = expression(paste("R(", a, ")")),
      yaxs = "i"
    )
    # R*(a) in dashes:
    lines(Ais, RpaSt, lty = 2, col = 1, lwd = 2)
    #----- PL Confidence intervals for a:
    lines(AInt[1, ], rep(c90, 2), lty = 1, col = 2, lwd = 2)
    lines(AInt[2, ], rep(c95, 2), lty = 1, col = 2, lwd = 2)
    lines(AInt[3, ], rep(c99, 2), lty = 1, col = 2, lwd = 2)
    dev.off()
    #------------------------------------- --
    # Computing Likelihood Regions of R(a,beta)
    # of confidence levels 90, 95, and 99%,
    png("WeibullPlots/LikeRegionsABeta.png", width = 600, height = 350)
    contour(AAS, BETAS, Rabetamat,
      levels = Ccontours, labels = Labelvec,
      labcex = 1, xlim = c(AA1, AA2), ylim = c(BBet1, BBet2),
      xlab = expression(a), ylab = expression(beta), cex.lab = 1, cex.axis = 1,
      main = expression(paste("Likelihood Regions of R(", a, ",", beta, ")"))
    )
    # The mles are marked with an asterisk over the contours:
    points(AMLE, BETAMLE, col = 1, pch = 8)
    dev.off()
    #--------------------------------- --
    # Plotting R(QWp), relative prof like. of Weibull quantile of prob. P:
    png("WeibullPlots/RWeiQuantile.png", width = 600, height = 350)
    shortP=round(P,digits=4)
    plot(qweis, RQgums,
      type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
      xlim = c(qweis[1], qweis[NQp]),
      main = paste("Profile Likelihood of Weibull Quantile of Probabilty",
                   shortP),
      xlab = expression(QWp),
      ylab = expression("R(QWp)"),
      yaxs = "i"
    )
    #----- PL Confidence intervals for Weibull Qp:
    lines(IntsQWP[1, ], rep(c90, 2), lty = 1, col = 2, lwd = 2)
    lines(IntsQWP[2, ], rep(c95, 2), lty = 1, col = 2, lwd = 2)
    lines(IntsQWP[3, ], rep(c99, 2), lty = 1, col = 2, lwd = 2)
    dev.off()
    #---------------------------- --
    # Plotting Weibull Probability Plot
    #-------------------- --
    # Weibull-Gumbel probability plot:
    # Transformed ordered data with Weibull distribution
    uwei = 1 - exp(-exp(BETAMLE * (X - AMLE)))
    enes = seq(1, N, by = 1)
    aenes = N + 1 - enes
    # Theoretical standard uniform quantiles:
    Ebetas = enes / (N + 1)
    ## 95% approximate confidence band for the P plot:
    tau = 1 - (.05 / N)
    tau1 = (1 - tau) / 2
    tau2 = (1 + tau) / 2
    Ban1 = qbeta(tau1, enes, aenes)
    Ban2 = qbeta(tau2, enes, aenes)
    #----------------- --
    png("WeibullPlots/WeibullProbPlot.png", width = 600, height = 350)
    plot(Ebetas, uwei,
      pch = 19, cex = .5, ylab = "F(x) Transform",
      xlab = "Uniform (0,1) Quantiles", xlim = c(0, 1), ylim = c(0, 1),
      main = "Weibull Probability Plot", pty = "m", col = 1
    )
    lines(c(0, 1), c(0, 1), lty = 2, col = 1)
    lines(Ebetas, uwei, lty = 2, col = 1)
    lines(Ebetas, Ban1, lty = 2, col = 1)
    lines(Ebetas, Ban2, lty = 2, col = 1)
    dev.off()
    #----------------------------- --
    print("For results and plots, please check contents of folders ")
    print("WeibullResults and WeibullPlots, located in the working directory.")
    print(paste(
      "The mles of a, beta, sigma, and Weibull Qp of prob.",
      P, "are:"
    ))
    WeiMLES = list(AMLE, BETAMLE, exp(AMLE), QWPmle)
    names(WeiMLES) = c("Amle", "Betamle", "Sigmamle", "WeibullQPmle")
    return(WeiMLES)
    #------------"------------------ --
  } else {
    print("Finally, check the value of ModAD given here")
    print(ModAD)
    print("If this is 99, then you made one of the above mentioned errors")
    print("when providing your input to this function.")
    print("However, if 1 < ModAD < 99, then the Weibull distribution")
    print("is not reasonable for the data according to Stephens' statistic.")
    print("Consequently, no inferences or plots were computed.")
  } # Second IF cicle ending
} # WeibullInferences Function ending


# ______________________________________####
# Computing Exponential relative likelihoods of mean and quantile ####
# These were used for Example 4.1

# Exponential relative likelihood of the mean theta:
RthetaExp = function(THETAS, YY) {
  N = length(YY)
  TT = sum(YY)
  tetamle = TT / N
  Rtet = (tetamle^N) * exp(-(TT / THETAS) + N) / (THETAS^N)
  return(Rtet)
}
#--------------------------------- --
# Exponential relative likelihood of quantile QP:
RQpExp = function(QPS, YY, PP) {
  N = length(YY)
  tetamle = sum(YY) / N
  TETIS = QPS / (-log(1 - PP))
  Rqp = RthetaExp(TETIS, YY)
  return(Rqp)
}

# __________________________ ####
# Replicating simulations presented in the article: ####

# Simulating Gumbel of minima samples 
# of size N, with A0=mode and scale B0 = 1/BETA0
#  where Beta0 is the Weibull shape parameter.
rGumin = function(N, mode, scale) {
  # This is equivalent to the log of a Weibull sample.
  X = mode + (scale) * log(-log(1 - runif(N)))
  xsort = sort(X) # ordered logdata to save time later
  return(xsort)
}

## Computing Coverage frequencies of considered intervals ####
#   and regions in simulations:

#---------------------------- --
# Computing Agresti-Coull's interval of 95% confidence for a Binomial
# proportion when V is equal to the number of successes
# and M is the number of  Bernoulli trials.
# This interval is used to estimate the coverage probability of each type
# of interval considered in the article for estimating Weibull parameters
# and quantiles. As used here, V below is the number of intervals that covered
# the true value of the parameter out of the M simulated Gumbel samples.
# The output is a matrix of two columns containing the AC interval
# in each row that corresponds to the coverage frequency of each entry of V.

AgrestiCoul = function(V, M) {
  zp = qnorm(0.975, mean = 0, sd = 1)
  z2 = zp^2
  phat = (V + (z2 / 2)) / (M + z2)
  SAC = zp * sqrt((phat * (1 - phat)) / (M + z2))
  ACinterval = cbind(phat - SAC, phat + SAC)
  return(ACinterval)
}


#----------------------- --
# The following function calculates the coverage frequencies of
# the considered intervals and also of those intervals that under
# and over covered the true value of the parameter. In addition,
# the coverage frequencies of likelihood regions proposed and asymptotic,
# as well as of the proposed and asymptotic profile likelihood intervals for
# a quantile of interest of probability PP are also calculated here.

CoverFreq = function(A0, BETA0, N, P, M) {
  # Number of simulated samples is M
  # The true values of Weibull parameters are BETA0 (Weibull shape)
  # and A0 (log of Weibull scale)
  # The theoretical Gumbel quantile of probability P is:
  QGp0 = A0 + (log(-log(1 - P))) / BETA0
  #------------------------ --
  # The standard Normal Quantile of prob. (1-alfa/2)=0.975 will be used:
  Znor = qnorm(0.975, mean = 0, sd = 1)
  # The Asymptotic likelihood levels for 95% confidence intervals and regions:
  casympt = 0.1465
  casympt2 = 0.05
  # The proposed likelihood levels for likelihood intervals with 95% confidence:
  ces = clikelevels(N)$csig
  c95 = ces[2]
  # For both parameters, the likelihood level of the likelihood region
  # having 95% confidence level is:
  ctous = clikelevels(N)$cbetasig
  ct95 = ctous[2]
  #----------------- --
  # Simulation results will be stored at:
  Rpa0s = rep(0, M) # R(a) at trie va;ie
  RpEDFa0s = rep(0, M) # R*(a) at true value
  Rpbeta0s = rep(0, M) # R(beta) at true value
  RpEDFbeta0s = rep(0, M) # R*(beta) at true value
  RpQp0s = rep(0, M) # R(QGp) at true value
  Rabeta0s = rep(0, M) # R(a,beta) at true values
  amles = rep(0, M) # mle of a
  betamles = rep(0, M) # mle of beta
  # (WA) intervals of 95% confidence level:
  IntWald1_a0 = rep(0, M) # Left endpoint of WA interval for a
  IntWald2_a0 = rep(0, M) # Right endpoint of WA interval for a
  IntWald1_beta0 = rep(0, M) # Left endpoint of WA interval for beta
  IntWald2_beta0 = rep(0, M) # Right endpoint of WA interval for beta
  #------------------- --
  for (i in (1:M)) {
    # A sample of size N of Gumbel data is simulated here.
    # This is equivalent to the log of a Weibull sample.
    X = A0 + (1 / BETA0) * log(-log(1 - runif(N)))
    xsort = sort(X) # ordered logdata to save time later
    # The corresponding ordered Weibull data are:
    Y = exp(xsort)
    X = xsort
    T = sum(X)
    #----- mles
    B = Weibullmles(Y)
    # Mle of a, Gumbel location parameter:
    amle = B$AGumbel
    amles[i] = amle
    # Mle of beta, Weibull shape parameter:
    betamle = B$BetaShape
    betamles[i] = betamle
    # Maximum value of the loglikelihood at the mle:
    LVmax = lvabetaWei(c(amle, betamle), X, T, N)
    #----- infos:
    B = InformationABeta(Y, amle, betamle)
    Infoa = B$Inf_a
    Infobeta = B$Inf_beta
    # Profile relative likelihoods at true values here:
    #-----R(a0) and R*(a0):
    Rpa0s[i] = RpA(A0, X, T, N, LVmax, betamle)$Rpas
    RpEDFa0s[i] = RpStarA(A0, X, N, amle, betamle, Infoa)
    #-----R(beta0) and R*(beta0):
    Rpbeta0s[i] = RpbetaWei(BETA0, X, T, N, betamle)
    RpEDFbeta0s[i] = RpbetaStWei(BETA0, betamle, Infobeta)
    # ---- R(Qp0):
    RpQp0s[i] = RpQGp(QGp0, P, X, T, N, LVmax, betamle)$RpQ
    # (WA) Wald type intervals for a and beta:
    aux1 = Znor / sqrt(Infoa)
    aux2 = Znor / sqrt(Infobeta)
    IntWald1_a0[i] = amle - aux1
    IntWald2_a0[i] = amle + aux1
    IntWald1_beta0[i] = betamle - aux2
    IntWald2_beta0[i] = betamle + aux2
    # Relative likelihood R(a,beta) at true values:
    Rabeta0s[i] = Rabeta(c(A0, BETA0), X, T, N, LVmax)
  }
  # Computing coverage frequencies of PL, WA, and APL intervals
  # and counting those intervals that under and overestimated true parameters:
  #------------------ --
  # Intervals for A0:
  # coverage freqs of PL intervals for A0:
  CovEDFA = sum(RpEDFa0s >= c95)
  # under and over estimated:
  UnderEDFA = sum((RpEDFa0s < c95) * (amles < A0))
  OverEDFA = sum((RpEDFa0s < c95) * (amles > A0))
  #------------------ --
  # Cov freqs. of APL Intervals for A0:
  CovRpA = sum(Rpa0s >= c95)
  # under and over estimated:
  UnderRpA = sum((Rpa0s < c95) * (amles < A0))
  OverRpA = sum((Rpa0s < c95) * (amles > A0))
  #-------------------- --
  # Cov freqs. of WA Intervals for A0:
  CovWA = sum((IntWald1_a0 <= A0) * (IntWald2_a0 >= A0))
  # under and over estimated:
  UnderWA = sum((IntWald2_a0 < A0))
  OverWA = sum((IntWald1_a0 > A0))
  #---------------------------- --
  # Coverage freqs. for Beta:
  # coverage freqs. of PL intervals for  Beta:
  CovEDFBeta = sum(RpEDFbeta0s >= c95)
  # under estimated:
  UnderEDFBeta = sum((RpEDFbeta0s < c95) * (betamles < BETA0))
  OverEDFBeta = sum((RpEDFbeta0s < c95) * (betamles > BETA0))
  #------------------ --
  # Cov freqs of APL Intervals for Beta:
  CovRpBeta = sum(Rpbeta0s >= casympt)
  # under estimated:
  UnderRpBeta = sum((Rpbeta0s < casympt) * (betamles < BETA0))
  OverRpBeta = sum((Rpbeta0s < casympt) * (betamles > BETA0))
  #-------------------- --
  # Cov feqs of WA Intervals for Beta:
  CovWBeta = sum((IntWald1_beta0 <= BETA0) * (IntWald2_beta0 >= BETA0))
  # under and over estimated:
  UnderWBeta = sum((IntWald2_beta0 < BETA0))
  OverWBeta = sum((IntWald1_beta0 > BETA0))
  #-----------------------  --
  # Coverage freqs for quantile Gumbel QGp:
  # Profile likelihood intervals for QGp with c95:
  CovQp = sum(RpQp0s >= c95)
  # APL Intervals for QGp:
  CovQpAsympt = sum(RpQp0s >= casympt)
  #---------------------------- --
  # Coverage freqs for likelihood regions of R(a,beta):
  CovRegAbeta = sum(Rabeta0s >= ct95)
  # Cov freqs of regions with Asymptotic likelihood level:
  CovRegAsympt = sum(Rabeta0s >= 0.05)
  #---------------------------------------- --
  # Generating Output:
  AcovPLs = c(UnderEDFA, CovEDFA, OverEDFA)
  AcovAPLs = c(UnderRpA, CovRpA, OverRpA)
  AcovWAs = c(UnderWA, CovWA, OverWA)
  #------ --
  BetacovPLs = c(UnderEDFBeta, CovEDFBeta, OverEDFBeta)
  BetacovAPLs = c(UnderRpBeta, CovRpBeta, OverRpBeta)
  BetacovWAs = c(UnderWBeta, CovWBeta, OverWBeta)
  #------ --
  covQPs = c(CovQp, CovQpAsympt)
  covRegs = c(CovRegAbeta, CovRegAsympt)
  #------ --
  sali = list(
    "AcovPL" = AcovPLs, "AcovAPL" = AcovAPLs, "AcovWA" = AcovWAs,
    "BetacovPL" = BetacovPLs, "BetacovAPL" = BetacovAPLs, 
    "BetacovWA" = BetacovWAs,
    "CovQGp" = covQPs, "CovReg" = covRegs
  )
  return(sali)
}  # END of function CoverFreq

#--------------------- --
# The following are functions that print the coverage frequencies and their AC
# intervals for the considered intervals: PL, APL, and WA for A and for BETA.
# In addition, the coverage frequencies of Qp intervals and regions(a,beta),
#  with their corresponding Agresti-Coull's intervals are also given.

# Coverage frequencies results for intervals for a:
CovResulA = function(HHH, M) {
  # PL For a:
  CovPLa = HHH$AcovPL
  AC_PLa = AgrestiCoul(CovPLa, M)
  OutputPLa = matrix(cbind(CovPLa, AC_PLa), nrow = 3, ncol = 3)
  colnames(OutputPLa) = c("CovFreqs", "PL_a_Lower", "PL_a_Upper")
  rownames(OutputPLa) = c("a Under", "a Covered", "a Over")
  # APL For a:
  CovAPLa = HHH$AcovAPL
  AC_APLa = AgrestiCoul(CovAPLa, M)
  OutputAPLa = matrix(cbind(CovAPLa, AC_APLa), nrow = 3, ncol = 3)
  colnames(OutputAPLa) = c("CovFreqs", "APL_a_Lower", "APL_a_Upper")
  rownames(OutputAPLa) = c("a Under", "a Covered", "a Over")
  # WA For a:
  CovWAa = HHH$AcovWA
  AC_WAa = AgrestiCoul(CovWAa, M)
  OutputWAa = matrix(cbind(CovWAa, AC_WAa), nrow = 3, ncol = 3)
  colnames(OutputWAa) = c("CovFreqs", "WA_a_Lower", "WA_a_Upper")
  rownames(OutputWAa) = c("a Under", "a Covered", "a Over")
  sali = list(
    "OutputPLa" = OutputPLa, "OutputAPLa" = OutputAPLa,
    "OutputWAa" = OutputWAa
  )
  return(sali)
}

#-----------------------  --
# Coverage frequencies results for intervals for beta and their
# Agresti-Coull's intervals for estimating these proportions:

CovResulBeta = function(HHH, MM) {
  # PL For a:
  CovPLbeta = HHH$BetacovPL
  AC_PLbeta = AgrestiCoul(CovPLbeta, MM)
  OutputPLbeta = matrix(cbind(CovPLbeta, AC_PLbeta), nrow = 3, ncol = 3)
  colnames(OutputPLbeta) = c("CovFreqs", "PL_beta_Lower", "PL_beta_Upper")
  rownames(OutputPLbeta) = c("Beta Under", "Beta Covered", "Beta Over")
  # APL For a:
  CovAPLbeta = HHH$BetacovAPL
  AC_APLbeta = AgrestiCoul(CovAPLbeta, MM)
  OutputAPLbeta = matrix(cbind(CovAPLbeta, AC_APLbeta), nrow = 3, ncol = 3)
  colnames(OutputAPLbeta) = c("CovFreqs", "APL_beta_Lower", "APL_beta_Upper")
  rownames(OutputAPLbeta) = c("Beta Under", "Beta Covered", "Beta Over")
  # WA For a:
  CovWAbeta = HHH$BetacovWA
  AC_WAbeta = AgrestiCoul(CovWAbeta, MM)
  OutputWAbeta = matrix(cbind(CovWAbeta, AC_WAbeta), nrow = 3, ncol = 3)
  colnames(OutputWAbeta) = c("CovFreqs", "WA_beta_Lower", "WA_beta_Upper")
  rownames(OutputWAbeta) = c("Beta Under", "Beta Covered", "Beta Over")
  sali = list(
    "OutputPLbeta" = OutputPLbeta,
    "OutputAPLbeta" = OutputAPLbeta, "OutputWAbeta" = OutputWAbeta
  )
  return(sali)
}

#----------------------------------- --
# Coverage freqs results for PL and APL intervals of quantile QGp:

CovResulQpRegion = function(HHH, MM) {
  # For Qp profile likelihood intervals and asymptotic 95% :
  CovQp = HHH$CovQGp
  AC_Qp = AgrestiCoul(CovQp, MM)
  MatQp = matrix(cbind(CovQp, AC_Qp), nrow = 2, ncol = 3)
  colnames(MatQp) = c("CovFreqs", "AC_Qp_Lower", "AC_Qp_Upper")
  rownames(MatQp) = c("Qp_PL", "Qp_APL")
  # For Likelihood Region exact LR and asymptotic AR of 95% confidence:
  CovLRAR = HHH$CovReg
  AC_LRAR = AgrestiCoul(CovLRAR, MM)
  MatReg = matrix(cbind(CovLRAR, AC_LRAR), nrow = 2, ncol = 3)
  colnames(MatReg) = c("CovFreqs", "AC_Reg_Lower", "AC_Reg_Upper")
  rownames(MatReg) = c("ProposedRegion", "AsymptoticRegion")
  sali = list("OutputQp" = MatQp, "OutputRegion" = MatReg)
  return(sali)
}

# ______________________________________####
# Simulations for a NEW setting  ####
## and a SINGLE sample size ####
# Computing coverage probabilities for a NEW valid single setting 
# of true parameters (A0, BETA0,N,P), (as defined in 
#   the readme for users file), by simulating M Gumbel samples.
#----------------------------------- --

# The next function estimates the coverage probabilities of 95% 
# confidence intervals PL,APL,WA for M simulated Gumbel samples,
# of sample size N with their corresopnding Agresti-Coull's intervals.

Coverprobintervals = function(A0 = 2, BETA0 = 3, N = 20, P = 0.25) {
  Z = 5 # Control value
  if ((P < 1 / (N + 1)) || (P > (N / (N + 1)))) {
    print("Error: Probability P must be a value within [1/(N+1),N/(N+1)]")
    Z = 2
  }
  if (N < 10) {
    print("Sample size N must be larger or equal to 10")
    Z = 2
  }
  if ((BETA0 < 0.3) || (BETA0 > 30)) {
    print("Weibull Shape parameter beta must be within [0.3,30]")
    Z = 2
  }
  if ((A0 < -3) || (A0 > 10)) {
    print("Gumbel location parameter A must be within [-3,10]")
    Z = 2
  }
  if (Z > 2) {
    M = 10000 # Recommended number of simulated Gumbel samples
    HH = CoverFreq(A0, BETA0, N, P, M)
    cat(paste("For sample size",N))
    cat("\n")
    cat(paste("true parameters A =",A0,"BETA =",BETA0))
    cat("\n")
    cat(paste("and probability of the quantile of interest, P=",P))
    cat("\n")
    cat("the coverage freqs. and the Agresti Coull intervals of the")
    cat("\n")
    cat("corresponding coverage probabilities were:")
    cat("\n")
    print(CovResulA(HH, M))
    print(CovResulBeta(HH, M))
    print(CovResulQpRegion(HH, M))
  } else {
    print("Please correct the indicated value and run again")
  }
} #END of function Coverprobintervals

# ______________________________________####
# Simulations for a NEW setting of Weibull parameters  ####  
## & SEVERAL sample sizes  ####
#  (same 19 considered in the article).
# The coverage frequencies of all considered intervals,
# PL, APL, and WA intervals are computed for M=10,000 
# simulated samples of 19 different sizes n={10,15,20,...,100}.
# Coverage freqs. of profile likelihood intervals for a quantile of interest
# are computed too as well as for the proposed likelihood regions.
# for simultaneous estimation of both Weibull parameters. 

# Results are exported to the output text file "SimulationResults.txt"
# within the folder WeibullSimulationResults generated in the working
# directory. Similar plots as those shown in Figures 1 and 2 of the article
# are exported to the same folder.
WeibullSimulations = function(a0 = 2, beta0 = 3, p = 0.25) {
  Z = 5 # Control values
  if ((p < 0.09) || (p > 0.91)) {
    print("Error: Probability P must be a value within [0.09,0.91]")
    Z = 2
  }
  if ((beta0 < 0.3) || (beta0 > 30)) {
    print("Weibull Shape parameter beta must be within [0.3,30]")
    Z = 2
  }
  if ((a0 < -3) || (a0 > 10)) {
    print("Gumbel location parameter A must be within [-3,10]")
    Z = 2
  }
  if (Z > 2) {
    MM = 10000 # Recommended number of simulated Gumbel samples
    # The true value of the quantile of interest is:
    Qp0 = a0 + (log(-log(1 - p)) / beta0)
    # The following are the 19 sample sizes considered for
    #  Figures 1 and 2 in the article.
    SampleSizes = seq(10, 100, by = 5) 
    nk = length(SampleSizes)
    Nss = SampleSizes
    # Creating matrices that will contain coverage freqs results for all considered
    # intervals and sample sizes for the given setting of Gumbel parameters.
    # The following matrices have 3 columns, each contains the number of intervals
    # that fell under, covered, and over the true value of the parameter.
    # The nk rows correspond each to a given sample size of the nk considered.
    # There is a matrix for each type of interval considered: PL, APL, and WA.
    # For parameter a:
    CovmatPLa = matrix(0, nrow = nk, ncol = 3)
    CovmatAPLa = matrix(0, nrow = nk, ncol = 3)
    CovmatWAa = matrix(0, nrow = nk, ncol = 3)
    # For parameter Beta:
    CovmatPLBeta = matrix(0, nrow = nk, ncol = 3)
    CovmatAPLBeta = matrix(0, nrow = nk, ncol = 3)
    CovmatWABeta = matrix(0, nrow = nk, ncol = 3)
    # For quantile Qp and Likelihood regions for (a,beta):
    # only cov freq of PL & APL intervals are kept.
    CovmatQp = matrix(0, nrow = nk, ncol = 2) 
    # Cov freq of proposed and asymptotic likelihood regions are kept, per column.
    CovmatReg = matrix(0, nrow = nk, ncol = 2) 
    # Matrices for the corresponding AGresti-Coull's intervals of each cov. freq:
    ACmatPLa = matrix(0, nrow = nk, ncol = 6)
    ACmatAPLa = matrix(0, nrow = nk, ncol = 6)
    ACmatWAa = matrix(0, nrow = nk, ncol = 6)
    ACmatPLBeta = matrix(0, nrow = nk, ncol = 6)
    ACmatAPLBeta = matrix(0, nrow = nk, ncol = 6)
    ACmatWABeta = matrix(0, nrow = nk, ncol = 6)
    ACmatQp = matrix(0, nrow = nk, ncol = 4)
    ACmatLR = matrix(0, nrow = nk, ncol = 4)
   
    for (i in (1:nk)) {
      n = SampleSizes[i]
      print(paste("Simulating M=10,000 samples of size:", n))
      # Computing the coverage freqs of intervals for given setting and n:
      HHH = CoverFreq(A0 = a0, BETA0 = beta0, N = n, P = p, M = MM)
      
      # PL For a:
      CovPLa = HHH$AcovPL # this is a vector of length 3
      AC_PLa = AgrestiCoul(CovPLa, MM) # this is a matrix of 3x2
      CovmatPLa[i, ] = CovPLa # The i-th row having 3 columns is filled
      ACmatPLa[i, 1:2] = AC_PLa[1, ]
      ACmatPLa[i, 3:4] = AC_PLa[2, ]
      ACmatPLa[i, 5:6] = AC_PLa[3, ]
      # APL For a:
      CovAPLa = HHH$AcovAPL
      AC_APLa = AgrestiCoul(CovAPLa, MM)
      CovmatAPLa[i, ] = CovAPLa # The i-th having 3 columns is filled
      ACmatAPLa[i, 1:2] = AC_APLa[1, ]
      ACmatAPLa[i, 3:4] = AC_APLa[2, ]
      ACmatAPLa[i, 5:6] = AC_APLa[3, ]
      # WA For a:
      CovWAa = HHH$AcovWA
      AC_WAa = AgrestiCoul(CovWAa, MM)
      CovmatWAa[i, ] = CovAPLa # The i-th row having 3 columns is filled
      ACmatWAa[i, 1:2] = AC_WAa[1, ]
      ACmatWAa[i, 3:4] = AC_WAa[2, ]
      ACmatWAa[i, 5:6] = AC_WAa[3, ]
      #------------------------------- --
      # PL For Beta:
      CovPLBeta = HHH$BetacovPL # this is a vector of length 3
      AC_PLBeta = AgrestiCoul(CovPLBeta, MM) # this is a matrix of 3x2
      CovmatPLBeta[i, ] = CovPLBeta # The i-th row having 3 columns is filled
      ACmatPLBeta[i, 1:2] = AC_PLBeta[1, ]
      ACmatPLBeta[i, 3:4] = AC_PLBeta[2, ]
      ACmatPLBeta[i, 5:6] = AC_PLBeta[3, ]
      # APL For Beta:
      CovAPLBeta = HHH$BetacovAPL
      AC_APLBeta = AgrestiCoul(CovAPLBeta, MM)
      CovmatAPLBeta[i, ] = CovAPLBeta # The i-th having 3 columns is filled
      ACmatAPLBeta[i, 1:2] = AC_APLBeta[1, ]
      ACmatAPLBeta[i, 3:4] = AC_APLBeta[2, ]
      ACmatAPLBeta[i, 5:6] = AC_APLBeta[3, ]
      # WA For Beta:
      CovWABeta = HHH$BetacovWA
      AC_WABeta = AgrestiCoul(CovWABeta, MM)
      CovmatWABeta[i, ] = CovAPLBeta # The i-th row having 3 columns is filled
      ACmatWABeta[i, 1:2] = AC_WABeta[1, ]
      ACmatWABeta[i, 3:4] = AC_WABeta[2, ]
      ACmatWABeta[i, 5:6] = AC_WABeta[3, ]
      #------------------------------- --
      # For Qp profile likelihood intervals and asymptotic 95% :
      CovQp = HHH$CovQGp # This is a vector of length 2.
      AC_Qp = AgrestiCoul(CovQp, MM) # This is a matrix of 2x2.
      CovmatQp[i, ] = CovQp
      ACmatQp[i, 1:2] = AC_Qp[1, ]
      ACmatQp[i, 3:4] = AC_Qp[2, ]
      
      # For proposed and asymptotic ikelihood regions of 95% confidence:
      CovLR = HHH$CovReg # This is a vector of length 2.
      CovmatReg[i, ] = CovLR
      AC_R = AgrestiCoul(CovLR, MM) # This is a matrix of 2x2.
      ACmatLR[i, 1:2] = AC_R[1, ]
      ACmatLR[i, 3:4] = AC_R[2, ]
    }
    #-------------------------- --
    # Generating Table 1 for new settings:
    Z = cbind(ACmatQp, ACmatLR)
    Z1 = round(Z, digits = 4)
    T1 = cbind(SampleSizes, Z1)
    T1nomis = c(
      "n", "LRCovLower", "LRCovUpper", "ARCovLower", "ARCovUpper",
      "PLQpCovLower", "PLQpCovUpper", "APLQpCovLower", "APLQpCovUpper"
    )
    colnames(T1) = T1nomis
    newTable1 = data.frame(T1, row.names = NULL)
    #-------------------------- --
    # Generating Table 3 for new settings
    Z = cbind(ACmatPLa, ACmatPLBeta)
    Z3 = round(Z, digits = 4)
    T3 = cbind(SampleSizes, Z3)
    T3nomis = c(
      "n", "PLaUndLo", "PLaUndUp", "PLaCovLo", "PLaCovUp", "PLaOverLo",
      "PLaOverUp", "PLBetUndLo", "PLBetUndUp", "PLBetCovLo", "PLBetCovUp",
      "PLBetOverLo", "PLBetOverUp"
    )
    colnames(T3) = T3nomis
    newTable3 = data.frame(T3, row.names = NULL)
    #------------------------------------- --
    # Generating Table 4 for new settings:
    Z = cbind(ACmatAPLa, ACmatAPLBeta)
    Z4 = round(Z, digits = 4)
    T4 = cbind(SampleSizes, Z4)
    T4nomis = c(
      "n", "APLaUndLo", "APLaUndUp", "APLaCovLo", "APLaCovUp", "APLaOverLo",
      "APLaOverUp", "APLBetUndLo", "APLBetUndUp", "APLBetCovLo", "APLBetCovUp",
      "APLBetOverLo", "APLBetOverUp"
    )
    colnames(T4) = T4nomis
    newTable4 = data.frame(T4, row.names = NULL)
    #------------------------------------- --
    # Generating Table 5 for new settings: 
    Z = cbind(ACmatWAa, ACmatWABeta)
    Z5 = round(Z, digits = 4)
    T5 = cbind(SampleSizes, Z5)
    T5nomis = c(
      "n", "WAaUndLo", "WAaUndUp", "WAaCovLo", "WAaCovUp", "WAaOverLo",
      "WAaOverUp", "WABetUndLo", "WABetUndUp", "WABetCovLo", "WABetCovUp",
      "WABetOverLo", "WABetOverUp"
    )
    colnames(T5) = T5nomis
    newTable5 = data.frame(T5, row.names = NULL)
    #------------------------------------- --
    # Generating Output file:  WeibullSimulationResults.txt"  
    dir.create("./WeibullSimulationResults", showWarnings = FALSE)
    sink(file = "WeibullSimulationResults/SimulationResults.txt")
    cat("The following number of Gumbel samples were simulated:")
    cat("\n")
    cat(MM)
    cat("\n")
    cat("for each of the following sample sizes, n:")
    cat("\n")
    cat(SampleSizes)
    cat("\n")
    cat("having the following  parameters (as their true values):")
    cat("\n")
    cat("The Gumbel location parameter a was:")
    cat("\n")
    cat(a0)
    cat("\n")
    cat("where the Weibull scale parameter sigma=exp(a) is:")
    cat("\n")
    cat(round(exp(a0),digits=2))
    cat("\n")
    cat("The Weibull shape parameter beta was:")
    cat("\n")
    cat(beta0)
    cat("\n")
    cat("The corresponding Gumbel scale parameter b=1/beta is:")
    cat("\n")
    cat(round((1/beta0),digits=2))
    cat("\n")
    cat("For probability:")
    cat("\n")
    cat(p)
    cat("\n")
    cat("the Gumbel quantile QGp was:")
    cat("\n")
    cat(Qp0)
    cat("\n")
    cat("and the corresponding Weibull quantile QWp was:")
    cat("\n")
    cat(exp(Qp0))
    cat("\n")
    cat("________________________________")
    cat("\n")
    cat("\n")
    cat("The new Table 1 with the Agresti-Coull intervals corresponding to")
    cat("\n")
    cat(" the Coverage Probabilities for the new M simulated samples")
    cat("\n")
    cat("of the proposed (LR) and asymptotic (AR) likelihood regions for (a,beta)")
    cat("\n")
    cat("and of the proposed and asymptotic likelihood intervals for quantile Qp")
    cat("\n")
    cat("is given by:")
    cat("\n")
    print(newTable1[, 1:7])
    cat("\n")
    print(newTable1[, 8:9])
    cat("\n")
    cat("________________________________")
    cat("\n")
    cat("For the provided setting, and new M simulated Gumbel samples, ")
    cat("\n")
    cat("the following tables are equivalent to Tables 3 to 5.")
    cat("\n")
    cat("The new Table 3 with the Agresti-Coull intervals for Coverage Frequencies")
    cat("\n")
    cat("for PL intervals for a & beta is given by:")
    cat("\n")
    print(newTable3[, 1:7])
    cat("\n")
    print(newTable3[, 8:13])
    cat("\n")
    cat("-----------------------------")
    cat("\n")
    cat("The new Table 4 with the Agresti-Coull intervals for Coverage Frequencies")
    cat("\n")
    cat("for APL intervals for a & beta is given by:")
    cat("\n")
    print(newTable4[, 1:7])
    cat("\n")
    print(newTable4[, 8:13])
    cat("\n")
    cat("-----------------------------")
    cat("\n")
    cat("The new Table 5 with the Agresti-Coull intervals for Coverage Frequencies")
    cat("\n")
    cat("for WA intervals for a & beta is given by:")
    cat("\n")
    print(newTable5[, 1:7])
    cat("\n")
    print(newTable5[, 8:13])
    cat("\n")
    cat("________________________________")
    cat("\n")
    sink(file = NULL)
    #------------------------------------- --
    # Generating Figure 1 for the NEW simulation settings  
    # Figure 1. Agresti-Coull's 95% confidence intervals for the coverage
    # probabilities of 95% confidence intervals for a, for selected sample sizes n
    # and selected simulation setting with M=10,000 samples. PL intervals in solid
    # line and within dots, APL intervals in dashes within asterisks, and WA
    # intervals in dashes within crosses.
    png("WeibullSimulationResults/NewFigure1.png", width = 600, height = 350)
    plot(c(5, 105), c(0.95, 0.95),
         type = "l", lty = 1,
         ylab = "Coverage Prob. (a)",
         xlab = "Sample Size n", xlim = c(5, 105), ylim = c(0.85, 0.962), yaxs = "i",
         xaxs = "i",
         main = paste("Figure 1. Intervals for a for true values a=",a0,",beta=",
                      beta0,",p=",p)
    )
    axis(side = 1, at = nk, tick = TRUE)
    for (i in (1:nk)) {
      lines(c(Nss[i], Nss[i]), T3[i, 4:5], lty = 1, col = 1)
      points(c(Nss[i], Nss[i]), T3[i, 4:5], pch = 19, col = 1)
      points(c(Nss[i], Nss[i]), T4[i, 4:5], pch = 8, col = 1)
      lines(c(Nss[i], Nss[i]), T4[i, 4:5], lty = 2, col = 1)
      points(c(Nss[i], Nss[i]), T5[i, 4:5], pch = 4, col = 1)
      lines(c(Nss[i], Nss[i]), T5[i, 4:5], lty = 2, col = 1)
    }
    lines(c(5, 105), c(0.96, 0.96), lty = 2, col = 1)
    lines(c(5, 105), c(0.94, 0.94), lty = 2, col = 1)
    dev.off()
    #------------------------------------- --
    # Generating Figure 2 for NEW simulation settings 
    # Figure 2. Agresti-Coull's 95% confidence intervals for the coverage
    # probabilities of 95% confidence intervals for Beta, for selected 
    # sample sizes n and simulation setting with M=10,000 samples. 
    # PL intervals in solid line and within dots, 
    # APL intervals appear in dashes within asterisks, and 
    # WA intervals in dashes within crosses.
    png("WeibullSimulationResults/NewFigure2.png", width = 600, height = 350)
    plot(c(5, 105), c(0.95, 0.95),
         type = "l", lty = 1,
         ylab = "Coverage Prob. (a)",
         xlab = "Sample Size n", xlim = c(5, 105), ylim = c(0.89, 0.962), yaxs = "i",
         xaxs = "i", 
         main = paste("Figure 2. Intervals for beta for true values a=",a0,",beta=",
                      beta0,",p=",p)
    )
    axis(side = 1, at = nk, tick = TRUE)
    for (i in (1:nk)) {
      lines(c(Nss[i], Nss[i]), T3[i, 10:11], lty = 1, col = 1)
      points(c(Nss[i], Nss[i]), T3[i, 10:11], pch = 19, col = 1)
      points(c(Nss[i], Nss[i]), T4[i, 10:11], pch = 8, col = 1)
      lines(c(Nss[i], Nss[i]), T4[i, 10:11], lty = 2, col = 1)
      points(c(Nss[i], Nss[i]), T5[i, 10:11], pch = 4, col = 1)
      lines(c(Nss[i], Nss[i]), T5[i, 10:11], lty = 2, col = 1)
    }
    lines(c(5, 105), c(0.96, 0.96), lty = 2, col = 1)
    lines(c(5, 105), c(0.94, 0.94), lty = 2, col = 1)
    dev.off()
  } else {
    print("Please correct the indicated value and run again")
  } # End of conditional clause: "if (Z>2)..."
}  # Ebd if function WeibullSimulations


