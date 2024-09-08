# Generating Figure 5 and Example 4.3 results  ####
# Tech23_068
# _______________________________________________ ####

# #  WARNING: SOURCE BEFORE: WeibullFunctions.R  ####
# #  (It contains functions tO be used.) ####

# _______________________________________________ ####
# ##  ........................................................
# #   Weibull Data: Example 4.3 Lawless p. 240  ####
# #  Comparison of two types of electrical cable insulation
## Type 1:
Y1 = c(
  32, 35.4, 36.2, 39.8, 41.2, 43.3, 45.5, 46, 46.2, 46.4, 46.5, 46.8, 47.3,
  47.3, 47.6, 49.2, 50.4, 50.9, 52.4, 56.3
)
Y = Y1

Y = sort(Y) # ordered Weibull data
X = log(Y) # ordered Gumbel logdata
n = length(Y)
t = sum(X)
#-------------------------- --
## Type 2:
Y2 = c(
  39.4, 45.3, 49.2, 49.4, 51.3, 52, 53.2, 53.2, 54.9, 55.5, 57.1, 57.2, 57.5,
  59.2, 61, 62.4, 63.8, 64.3, 67.3, 67.7
)
Y2 = sort(Y2) # ordered Weibull data
X2 = log(Y2) # ordered Gumbel logdata
n2 = length(Y2)
t2 = sum(X2)

## ------------------------------ --
## Type 1: ####
# Computing Weibull mles:  ####

B = Weibullmles(Y)
# Mle of a, Gumbel location parameter:
amle = B$AGumbel
# Mle of beta, Weibull shape parameter:
betamle = B$BetaShape
# Mle of sigma, WEibull scale parameter:
sigmamle = B$SigmaScale
# Mle of b, Gumbel scale parameter:
bmle = 1 / betamle
# Maximum value of the loglikelihood at the mle:
LVmax = lvabetaWei(c(amle, betamle), X, t, n)

#------------------------------ --
# Computing Fisher's Observed Informations:  ####

B = InformationABeta(Y, amle, betamle)
Infoa = B$Inf_a
Infobeta = B$Inf_beta
Infoabeta = B$Infabeta

# Quantiles mles and intervals: ####
#  Mles of quantiles 0.10 and 0.90
P = 0.10
Wp = log(-log(1 - P))
Qg10mle = amle + Wp / betamle
Qw10mle = exp(Qg10mle)

P = 0.20
Wp = log(-log(1 - P))
Qg20mle = amle + Wp / betamle
Qw20mle = exp(Qg20mle)


P = 0.90
Wp = log(-log(1 - P))
Qg90mle = amle + Wp / betamle
Qw90mle = exp(Qg90mle)



#-------------------------- --
# Quantiles QGp and QWp mles and intervals:
# Profile Likelihood intervals for quantile p=0.10:
IntQW10 = IntervalQP(0.10, Y, n, amle, betamle)
# Profile Likelihood intervals for quantile p=0.20:
IntQW20 = IntervalQP(0.20, Y, n, amle, betamle)
# Profile Likelihood intervals for quantile p=0.90:
IntQW90 = IntervalQP(0.90, Y, n, amle, betamle)


## ------------------------------ --
## Type 2: ####
# Computing Weibull mles:  ####
B2 = Weibullmles(Y2)
# Mle of a, Gumbel location parameter:
amle2 = B2$AGumbel
# Mle of beta, Weibull shape parameter:
betamle2 = B2$BetaShape
# Mle of sigma, WEibull scale parameter:
sigmamle2 = B2$SigmaScale
# Mle of b, Gumbel scale parameter:
bmle2 = 1 / betamle2
# Maximum value of the loglikelihood at the mle:
LVmax2 = lvabetaWei(c(amle2, betamle2), X2, t2, n2)

#------------------------------ --
# Computing Fisher's Observed Informations:  ####
B2 = InformationABeta(Y2, amle2, betamle2)
Infoa2 = B2$Inf_a
Infobeta2 = B2$Inf_beta
Infoabeta2 = B2$Infabeta

# Quantiles mles and intervals: ####
#  Mles of quantiles 0.10 and 0.90
P = 0.10
Wp = log(-log(1 - P))
Qg10mle2 = amle2 + Wp / betamle2
Qw10mle2 = exp(Qg10mle2)

P = 0.20
Wp = log(-log(1 - P))
Qg20mle2 = amle2 + Wp / betamle2
Qw20mle2 = exp(Qg20mle2)

P = 0.90
Wp = log(-log(1 - P))
Qg90mle2 = amle2 + Wp / betamle2
Qw90mle2 = exp(Qg90mle2)


# Type 2: Quantiles QGp and QWp mles and intervals:
IntQW10_2 = IntervalQP(0.10, Y2, n2, amle2, betamle2)
# Profile Likelihood intervals for quantile p=0.20:
IntQW20_2 = IntervalQP(0.20, Y2, n2, amle2, betamle2)
# Profile Likelihood intervals for quantile p=0.90:
IntQW90_2 = IntervalQP(0.90, Y2, n2, amle2, betamle2)


# ________________________________________ ####
# Making plots of Rp(beta) & Rp*(beta) for both Types: ####

# Profile likelihood of beta Rp(beta) for both types.
Nbeta = 250
betis = seq(4, 16, length.out = Nbeta)
# For Type 1:
Rpbeta1 = rep(0, Nbeta)
# For Type 2:
Rpbeta2 = rep(0, Nbeta)
for (i in (1:Nbeta)) {
  BETA = betis[i]
  # Rp(beta) for Type 1:
  Rpbeta1[i] = RpbetaWei(BETA, X, t, n, betamle)
  # Rp(beta) for Type 2:
  Rpbeta2[i] = RpbetaWei(BETA, X2, t2, n2, betamle2)
}

# Rp*(beta) Approximating functions for both types:
RpbetaSt1 = RpbetaStWei(betis, betamle, Infobeta)
RpbetaSt2 = RpbetaStWei(betis, betamle2, Infobeta2)

# PL 95% confidence intervals for beta for both types are obtained:
# For Type 1:
Z1 = IntervalsBeta(n, betamle, Infobeta)$BetaIntervals
PLintBeta = Z1
# For Type 2:
Z2 = IntervalsBeta(n2, betamle2, Infobeta2)$BetaIntervals
PLintBeta2 = Z2
CS = clikelevels(n)$csig

# Figure 5a R*(beta) for two types:  ####
#-------------------------------------- --
# Plot of Rp*(beta) with 95% confidence interval for Types 1 and 2:
# this is Figure 5a:
X11()
plot(betis, RpbetaSt1,
  type = "l", ylim = c(0, 1),
  main = "(5a)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 1, lty = 1
)
# lines(c(1,1),c(0,1),lty=2,col=2)
lines(betis, RpbetaSt2, type = "l", lty = 2, col = 1, lwd = 2)
# PL intervals of 95% confidence:
# Type 1:
lines(Z1[2, ], c(CS[2], CS[2]), lty = 1, col = 1, lwd = 2)
lines(Z2[2, ], c(CS[2], CS[2]), lty = 2, col = 1, lwd = 2)
#-------------------------------------- --
#------------------------------- --
# Making plots of Rp*(sigma) for both Types: ####
#-------------------------------------- --
# For Type 1
B = IntervalsA(n, amle, betamle, Infoa)
Aendpt = B$EndpointsA
SigmaInt = B$SigmaIntervals
AInt = B$AIntervals
AA1 = Aendpt[1]
AA2 = Aendpt[2]
cs = B$LikeLevels
c90 = cs[1]
c95 = cs[2]
c99 = cs[3]
# Setting range for plot of Rp(a) and Rp*(a) that includes amle:
nhalf = 100
As1 = seq(Aendpt[1], amle, length.out = nhalf)
difi = As1[2] - As1[1]
As2 = seq((amle + difi), Aendpt[2], length.out = nhalf)
Ais1 = c(As1, As2)
nk = length(Ais1)
#--- Calculating Rp(a) and Rp*(a) over plot range here:
RpaExact = rep(0, nk)
BetarmleA = rep(0, nk)
RpaSt = rep(0, nk)
LVmax1 = lvabetaWei(c(amle, betamle), X, t, n)

for (i in 1:nk) {
  AA = Ais1[i]
  B = RpA(AA, X, t, n, LVmax1, betamle)
  RpaExact[i] = B$Rpas
  BetarmleA[i] = B$Betarmle
  RpaSt[i] = RpStarA(AA, X, n, amle, betamle, Infoa)
}

#------------------------------ --
# For Type 2:
B2 = IntervalsA(n2, amle2, betamle2, Infoa2)
Aendpt2 = B2$EndpointsA
SigmaInt2 = B2$SigmaIntervals
AInt2 = B2$AIntervals
AA1 = Aendpt2[1]
AA2 = Aendpt2[2]
cs2 = B2$LikeLevels
c90 = cs2[1]
c95 = cs2[2]
c99 = cs2[3]
# Setting range for plot of Rp(a) and Rp*(a) that includes amle:
nhalf = 100
As1 = seq(Aendpt2[1], amle2, length.out = nhalf)
difi = As1[2] - As1[1]
As2 = seq((amle2 + difi), Aendpt2[2], length.out = nhalf)
Ais2 = c(As1, As2)
nk = length(Ais2)
#--- Calculating Rp(a) and Rp*(a) over plot range here:
RpaExact2 = rep(0, nk)
BetarmleA2 = rep(0, nk)
RpaSt2 = rep(0, nk)
LVmax2 = lvabetaWei(c(amle2, betamle2), X2, t2, n2)

for (i in 1:nk) {
  AA = Ais2[i]
  B2 = RpA(AA, X2, t2, n2, LVmax2, betamle2)
  RpaExact2[i] = B2$Rpas
  BetarmleA2[i] = B2$Betarmle
  RpaSt2[i] = RpStarA(AA, X2, n2, amle2, betamle2, Infoa2)
}

# Figure 5b R*(sigma) for two types:  ####
#----------------------- --
# Plot of Rp(sigma), Rp*(sigma) and proposed intervals:
X11()
# Type 1: R*(a) in solid black line:
plot(exp(Ais1), RpaSt,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = exp(c(Aendpt[1], Aendpt2[2])),
  main = "(5b)", xlab = expression(sigma),
  ylab = expression(paste("R(", sigma, ")")),
  yaxs = "i"
)
# Type 2: R*(a) in dashes:
lines(exp(Ais2), RpaSt2, lty = 2, col = 1, lwd = 2)
#----- 95% confidence interval for Type 1:
lines(exp(AInt[2, ]), rep(c95, 2), lty = 1, col = 1, lwd = 2)
#----- 95% confidence interval for Type 2:
lines(exp(AInt2[2, ]), rep(c95, 2), lty = 1, col = 1, lwd = 2)


# _____________________ ####
# Figures 5 a & b together: ####
X11()
split.screen(c(1, 2))
screen(1)
plot(betis, RpbetaSt1,
  type = "l", ylim = c(0, 1),
  main = "(5a)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 1, lty = 1
)
# lines(c(1,1),c(0,1),lty=2,col=2)
lines(betis, RpbetaSt2, type = "l", lty = 2, col = 1, lwd = 2)
# PL intervals of 95% confidence:
# Type 1:
lines(Z1[2, ], c(CS[2], CS[2]), lty = 1, col = 1, lwd = 2)
lines(Z2[2, ], c(CS[2], CS[2]), lty = 2, col = 1, lwd = 2)
# PL intervals of 95% confidence:
# Type 1:
lines(Z1[2, ], c(CS[2], CS[2]), lty = 1, col = 1)
lines(Z2[2, ], c(CS[2], CS[2]), lty = 2, col = 1)
#----------------- --
screen(2)
plot(exp(Ais1), RpaSt,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = exp(c(Aendpt[1], Aendpt2[2])),
  main = "(5b)", xlab = expression(sigma),
  ylab = expression(paste("R(", sigma, ")")),
  yaxs = "i"
)
# Type 2: R*(a) in dashes:
lines(exp(Ais2), RpaSt2, lty = 2, col = 1, lwd = 3)
#----- 95% confidence interval for Type 1:
lines(exp(AInt[2, ]), rep(c95, 2), lty = 1, col = 1, lwd = 2)
#----- 95% confidence interval for Type 2:
lines(exp(AInt2[2, ]), rep(c95, 2), lty = 1, col = 1, lwd = 2)

# ______________________________________ ####
#       Generating outpuf file with results:                 ####
sink(file = "Example4_3Output.txt")
cat("Example 4.3. Comparison of Two Treatments and Figure 5 (a and b).")
cat("\n")
cat("____________________________________________")
cat("\n")
cat("\n")
cat("TYPE 1: Sample size was:")
cat(n)
cat("\n")
cat("The mles of a, beta, sigma, and b are:")
cat("\n")
cat(amle)
cat("\n")
cat(betamle)
cat("\n")
cat(sigmamle)
cat("\n")
cat(bmle)
cat("\n")
cat("The maximum value of the loglikelihood (a,beta)at the mle was:")
cat("\n")
cat(LVmax)
cat("\n")
cat("\n")
cat("------------------------------------")
cat("\n")
cat("The observed Fisher's informations for a and beta are:")
cat("\n")
cat(Infoa)
cat("\n")
cat(Infobeta)
cat("\n")
cat("Fisher's observed information matrix for (a,beta) is:")
cat("\n")
print(Infoabeta)
cat("\n")
cat("------------------------------------")
cat("\n")
cat("\n")
cat("Type 1: PL intervals for beta:")
cat("\n")
print(PLintBeta)
cat("\n")
cat("Type 1: PL intervals for sigma:")
cat("\n")
print(SigmaInt)
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Type 1: Mles of Gumbel Quantiles 10, 20 and 90 are:")
cat("\n")
cat(paste(Qg10mle, ",", Qg20mle, ",", Qg90mle))
cat("\n")
cat("Type 1: Mles of Weibull Quantiles 10, 20 and 90 are:")
cat("\n")
cat(paste(Qw10mle, ",", Qw20mle, ",", Qw90mle))
cat("\n")
cat("------------------------------------")
cat("\n")
cat("\n")
cat("Type 1: confidence intervals for WEibull Q10:")
cat("\n")
print(IntQW10$WeibullQPIntervals)
cat("\n")
cat("Type 1: confidence intervals for Weibull Q20:")
cat("\n")
print(IntQW20$WeibullQPIntervals)
cat("\n")
cat("Type 1: confidence intervals for Weibull Q90:")
cat("\n")
print(IntQW90$WeibullQPIntervals)
cat("\n")
cat("Likelihood levels of intervals of a single parameter for n:")
cat("\n")
print(IntQW90$LikeLevels)
cat("\n")
cat("------------------------------------")
cat("\n")
cat("TYPE 2: Sample size was:")
cat(n2)
cat("\n")
cat("TYPE 2: The mles of a, beta, sigma, and b are:")
cat("\n")
cat(amle2)
cat("\n")
cat(betamle2)
cat("\n")
cat(sigmamle2)
cat("\n")
cat(bmle2)
cat("\n")
cat("The maximum value of the loglikelihood (a,beta)at the mle was:")
cat("\n")
cat(LVmax2)
cat("\n")
cat("\n")
cat("------------------------------------")
cat("\n")
cat("The observed Fisher's informations for a and beta are:")
cat("\n")
cat(Infoa2)
cat("\n")
cat(Infobeta2)
cat("\n")
cat("Fisher's observed information matrix for (a,beta) is:")
cat("\n")
print(Infoabeta2)
cat("\n")
cat("------------------------------------")
cat("\n")
cat("\n")
cat("Type 2: PL intervals for beta:")
cat("\n")
print(PLintBeta2)
cat("\n")
cat("Type 2: PL intervals for sigma:")
cat("\n")
print(SigmaInt2)
cat("\n")
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Type 2: Mles of Gumbel Quantiles 10, 20 and 90 are:")
cat("\n")
cat(paste(Qg10mle2, ",", Qg20mle2, ",", Qg90mle2))
cat("\n")
cat("Type 2: Mles of Weibull Quantiles 10, 20 and 90 are:")
cat("\n")
cat(paste(Qw10mle2, ",", Qw20mle2, ",", Qw90mle2))
cat("\n")
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Type 2: confidence intervals for Weibull Q10:")
cat("\n")
print(IntQW10_2$WeibullQPIntervals)
cat("\n")
cat("Type 2: confidence intervals for Weibull Q20:")
cat("\n")
print(IntQW20_2$WeibullQPIntervals)
cat("\n")
cat("Type 2: confidence intervals for Weibull Q90:")
cat("\n")
print(IntQW90_2$WeibullQPIntervals)
cat("\n")
cat("Likelihood levels for n2:")
cat("\n")
print(IntQW90$LikeLevels)
sink(file = NULL)
