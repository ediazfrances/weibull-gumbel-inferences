# Generating Example 4.1 results and Figure 3 ####
# Tech23_068 ####
# _______________________________________________ ####

# #  IMPORTANT: SOURCE FIRST: WeibullFunctions.R  ####
# #  (It contains functions to be used .) ####

# _______________________________________________ ####
# ##  ............................................................................
# #   Example 4.1 Failure times of Air Conditioning Equipment  ####
# ##   in hours of operation. Proschan (1963, Airplane 7911)
Y = c(55, 320, 56, 104, 220, 239, 47, 246, 176, 182, 33, 15, 104, 35)
Y = sort(Y) # ordered Weibull data
X = log(Y) # ordered Gumbel logdata
n = length(Y)
t = sum(X)
## ------------------------------ --
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


# Validation of the Weibull model.  ####
## Anderson Darling Modified Statistic (Stephens 1977) ####
ModAD = AndersonDarlingModif(Y, amle, betamle)
ModADQuantiles = c(0.637, 0.757, 0.877, 1.038)
names(ModADQuantiles) = c("Q90", "Q95", "Q975", "Q99")


#------------------------------- ---
# Figure 3a R*(beta) and R(beta):  ####
#------------------------------------ --
#  PL Intervals for beta
#  of likelihood levels c90, c95, c99 and c1=0.001
BB = IntervalsBeta(n, betamle, Infobeta)
Betaint = BB$BetaIntervals
Betaend = BB$EndpointsBeta
cs = BB$LikeLevels
c90 = cs[1]
c95 = cs[2]
c99 = cs[3]
Nbet = 150
betis = seq(Betaend[1], Betaend[2], length.out = Nbet)
# R*(beta) Approximating function fpr R(beta)
RpbetaSt = RpbetaStWei(betis, betamle, Infobeta)
# Exact relative profile likelihood of beta, R(beta),
# for comparison
Rpbetas = rep(0, Nbet)
for (i in (1:Nbet)) {
  BETA = betis[i]
  Rpbetas[i] = RpbetaWei(BETA, X, t, n, betamle)
}

#-------------------------------------- --
X11()
plot(betis, RpbetaSt,
  type = "l", ylim = c(0, 1), lty = 2,
  main = "(3a)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 2
)
lines(c(1, 1), c(0, 1), lty = 2, col = 1, lwd = 2)
lines(betis, Rpbetas, type = "l", lty = 1, col = 1, lwd = 2)
# PL intervals :
lines(Betaint[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

#------------------------------- --
# Figure 3b R(sigma), R*(sigmaa) and R(theta) where theta=Exponential mean: ####
#-------------------------------------- --
B = IntervalsA(n, amle, betamle, Infoa)
Aendpt = B$EndpointsA
SigmaInt = B$SigmaIntervals
AInt = B$AIntervals
cs = B$LikeLevels
c90 = cs[1]
c95 = cs[2]
c99 = cs[3]
# Setting range for plot of Rp(a) and Rp*(a) that includes amle:
nhalf = 100
As1 = seq(Aendpt[1], amle, length.out = nhalf)
difi = As1[2] - As1[1]
As2 = seq((amle + difi), Aendpt[2], length.out = nhalf)
Ais = c(As1, As2)
Sigis = exp(Ais)
nk = length(Ais)
#--- Calculating Rp(a) and Rp*(a) over plot range here:
RpaExact = rep(0, nk)
BetarmleA = rep(0, nk)
RpaSt = rep(0, nk)
LVmax = lvabetaWei(c(amle, betamle), X, t, n)

for (i in 1:nk) {
  AA = Ais[i]
  B = RpA(AA, X, t, n, LVmax, betamle)
  RpaExact[i] = B$Rpas
  BetarmleA[i] = B$Betarmle
  RpaSt[i] = RpStarA(AA, X, n, amle, betamle, Infoa)
}

# Relative likelihood of Exponential mean:
Expmle = mean(Y) # mean mle
Rtheta = RthetaExp(Sigis, Y)
ExpMedian = -Expmle * log(1 - 0.5)
ExpQ90 = -Expmle * log(1 - 0.9)
ExpQ10 = -Expmle * log(1 - 0.1)



# Plot of R(sigma), Rp*(sigma), PL intervals and R(theta) Exponential:
X11()
# R(sigmaa) in solid line, R*(sigma) in dashes, R(theta) exponent. in dashes:
plot(Sigis, RpaExact,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = exp(c(Aendpt[1], Aendpt[2])),
  main = "(3b)", xlab = expression(paste("QW63=", sigma)),
  ylab = expression(paste("R(", sigma, ")")),
  yaxs = "i"
)
# R*(a) in dashes:
lines(Sigis, RpaSt, lty = 2, col = 1, lwd = 2)
# R(theta) Exponential in dashes:
lines(Sigis, Rtheta, lty = 2, col = 1, lwd = 2)


#------------------------------- --
# Figure 3c R(Q10) for Weibull and Exponential Q10: ####
#-------------------------------------- --
# Weibull Quantiles, Mles :
#  Mles of quantiles 0.10 and 0.90
P = 0.10
Wp = log(-log(1 - P))
QG10mle = amle + Wp / betamle # Gumbel Q10 mle
QW10mle = exp(QG10mle) # Weibull Q10 mle
B = IntervalQP(P, Y, n, amle, betamle)
QGendpt = B$EndpointsQG
GumInt = B$GumbelQPIntervals
WeiInt = B$WeibullQPIntervals
cs = B$LikeLevels
c90 = cs[1]
c95 = cs[2]
c99 = cs[3]
# Setting range for plotting R(QgumP)
nk = 150
Qs = seq(QGendpt[1], QGendpt[2], length.out = (nk - 1))
Qs = sort(c(Qs, QG10mle)) # The quantile mle is included
Rq10s = rep(0, nk)
betars = rep(0, nk)
# Evaluating Rp(Qp) at points to be plotted:
for (i in (1:nk)) {
  QQ = Qs[i]
  Sali = RpQGp(QQ, P, X, t, n, LVmax, betamle)
  Rq10s[i] = Sali$RpQ
  betars[i] = Sali$betarmle
}

QWeis10 = exp(Qs)
# Calculating R(Qp) of exponential distribution:
RQ10exp = RQpExp(QWeis10, Y, P)

#  --------------------- --
X11() # Rp(QGp) in black solid line:
plot(QWeis10, Rq10s,
  type = "l", ylim = c(0, 1),
  main = "(3c)", lwd = 2,
  xlab = "QWp", yaxs = "i"
)
lines(QWeis10, RQ10exp, lty = 2, col = 1, lwd = 2)
#


#------------------------------- --
# Figure 3d R(Q90) for Weibull and Exponential Q90: ####
#-------------------------------------- --
P = 0.90
Wp = log(-log(1 - P))
Qg90mle = amle + Wp / betamle # Gumbel Q90 mle
Qw90mle = exp(Qg90mle) # Weibull Q90 mle
Q90Expmle = Expmle * (-log(1 - 0.90))

B = IntervalQP(P, Y, n, amle, betamle)
QGendpt = B$EndpointsQG
GumInt = B$GumbelQPIntervals
WeiInt = B$WeibullQPIntervals
cs = B$LikeLevels
c90 = cs[1]
c95 = cs[2]
c99 = cs[3]

# Setting range for plotting R(QgumP)
nk = 150
Qs = seq(QGendpt[1], QGendpt[2], length.out = (nk - 1))
Qs = seq(log(90), log(800), length.out = (nk - 1))
Qs = sort(c(Qs, QG10mle)) # The quantile mle is included
Rqs90 = rep(0, nk)
betars = rep(0, nk)
# Evaluating Rp(Qp) at points to be plotted:
for (i in (1:nk)) {
  QQ = Qs[i]
  Sali = RpQGp(QQ, P, X, t, n, LVmax, betamle)
  Rqs90[i] = Sali$RpQ
  betars[i] = Sali$betarmle
}

QWeis90 = exp(Qs)
QWeis90 = exp(Qs)
# Calculating R(Qp) of exponential distribution:
RQ90exp = RQpExp(QWeis90, Y, P)

#  --------------------- --
X11() # For Q90, Rp(QGp) in black solid line,
#       and of exponential in dashes:
plot(QWeis90, Rqs90,
  type = "l", ylim = c(0, 1),
  main = "(3d)", lwd = 2,
  xlab = "QW90, p=0.90", yaxs = "i"
)
lines(QWeis90, RQ90exp, lty = 2, col = 1, lwd = 2)
#


# __________________________ ####
# All Figures 3 a,b,c,d together: ####
X11()
split.screen(c(2, 2))
screen(1)
plot(betis, RpbetaSt,
  type = "l", ylim = c(0, 1), lty = 2,
  main = "(3a)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 1
)
lines(c(1, 1), c(0, 1), lty = 2, col = 1, lwd = 2)
lines(betis, Rpbetas, type = "l", lty = 1, col = 1, lwd = 2)
# PL intervals :
lines(Betaint[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

screen(2)
# R(sigma) in solid line:
plot(Sigis, RpaExact,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = exp(c(Aendpt[1], Aendpt[2])),
  main = "(3b)", xlab = expression(paste("QW63=", sigma, "p=0.63")),
  ylab = expression(paste("R(", sigma, ")")),
  yaxs = "i"
)
# R*(sigma) in dashes:
lines(Sigis, RpaSt, lty = 2, col = 1, lwd = 2)
# R(theta) Exponential in dashes:
lines(Sigis, Rtheta, lty = 2, col = 1, lwd = 2)

screen(3)
# Relative profile likelihoods of Weibull and Exponential distribs.
plot(QWeis10, Rq10s,
  type = "l", ylim = c(0, 1),
  main = "(3c)", lwd = 2,
  xlab = "QW10, p=0.10", yaxs = "i"
)
lines(QWeis10, RQ10exp, lty = 2, col = 1, lwd = 2)

screen(4)
plot(QWeis90, Rqs90,
  type = "l", ylim = c(0, 1),
  main = "(3d)", lwd = 2,
  xlab = "QW90, p=0.90", yaxs = "i"
)
lines(QWeis90, RQ90exp, lty = 2, col = 1, lwd = 2)


# ______________________________________ ####
#       Generating outpuf file with results:                 ####
sink(file = "Example4_1Output.txt")
cat("Example 4.1. Failure times of air conditioning equipment")
cat("\n")
cat("and Figure 3 (a,b,c,d).")
cat("\n")
cat("____________________________________________")
cat("\n")
cat("\n")
cat("Sample size was:")
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
cat("The Modified Anderson Darling (ModAD) Statistic for testing Weibull was:")
cat("\n")
print(ModAD)
cat("\n")
cat("to be compared with the following quantiles of its distribution")
cat("\n")
print(ModADQuantiles)
cat("\n")
cat("Overall, if ModAD<0.637, then the Weibull is reasonable for the data.")
cat("\n")
cat("\n")
cat("-------------------------------------")
cat("\n")
cat("\n")
cat(" PL intervals for beta:")
cat("\n")
print(Betaint)
cat("\n")
cat("Since all intervals include beta=1, ")
cat("\n")
cat("the exponential distribution is also reasonable.")
cat("\n")
cat("\n")
cat("\n")
cat("PL intervals for sigma:")
cat("\n")
print(SigmaInt)
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Mle of Exponential mean is:")
cat("\n")
print(Expmle)
cat("\n")
cat("This is also the quantile of p= 0.632")
cat("\n")
cat("Thus comparable with the same Weibull quantile, sigmamle:")
cat("\n")
print(sigmamle)
cat("\n")
cat("------------------------------------")
cat("\n")
cat(" Mles of Gumbel, Weibull, and Exponential Qp90 are:")
cat("\n")
cat(paste(Qg90mle, ",", Qw90mle, ",", Q90Expmle))
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Likelihood levels of intervals of a single parameter for n:")
cat("\n")
print(cs)
cat("\n")
cat("------------------------------------")
sink(file = NULL)
