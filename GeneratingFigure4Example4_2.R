# Generating Example 4.2 results and Figure 4 ####
# _______________________________________________ ####
# Tech23_068
# _______________________________________________ ####

# #  IMPORTANT: SOURCE FIRST: WeibullFunctions.R  ####
# #  (It contains functions tO be used.) ####

# _______________________________________________ ####
# ##  ............................................................................
# #   Example 4.2 Interarrival times to the Mummies Museum ####
Y = c(
  54, 24, 54, 6, 14, 87, 42, 30, 21, 4, 43, 34, 25, 26, 79, 35, 15, 25,
  33, 20, 12, 34, 24, 29, 28
)
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

#--------------- --
## Figure 4a. Weibull Probability plot ####
# Transformed ordered data with Weibull distribution
uwei = 1 - exp(-exp(betamle * (X - amle)))
enes = seq(1, n, by = 1)
aenes = n + 1 - enes
# Theoretical standard uniform quantiles:
Ebetas = enes / (n + 1)
## 95% approximate confidence band for the PP plot:
tau = 1 - (.05 / n)
tau1 = (1 - tau) / 2
tau2 = (1 + tau) / 2
Ban1 = qbeta(tau1, enes, aenes)
Ban2 = qbeta(tau2, enes, aenes)

X11()
plot(Ebetas, uwei,
  pch = 19, cex = .5, ylab = "F(x) Transform",
  xlab = "Uniform (0,1) Quantiles", xlim = c(0, 1), ylim = c(0, 1),
  main = "(4a)", pty = "m", col = 1
)
lines(c(0, 1), c(0, 1), lty = 2, col = 1)
lines(Ebetas, uwei, lty = 2, col = 1)
lines(Ebetas, Ban1, lty = 2, col = 1)
lines(Ebetas, Ban2, lty = 2, col = 1)

#------------------------------- --
# Figure 4b R(a) and R*(a): ####
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

# Plot of R(a), Rp*(a) and PL intervals:
X11()
# R(a) in solid line:
plot(Ais, RpaExact,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = c(Aendpt[1], Aendpt[2]),
  main = "(5b)", xlab = expression(a),
  ylab = expression(paste("R(", a, ")")),
  yaxs = "i"
)
# R*(a) in dashes:
lines(Ais, RpaSt, lty = 2, col = 1, lwd = 2)
#----- PL Confidence intervals
lines(AInt[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(AInt[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(AInt[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

#------------------------------- ---
# Figure 4c R*(beta) and R(beta):  ####
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
  main = "(4c)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 2
)
lines(c(1, 1), c(0, 1), lty = 2, col = 1, lwd = 2)
lines(betis, Rpbetas, type = "l", lty = 1, col = 1, lwd = 2)
# PL intervals :
lines(Betaint[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

#--------------------- --
# Figure 4d. Contours of R(a,beta) ####
Ccontours = clikelevels(n)$cbetasig
NC = 50 # number of points in grid
labetamat = matrix(rep(0, NC * NC), nrow = NC, ncol = NC)
# Values for axis limits for contour plot:
AA1 = Aendpt[1]
AA2 = Aendpt[2]
AAS = seq(AA1, AA2, length = NC)
BBet1 = Betaend[1]
BBet2 = Betaend[2]
BETAS = seq(BBet1, BBet2, length = NC)
# Log relative likelihood r(a,beta) is evaluated over grid:
for (i in 1:NC) {
  for (j in 1:NC) {
    bivec = c(AAS[i], BETAS[j])
    labetamat[i, j] = lvabetaWei(bivec, X, t, n)
  }
}

Rabetamat = exp(labetamat - LVmax) # R(a,beta) calculated over grid
Labelvec = c("90%", "95%", "99%") # contour labels

X11()
contour(AAS, BETAS, Rabetamat,
  levels = Ccontours, labels = Labelvec,
  labcex = 1, xlim = c(AA1, AA2), ylim = c(BBet1, BBet2),
  xlab = expression(a), ylab = expression(beta), cex.lab = 1, cex.axis = 1,
  main = "(4d)"
)
# The mles are marked with an asterisk over the contours:
points(amle, betamle, col = 1, pch = 8)

#--------------------- --
# Computing Weibull Median mles and intervals ####
# Median mles and profile likelihood intervals:
IntQW50 = IntervalQP(0.50, Y, n, amle, betamle)
#  Mles of quantile 0.50:
QG50mle = IntQW50$QGpmle
QW50mle = IntQW50$QWpmle
# Weibull profile likelihood intervals for the median:
IntsQW50 = IntQW50$WeibullQPIntervals


# __________________________ ####
# All Figures 4 a,b,c,d together: ####

X11()
split.screen(c(2, 2))
screen(1)
plot(Ebetas, uwei,
  pch = 19, cex = .5, ylab = "F(x) Transform",
  xlab = "Uniform (0,1) Quantiles", xlim = c(0, 1), ylim = c(0, 1),
  main = "(4a)", pty = "m", col = 1
)
lines(c(0, 1), c(0, 1), lty = 2, col = 1)
lines(Ebetas, uwei, lty = 2, col = 1)
lines(Ebetas, Ban1, lty = 2, col = 1)
lines(Ebetas, Ban2, lty = 2, col = 1)
#------------------- --
screen(2)
# R(a) in solid line:
plot(Ais, RpaExact,
  type = "l", lty = 1, ylim = c(0, 1), lwd = 2,
  xlim = c(Aendpt[1], Aendpt[2]),
  main = "(5b)", xlab = expression(a),
  ylab = expression(paste("R(", a, ")")),
  yaxs = "i"
)
# R*(a) in dashes:
lines(Ais, RpaSt, lty = 2, col = 1, lwd = 2)
#----- PL Confidence intervals
lines(AInt[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(AInt[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(AInt[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

#------------------- --
screen(3)
plot(betis, RpbetaSt,
  type = "l", ylim = c(0, 1), lty = 2,
  main = "(4c)", lwd = 2, ylab = expression(paste("R(", beta, ")")),
  xlab = expression(beta), yaxs = "i", col = 2
)
lines(c(1, 1), c(0, 1), lty = 2, col = 1, lwd = 2)
lines(betis, Rpbetas, type = "l", lty = 1, col = 1, lwd = 2)
# PL intervals :
lines(Betaint[1, ], rep(c90, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[2, ], rep(c95, 2), lty = 1, col = 1, lwd = 2)
lines(Betaint[3, ], rep(c99, 2), lty = 1, col = 1, lwd = 2)

#------------------- --
screen(4)
contour(AAS, BETAS, Rabetamat,
  levels = Ccontours, labels = Labelvec,
  labcex = 1, xlim = c(AA1, AA2), ylim = c(BBet1, BBet2),
  xlab = expression(a), ylab = expression(beta), cex.lab = 1, cex.axis = 1,
  main = "(4d)"
)
# The mles are marked with an asterisk over the contours:
points(amle, betamle, col = 1, pch = 8)


# ______________________________________ ####
#       Generating outpuf file with results:                 ####
sink(file = "Example4_2Output.txt")

cat("Example 4.2. Interrarival times between customers to a Museum")
cat("\n")
cat("and Figure 4 (a,b,c,d).")
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
cat("The Modified Anderson Darling Statistic for testing Weibull was:")
cat("\n")
print(ModAD)
cat("\n")
cat("to be compared with the following quantiles of its distribution")
cat("\n")
cat("when the Weibull is reasonable:")
cat("\n")
print(ModADQuantiles)
cat("\n")
cat("If ModAD is larger than some of these percentiles,")
cat("\n")
cat("then the Weibull/Gumbel are not reasonable models for the data/logdata")
cat("-------------------------------------")
cat("\n")
cat("\n")
cat(" PL intervals for beta:")
cat("\n")
print(Betaint)
cat("\n")
cat("PL intervals for sigma:")
cat("\n")
print(SigmaInt)
cat("\n")
cat("------------------------------------")
cat("\n")
cat(" Mles of Gumbel and Weibull Median are:")
cat("\n")
cat(paste(QG50mle, ",", QW50mle))
cat("\n")
cat("------------------------------------")
cat("\n")
cat("\n")
cat("confidence intervals for WEibull Median Q50:")
cat("\n")
print(IntsQW50)
cat("\n")
cat("Likelihood levels of intervals of a single parameter for n:")
cat("\n")
print(IntQW50$LikeLevels)
cat("\n")
cat("------------------------------------")
sink(file = NULL)
