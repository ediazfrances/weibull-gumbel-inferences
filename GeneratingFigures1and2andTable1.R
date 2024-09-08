# Generating Figures 1 and 2 as well as Table 1.  ####
# After running this scirpt file check text file Table1.txt

# _______________________________________________ ####
# Cleaning Workspace (recommended):
# rm(list=ls(all=TRUE))
#---------------------------------------- --

# Functions to be used: ####
#-------------------------------- --
# Agresti-Coul interval for binomial proportion ####

AgrestiCoul = function(V, MM) {
  zp = qnorm(0.975, mean = 0, sd = 1)
  z2 = zp^2
  phat = (V + (z2 / 2)) / (MM + z2)
  SAC = zp * sqrt((phat * (1 - phat)) / (MM + z2))
  ACinterval = cbind(phat - SAC, phat + SAC)
  return(ACinterval)
}

# __________________________________________________ ####

#--------------------------------------- --

#--- FIgs 1  Plot: ####
# Coverage frequencies of PL, APL, and WA intervals for a
# in M simulated samples.
# Considered sample sizes:
enitas = seq(10, 100, by = 10)
# number of total simulated samples
M = 10000
Ns = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
# n=Ns[8]
# Coverage probs of proposed intervals PL for a Gumbel:
PLa = c(
  9498, 9496, 9506, 9489, 9534, 9519, 9510, 9510, 9490, 9472,
  9469, 9491, 9488, 9488, 9537, 9514, 9489, 9480, 9495
)
# van intervalos AGrestiCoul para proporción binomial:
PLaBin = AgrestiCoul(PLa, M)
#---------------------- --
# coverage probabilities for APL asymptotic prof like intervals for a:
APLa = c(
  9474, 9484, 9493, 9478, 9529, 9511, 9485, 9508, 9484, 9469,
  9465, 9489, 9483, 9487, 9534, 9511, 9488, 9480, 9493
)
APLaBin = AgrestiCoul(APLa, M)
#---------------------- --
# coverage probabilities for WALD For a, denoted as WA:
WAa = c(
  8936, 9194, 9278, 9275, 9381, 9378, 9385, 9412, 9413, 9393, 9388, 9422,
  9414, 9422, 9479, 9468, 9435, 9433, 9455
)
WAaBin = AgrestiCoul(WAa, M)

# Figure 1 Cov probs of proposed intervals for a vs sample size
X11()
plot(c(5, 105), c(0.95, 0.95),
  type = "l", lty = 1,
  ylab = "Coverage Prob. (a)",
  xlab = "Sample Size n", xlim = c(5, 105), ylim = c(0.885, 0.962), yaxs = "i",
  xaxs = "i"
)
axis(side = 1, at = enitas, tick = TRUE)
for (i in (1:19)) {
  lines(c(Ns[i], Ns[i]), PLaBin[i, ], lty = 1, col = 1)
  points(c(Ns[i], Ns[i]), APLaBin[i, ], pch = 8, col = 2)
  points(c(Ns[i], Ns[i]), WAaBin[i, ], pch = 4, col = 4)
  lines(c(Ns[i], Ns[i]), WAaBin[i, ], lty = 2, col = 2)
}
lines(c(5, 105), c(0.96, 0.96), lty = 2, col = 1)
lines(c(5, 105), c(0.94, 0.94), lty = 2, col = 1)

#------------------------------------ --
#--- FIgs 2  Plot: ####
# Beta: Coverage frequencies of PL, APL, and WA intervals for beta
# in M simulated samples with different sample sizes n.
#------------------------------------ --
# considered sample sizes n
enitas = seq(10, 100, by = 10)
# number of total simulated Gumbel samples per n.
M = 10000
# Coverage frequencies of PL intervals for beta:
# number of intervals that included true value of beta:
PLBeta = c(
  9507, 9516, 9526, 9502, 9482, 9522, 9527, 9498, 9475, 9510,
  9504, 9502, 9498, 9513, 9514, 9476, 9481, 9481, 9493
)
# Agresti-Coull's intervals for the corresponding binomial proportion:
# These intervals estimate the confidence level of these intervals.
PLBetaBin = AgrestiCoul(PLBeta, M)
#---------------------- --
# coverage frequencies of APL intervals for beta:
APLBeta = c(
  9501, 9509, 9520, 9497, 9481, 9523, 9524, 9498, 9472, 9508,
  9503, 9501, 9498, 9515, 9512, 9472, 9482, 9478, 9493
)
APLBetaBin = AgrestiCoul(APLBeta, M)
#---------------------- --
# coverage probabilities for WALD For a, denoted as WA:
WABeta = c(
  9531, 9526, 9534, 9493, 9482, 9525, 9543, 9506, 9483, 9511, 9512,
  9510, 9496, 9517, 9509, 9488, 9483, 9477, 9490
)
WABetaBin = AgrestiCoul(WABeta, M)
#----------------- --
# Figure 2 plot. Coverage probabilities of PL, APL, and WA intervals for beta:
X11()
plot(c(5, 105), c(0.95, 0.95),
  type = "l", lty = 1,
  ylab = expression(paste("Coverage Prob. (", beta, ")")),
  xlab = "Sample Size n", xlim = c(5, 105), ylim = c(0.938, 0.962), yaxs = "i",
  xaxs = "i"
)
axis(side = 1, at = enitas, tick = TRUE)
for (i in (1:19)) {
  lines(c(Ns[i], Ns[i]), PLBetaBin[i, ], lty = 1, col = 1)
  points(c(Ns[i], Ns[i]), APLBetaBin[i, ], pch = 8, col = 2)
  points(c(Ns[i], Ns[i]), WABetaBin[i, ], pch = 4, col = 4)
  # lines(c(Ns[i],Ns[i]),BetaWBin[i,],lty=2,col=2)
}
lines(c(5, 105), c(0.96, 0.96), lty = 2, col = 1)
lines(c(5, 105), c(0.94, 0.94), lty = 2, col = 1)


# ______________________________ ####

# Table 1 is calculated here:  ####
# Sample sizes considered for Table 1:
NS = c(10, 15, 20, 25, 35, 50, 75, 100)
# Coverage frequencies of proposed likelihood region LR of (a,beta):
regprop = c(9484, 9489, 9499, 9474, 9496, 9492, 9482, 9478)
# Agresti-Coull's intervals for estimating their coverage probabilities:
RegpropBin = AgrestiCoul(regprop, M)

# Coverage frequencies of Asymptotic likelihood regions (AR) of (a,beta):
regasint = c(9328, 9384, 9436, 9427, 9457, 9464, 9467, 9468)
# Agresti-Coull's intervals for estimating their coverage probabilities:
RegasintBin = AgrestiCoul(regasint, M)

#------------------------------------ --
# Coverage frequencies of  profile likelihood intervals of
# Gumbel quantile QG10:
IntQpPL = c(9525, 9521, 9548, 9524, 9525, 9469, 9498, 9496)
# Agresti-Coull's intervals for estimating their coverage probabilities:
IntQpPLBin = AgrestiCoul(IntQpPL, M)

# Coverage frequencies of asymptotic profile likelihood (APL) intervals of
# Gumbel quantile QG10
IntQpAPL = c(9305, 9390, 9446, 9447, 9468, 9433, 9472, 9476)
# Agresti-Coull's intervals for estimating their coverage probabilities:
IntQpAPLBin = AgrestiCoul(IntQpAPL, M)

# Summary of coverage probabilities for Gumbel quantile QG10 PL and APL:
CovPRegionsQG10 = cbind(NS, RegpropBin, RegasintBin, IntQpPLBin, IntQpAPLBin)

# ______________________________________ ####
#       Generating outpuf file with results:                 ####
sink(file = "Table1.txt")
cat("Table 1. Agresti-Coull's confidence intervals for binomial proportions")
cat("\n")
cat("of proposed and asymptotic Likelihood Regions, LR and AR, and")
cat("\n")
cat("of PL intervals for Gumbel quantile QG10 and APL intervals.")
cat("\n")
cat("Col 1= N, Cols 2-3= LR, Cols 4-5= AR")
cat("\n")
cat(" Col 5-6=QG10 PL, Cols 7-8= QG10 APL.")
cat("\n")
print(CovPRegionsQG10)
cat("\n")
cat("____________________________________________")
cat("------------------------------------")
sink(file = NULL)
