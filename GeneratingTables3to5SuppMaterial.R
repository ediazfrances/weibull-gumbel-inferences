# Generating Tables 3 to 5 in the Supplementary Material OF TCH-23-068R1####
# Exactly those presented in the article are generated here:

#---------------------------------------- --
# Functions to be used: ####
#-------------------------------- --
# Agresti Coull interval for a Binomial proportion
AgrestiCoul = function(V, MM) {
  zp = qnorm(0.975, mean = 0, sd = 1)
  z2 = zp^2
  phat = (V + (z2 / 2)) / (MM + z2)
  SAC = zp * sqrt((phat * (1 - phat)) / (MM + z2))
  ACinterval = cbind(phat - SAC, phat + SAC)
  return(ACinterval)
}
# __________________________ ####
#--------------------------------- --
# Considered sample sizes:
Ns = seq(10, 100, by = 5)
#---------------------------- --
# Selected sample sizes ####
Indeces = c(1, 2, 3, 4, 6, 9, 14, 19)
NSelec = Ns[Indeces]
# number of total simulated samples
M = 10000

#-------------------------------- --
# Generating tables for Cov  Freqs of intervals for (a): ####
# Table 3: Coverage Freqa of PL(a)
# Table 4: Coverage Freqs of APL(a)
# Table 5: Coverage Freqs of WA(a)

#------------------------------ --
# Table 3 PL for (a): ####
#-------------------------------------- --
# Coverage probabilities of PL intervals for parameter (a) Gumbel,
# obtained in simulations:
PLa = c(
  9498, 9496, 9506, 9489, 9534, 9519, 9510, 9510, 9490, 9472,
  9469, 9491, 9488, 9488, 9537, 9514, 9489, 9480, 9495
)
# AgrestiCoull intervals for these binomial proportions, out of M:
PLaBin = AgrestiCoul(PLa, M)

# Selecting some sample sizes:
PLaSelec = PLa[Indeces]
PLaSelecBin = AgrestiCoul(PLaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderPLa = c(271, 291, 268, 268, 269, 269, 275, 270)
UnderPLaBin = AgrestiCoul(UnderPLa, M)
OverPLa = c(231, 213, 226, 243, 212, 241, 237, 235)
OverPLaBin = AgrestiCoul(OverPLa, M)


Table3PLa = cbind(NSelec, UnderPLaBin, PLaSelecBin, OverPLaBin)
Table3PLa = matrix(Table3PLa, nrow = 8, ncol = 7)
nomis = c(
  "N (a)", "UnderPL1", "UnderPL2", "CoverPL1", "CoverPL2",
  "OverPL1", "OverPL2"
)
colnames(Table3PLa) = nomis

#------------------------------ --
# Table 4 APL for (a): ####
#-------------------------------- --
# coverage probabilities for APL intervals for (a):
APLa = c(
  9474, 9484, 9493, 9478, 9529, 9511, 9485, 9508, 9484, 9469,
  9465, 9489, 9483, 9487, 9534, 9511, 9488, 9480, 9493
)
APLaBin = AgrestiCoul(APLa, M)
# Selecting some sample sizes:
APLaSelec = APLa[Indeces]
APLaSelecBin = AgrestiCoul(APLaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderAPLa = c(288, 298, 278, 275, 275, 272, 277, 271)
UnderAPLaBin = AgrestiCoul(UnderAPLa, M)
OverAPLa = c(238, 218, 229, 247, 214, 244, 236, 236)
OverAPLaBin = AgrestiCoul(OverAPLa, M)

Table4APLa = cbind(NSelec, UnderAPLaBin, APLaSelecBin, OverAPLaBin)
Table4APLa = matrix(Table4APLa, nrow = 8, ncol = 7)
nomis = c(
  "N (a)", "UnderAPL1", "UnderAPL2", "CoverAPL1", "CoverAPL2",
  "OverAPL1", "OverAPL2"
)
colnames(Table4APLa) = nomis

#------------------------------ --
# Table 5 for WA (a): ####
#-------------------------------- --
# coverage probabilities for WA intervals for (a):
WAa = c(
  8936, 9194, 9278, 9275, 9381, 9378, 9385, 9412, 9413, 9393, 9388, 9422,
  9414, 9422, 9479, 9468, 9435, 9433, 9455
)
WAaBin = AgrestiCoul(WAa, M)
# Selecting some sample sizes:
WAaSelec = WAa[Indeces]
WAaSelecBin = AgrestiCoul(WAaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderWAa = c(499, 409, 363, 349, 317, 285, 291, 274)
UnderWAaBin = AgrestiCoul(UnderWAa, M)
OverWAa = c(565, 397, 359, 376, 305, 302, 287, 271)
OverWAaBin = AgrestiCoul(OverWAa, M)

Table5WAa = cbind(NSelec, UnderWAaBin, WAaSelecBin, OverWAaBin)
Table5WAa = matrix(Table5WAa, nrow = 8, ncol = 7)
nomis = c(
  "N (a)", "UnderWA1", "UnderWA2", "CoverWA1", "CoverWA2",
  "OverWA1", "OverWA2"
)
colnames(Table5WAa) = nomis

# __________________________________________________ ####
#-------------------------------- --
# Generating tables for Cov Freqs of intervals for (Beta): ####
# Table 3b: Coverage Freqa of PL(Beta)
# Table 4b: Coverage Freqs of APL(Beta)
# Table 5b: Coverage Freqs of WA(Beta)

#------------------------------ --
# Table 3b: PL for Beta ####
#-------------------------------------- --
# Coverage probabilities of PL intervals for Beta,
# obtained in simulations:
PLBeta = c(
  9507, 9516, 9526, 9502, 9482, 9522, 9527, 9498, 9475, 9510,
  9504, 9502, 9498, 9513, 9514, 9476, 9481, 9481, 9493
)
# Agresti-Coull's intervals for the corresponding binomial proportion:
# These intervals estimate the confidence level of these intervals.
PLBetaBin = AgrestiCoul(PLBeta, M)

# Selecting some sample sizes:
PLBetaSelec = PLBeta[Indeces]
PLBetaSelecBin = AgrestiCoul(PLBetaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderPLBeta = c(78, 97, 105, 136, 131, 167, 177, 175)
UnderPLBetaBin = AgrestiCoul(UnderPLBeta, M)
OverPLBeta = c(415, 387, 369, 362, 347, 358, 310, 332)
OverPLBetaBin = AgrestiCoul(OverPLBeta, M)

Table3bPLBeta = cbind(NSelec, UnderPLBetaBin, PLBetaSelecBin, OverPLBetaBin)
Table3bPLBeta = matrix(Table3bPLBeta, nrow = 8, ncol = 7)
nomis = c(
  "N (Beta)", "UnderPL1", "UnderPL2", "CoverPL1", "CoverPL2",
  "OverPL1", "OverPL2"
)
colnames(Table3bPLBeta) = nomis

#------------------------------ --
# Table 4b: APL for Beta ####
#-------------------------------------- --
# coverage frequencies of APL intervals for beta:
APLBeta = c(
  9501, 9509, 9520, 9497, 9481, 9523, 9524, 9498, 9472, 9508,
  9503, 9501, 9498, 9515, 9512, 9472, 9482, 9478, 9493
)
APLBetaBin = AgrestiCoul(APLBeta, M)

# Selecting some sample sizes:
APLBetaSelec = APLBeta[Indeces]
APLBetaSelecBin = AgrestiCoul(APLBetaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderAPLBeta = c(73, 94, 103, 133, 129, 167, 176, 175)
UnderAPLBetaBin = AgrestiCoul(UnderAPLBeta, M)
OverAPLBeta = c(426, 397, 377, 370, 348, 361, 309, 332)
OverAPLBetaBin = AgrestiCoul(OverAPLBeta, M)

Table4bAPLBeta = cbind(NSelec, UnderAPLBetaBin, APLBetaSelecBin, OverAPLBetaBin)
Table4bAPLBeta = matrix(Table4bAPLBeta, nrow = 8, ncol = 7)
nomis = c(
  "N (Beta)", "UnderAPL1", "UnderAPL2", "CoverAPL1", "CoverAPL2",
  "OverAPL1", "OverAPL2"
)
colnames(Table4bAPLBeta) = nomis

#------------------------------ --
# Table 5b: WA for Beta ####
#-------------------------------------- --
# coverage probabilities for WALD For a, denoted as WA:
WABeta = c(
  9531, 9526, 9534, 9493, 9482, 9525, 9543, 9506, 9483, 9511, 9512,
  9510, 9496, 9517, 9509, 9488, 9483, 9477, 9490
)
WABetaBin = AgrestiCoul(WABeta, M)
# Selecting some sample sizes:
WABetaSelec = WABeta[Indeces]
WABetaSelecBin = AgrestiCoul(WABetaSelec, M)
# Number of intervals  that under and over estimated A0:
# for selected sample sizes: N= 10,15,20,25,30,35,50,75,100
UnderWABeta = c(187, 188, 185, 220, 209, 220, 220, 215)
UnderWABetaBin = AgrestiCoul(UnderWABeta, M)
OverWABeta = c(282, 286, 281, 287, 266, 297, 263, 295)
OverWABetaBin = AgrestiCoul(OverWABeta, M)

Table5bWABeta = cbind(NSelec, UnderWABetaBin, WABetaSelecBin, OverWABetaBin)
Table5bWABeta = matrix(Table5bWABeta, nrow = 8, ncol = 7)
nomis = c(
  "N (Beta)", "UnderWA1", "UnderWA2", "CoverWA1", "CoverWA2",
  "OverWA1", "OverWA2"
)
colnames(Table5bWABeta) = nomis


# ______________________________________ ####
# Generating output Tables 3 to 5 in Tables3to5.txt ####
sink(file = "Tables3to5.txt")
cat("Tables 3 to 5. Agresti-Coull's confidence intervals for binomial proportions")
cat("\n")
cat("of PL, APL, and WA intervals for (a) and (beta) for selected sample sizes:")
cat("\n")
cat("\n")
cat("\n")
cat("------------------------------------")
cat("\n")
cat("Tables 3a. AC intervals for PL intervals for a:")
cat("\n")
print(Table3PLa)
cat("\n")
cat("Table 3b. AC intervals for PL intervals for beta:")
cat("\n")
print(Table3bPLBeta)
cat("\n")
cat("Table 4a. AC intervals for APL intervals for a:")
cat("\n")
print(Table4APLa)
cat("\n")
cat("Table 4b. AC intervals for APL intervals for Beta:")
cat("\n")
print(Table4bAPLBeta)
cat("\n")
cat("Table 5a. AC intervals for WA intervals for a:")
cat("\n")
print(Table5WAa)
cat("\n")
cat("Table 5b. AC intervals for WA intervals for Beta:")
cat("\n")
print(Table5bWABeta)
cat("\n")
cat("___________________________________________________________________")
sink(file = NULL)
