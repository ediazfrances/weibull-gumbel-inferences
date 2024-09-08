Notes for general users of the proposed Weibull and Gumbel inferences, described in "Practical and Optimal Likelihood Intervals and Regions for Weibull and Gumbel Distributions" by E. Díaz-Francés (diazfran@cimat.mx) 

[TCH23_068R2] Updated by E. Díaz-Francés on September 7, 2024.

1) IMPORTANT: The user must source the script file WeibullFunctions.R first since it contains all functions that will be used. Each function is described in comments within this script file.

2) Fill an R numeric vector (called Data here) with the observed data. The entries in the data vector must all consist of strictly positive numbers. The sample size N must be greater or equal to ten. In case the user has a data vector X assumed to follow a Gumbel of minima distribution, then the corresponding exponentiated vector Y = exp(X) that follows a Weibull distribution must be fed as Data = exp(X) to the following functions detailed in points 3 and 4.

3) First, one must check if the Weibull distribution is a reasonable model for the considered data vector by running the function:  Weibullvalidation(Data)
This function produces the mles of the Weibull parameters, the Modified Anderson Darling (ModAD) statistic of Stephens (1977), and the proposed Weibull probability plot for the observed data. Overall, if the ModAD statistic is larger than 0.637 the Weibull distribution might not be reasonable for the data, but if it is larger than 1, the model is highly questionable. In both cases, the probability plot is also informative about any lack of fit. 

4) Provided the Weibull model is a reasonable fit for the data, then one can proceed to obtain the proposed inferences by running the function: 
Weibullinferences(Data, QpProb)
The first argument is the data vector. The second argument is the probability of the quantile of interest of the Weibull distribution to be estimated and recommended to be a value within the interval [1/(N+1), N/(N+1)], where N is the sample size. By default this function will yield inferences for the quantile of probability 0.50, the median, as well as all proposed inferences for Weibull parameters, PL intervals and likelihood regions, and plots of relative likelihoods —global and profile— , as described in the article. 
Two folders will be created in the current working directory named WeibullResults and WeibullPlots. All proposed inferences will be exported to the first folder within the output text file named: Weibull_Inferences_Output.txt. The proposed plots will be exported as png files to the second folder.

5) The coverage frequencies of the proposed PL intervals as well as those of APL and WA intervals for a desired single set of values of true parameters (A,BETA) and a quantile of interest of probability P with sample size N can be checked by running the function Coverprobintervals(A, BETA, N, P). The coverage frequencies of the proposed likelihood regions are computed too. The recommended ranges for the true parameters that describe many and most of real world situations, are: A within [-3,10], BETA within [0.3,30], the probability P of the Gumbel quantile of interest QGp (or Weibull quantile QWp) must lie within [1/(N+1),N/(N+1)], and the sample size must be larger or equal to 10. The default values for this function are A=2, BETA=3, N=20, P=0.25.

6) For a new given setting of the true parameters of interest (A,BETA,P) the coverage frequencies of the previous point can also be calculated for 19 sample sizes as those considered in the article by running the function WeibullSimulations(a0,beta0,p). The default values are: (a0 = 2, beta0 = 3, p = 0.25). Equivalent figures as Figures 1 and 2 are produced as well as Tables 1, 3,4,5 for the selected new setting.

7) By running the following script files the examples and figures presented in the article are reproduced. Inferences for each example are exported to files Example4_1Output.txt, Example4_2Output.txt, and Example4_3Output.txt. 
a) GeneratingFigures1and2andTable1.R
b) GeneratingFigure3_Example4_1.R
c) GeneratingFigure4Example4_2.R
d) GeneratingFigure5Example4_3.R
e) GeneratingTables3to5SuppMaterial.R
_________________________________

EXAMPLE:

1. The file WeibullFunctions.R must be sourced first.
       The following data vector Y of size n = 14 was used in Example 4.1: 
              Y = c (55, 320, 56, 104, 220, 239, 47, 246, 176, 182, 33, 15, 104, 35) 

Assuming that the quantile of probability 0.25 is of interest, in order to obtain the proposed Weibull inferences for this data vector and the plots shown here at the end, the user can run the functions described in the following points with Y.

2. First, the Weibull model is validated for the data set by running the following function: 

Weibullvalidation (Data = Y ).

This function yields the mles 4.96 for a, 1.36 for beta, and ModAD = 0.447 for the data Y. Therefore, the Weibull distribution appears to be reasonable for this data set. Only when ModAD is smaller than one should one proceed to the next step; otherwise the Weibull model is clearly not reasonable for the data. The mles and probability plot that were produced by this function will be presented again as part of the results of the function described next.

3. Weibull inferences results for the provided data vector Y are obtained by running: Weibullinferences ( Data = Y, QpProb = 0.25 ). Results can be checked as the contents of two folders that are created in the current working directory named WeibullResults and WeibullPlots. Inferences for the Weibull vector Y are equivalent to Gumbel inferences for the log data X=lnY.  

