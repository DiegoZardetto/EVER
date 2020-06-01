# EVER <img src="EVER_LOGO_small.png" align="right" alt="" />

**EVER** (Estimation of Variance by Efficient Replication) is an R package for calibration, estimation and sampling error assessment in complex sample surveys. 

EVER's sampling variance estimation is based on the extended **DAGJK (Delete-A-group Jackknife)** technique proposed by P. S. Kott.


## Installation
You can install the **development version** of EVER from [GitHub](https://github.com/DiegoZardetto/EVER) as follows:

```r
install.packages("devtools")
devtools::install_github("DiegoZardetto/EVER")
```

The **last released version** of EVER can be downloaded from [Istat website](https://www.istat.it/en/methods-and-tools/methods-and-it-tools/process/processing-tools/ever).


## Main Statistical Functions
  - #### Delete-A-Group Jackknife replication
  - #### Calibration of replicate weights
  - #### Estimates and Sampling Errors (standard error, variance, coefficient of variation, confidence interval, design effect) for:
    - Totals
    - Means
    - Absolute and relative frequency distributions
    - Ratios between totals
    - Multiple regression coefficients
    - Quantiles
  - #### Estimates and Sampling Errors for user-defined Complex Estimators (even non-analytic)
  - #### Estimates and Sampling Errors for Subpopulations (Domains)
    - All the analyses above can be carried out for arbitrary domains


## Sampling Variance Estimation Methodology
The advantage of the DAGJK method over the traditional jackknife is that, unlike the latter, it remains computationally manageable even when dealing with "complex and big" surveys (tens of thousands of PSUs arranged in a large number of strata with widely varying sizes). In fact, the DAGJK method is known to provide, for a broad range of sampling designs and estimators, (near) unbiased standard error estimates even with a "small" number (e.g. a few tens) of replicate weights.

Besides its peculiar computational efficiency, the DAGJK method takes advantage of the strong points it shares with the most common replication methods. As a remarkable example, EVER is designed to fully exploit DAGJK's versatility: the package provides the user with a user-friendly tool for calculating estimates, standard errors and confidence intervals for estimators defined by the user themselves (even non-analytic). This functionality makes EVER especially appealing whenever variance estimation by Taylor linearisation can be applied only at the price of crude approximations (e.g. poverty estimates).


## References
Kott, Phillip S. (1999) "The Extended Delete-A-Group Jackknife". Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) "The Delete-A-Group Jackknife". Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
