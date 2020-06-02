# EVER <img src="EVER_LOGO_small.png" align="right" alt="" />

**EVER (Estimation of Variance by Efficient Replication)** is an R package for calibration, estimation and sampling error assessment in complex sample surveys. 

EVER's sampling variance estimation is based on the extended **DAGJK (Delete-A-group Jackknife)** technique proposed by P. S. Kott (see References).


## Installation
You can install the **development version** of EVER from [GitHub](https://github.com/DiegoZardetto/EVER) as follows:

```r
install.packages("devtools")
devtools::install_github("DiegoZardetto/EVER")
```

The **last released version** of EVER can be downloaded from [Istat website](https://www.istat.it/en/methods-and-tools/methods-and-it-tools/process/processing-tools/ever).


## Website
A [pkgdown website](https://DiegoZardetto.github.io/EVER) for the EVER package is available here:
- <https://DiegoZardetto.github.io/EVER>


## References
Kott, Phillip S. (1999) "The Extended Delete-A-Group Jackknife". Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) "The Delete-A-Group Jackknife". Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
