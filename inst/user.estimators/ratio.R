`ratio` <-
function (d, w, num, den)
###########################################
#  Ratio estimator for totals (or means)  #
#  of quantitative variables.             #
###########################################
{
    sum(d[, w] * d[, num])/sum(d[, w] * d[, den])
}