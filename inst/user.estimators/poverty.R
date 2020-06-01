`poverty` <-
function (d, w, y, threshold)
######################################################################################
#  Population percentage with income below the poverty threshold. Suppose poverty    #
#  threshold ('threshold' in the function body) is defined as 0.6 times the average  #
#  income for the whole population.                                                  #
#  NOTE: since the poverty status of each final unit depends on a global value       #
#        (that is, the average income for the whole population) 'global' function    #
#        is used to prevent, whenever a sub-population poverty estimate is needed,   #
#        this global value being calculated locally i.e. within the sub-population   #
#        itself.                                                                     #
######################################################################################
{
    if (missing(threshold)) {
    # if I do want to take into account the variance of the poverty
    # threshold, I let it be re-calculated replicate by replicate.
        d.global = global(d)
        th.value = 0.6 * sum(d.global[, w] * d.global[, y])/sum(d.global[, w])
    }
    else {
    # if I do not want to take into account the variance of the poverty
    # threshold, I will supply its point estimate to the 'threshold' argument.
        th.value = threshold 
    }
    est = 100 * sum(d[d[, y] < th.value, w])/sum(d[, w])
    est
}