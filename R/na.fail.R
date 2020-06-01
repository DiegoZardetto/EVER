`na.fail` <-
function (data, varnames)
###########################################################
# Controlla che le variabili i cui nomi sono specificati  #
# da 'varnames', SUPPOSTE PRESENTI nel dataframe 'data'   #
# siano prive di valori mancanti (NA).                    #
###########################################################
{
    has.na <- sapply(varnames, function(var) any(is.na(data[, 
        var])))
    if (any(has.na)) 
        stop("Missing values in: ", paste(varnames[has.na], 
              collapse = ", "), "\n")
}