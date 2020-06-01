`few.obs` <- 
function(survey, by.char = NULL)
######################################################
# Controlla il numero di osservazioni contenute nel  #
# dataframe dei dati campionari 'survey' e nei suoi  #
# subset determinati dalle variabili 'by.char':      #
# se tale numero e' MINORE DI 10 genera un warning.  #
######################################################
{
if (identical(by.char, NULL))
    {
     if (nrow(survey) < 10) 
         warning("Less than 10 observations in survey data frame")
     return(invisible(NULL))
    }
survey.by <- split(survey, survey[, by.char], drop = TRUE)	
sapply(seq_along(along.with = survey.by),
       function(i) if (nrow(survey.by[[i]]) < 10) 
                       warning("Less than 10 observations from subpopulation: ",
                                names(survey.by)[i], call. = FALSE))
return(invisible(NULL))
}