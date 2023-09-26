#' Helper function for computing the difference in cell proportions within my1way()
#'
#' @param table a one way table
#' @param confidence the desired level of confidence
#' @param adjusted whether or not to do Bonferoni adjustment
#'
#' @return nothing - only prints to console
#' @export
#'
#' @examples
#' cell_diff_proportion_CIs(table = table(c(1,2),c(1,2)))
cell_diff_proportion_CIs <- function(table, confidence = .95, adjusted = FALSE){

  ######## Modifying Table ########

  # Constructing a copy with margins
  tab_margins <- addmargins(table)

  # Constructing a copy with proportions
  tab_prop <- prop.table(table)

  ######## Declaring Variables ########

  # Table Number Columns
  numcols <- length(table)

  # Table Row Sum
  sum <- tab_margins[numcols + 1]

  # Alpha Value
  alpha <- 1 - confidence

  # Z Score
  z <- qnorm(1-alpha/2,0,1)

  # Bonferroni Corrected Z Score
  m <- numcols* (numcols-1) / 2
  Bz <- qnorm(1-alpha/(2*m),0,1)

  ######## Computing Intervals (Adjusted == FALSE) ########

  if (adjusted == FALSE){

    for (i in 1:(numcols-1)){

      for (j in (i+1):(numcols)){

        if (i == 1 & j == 2){
          print('----------------------------------------------------------')
          print(paste(confidence, '% Confidence Intervals For Cell Proportion Differences', sep = ""))
          print('----------------------------------------------------------')
        }
        p1 <- tab_prop[i]
        p2 <- tab_prop[j]
        diff <- tab_prop[i] - tab_prop[j]
        numerator <- (p1 * (1 - p1)) + (p2 * (1 - p2)) + (2 * p1 * p2)
        denominator <- sum
        lower = diff - (z * sqrt(numerator / denominator))
        upper = diff + (z * sqrt(numerator / denominator))
        print(paste(rownames(table)[i],
                    rownames(table)[j],
                    ' (', round(lower, 5),
                    ', ',
                    round(upper, 5),
                    ')',
                    sep =""))
      }
    }
  }
  ######## Computing Intervals (Adjusted == TRUE) ########

  if (adjusted == TRUE) {

    for (i in 1:(numcols-1)){

      for (j in (i+1):(numcols)){

        if (i == 1 & j == 2){
          print('----------------------------------------------------------')
          print(paste(confidence, '% Confidence Intervals For Cell Proportion Differences', sep = ""))
          print(paste('         (Bonferonni Adjusted For ', m, ' Hypothesis)          ', sep =""))
          print('----------------------------------------------------------')
        }
        p1 <- tab_prop[i]
        p2 <- tab_prop[j]
        diff <- tab_prop[i] - tab_prop[j]
        numerator <- (p1 * (1 - p1)) + (p2 * (1 - p2)) + (2 * p1 * p2)
        denominator <- sum
        lower = diff - (Bz * sqrt(numerator / denominator))
        upper = diff + (Bz * sqrt(numerator / denominator))
        print(paste(rownames(table)[i],
                    rownames(table)[j],
                    ' (', round(lower, 5),
                    ', ',
                    round(upper, 5),
                    ')',
                    sep =""))
      }
    }
  }
}
