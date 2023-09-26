#' twofish - Fisher's Exact Test on 2X2 Contingency Tables
#'
#' @param table a 2X2 contingency table
#' @param alpha a confidence coefficient (default = .05)
#' @param alt the type of alternative hypothesis {'less', 'two.sided', 'greater'} ... default = two.sided
#' @param verbose A Boolean for whether or not to print extra information (default = TRUE)
#' @param visual A Boolean for whether or not to construct plots (default = TRUE)
#'
#' @return an invisible list called 'tf' that contains various tables,
#' the calculated p-value,and the sample estimate of the odds ration
#' @export
#' @importFrom stats chisq.test dhyper
#'
#' @examples
#' twofish(table = table(c(1,2),c(1,2)))
twofish <- function(table,
                    alpha = .05,
                    alt = 'greater',
                    verbose = TRUE,
                    visual = TRUE){

  # Constructing relevant tables

  mtab <- addmargins(table)
  ptab <- addmargins(prop.table(table))
  extab <- suppressWarnings(chisq.test(table)$expected)

  # Storing expected values in a vector

  ex <- c(extab[1,1],extab[2,1],extab[1,2],extab[2,2])

  # Printing additional text based on user's choice

  if(verbose){

    print("---------------------------------------------")
    print("    Observed Counts for Contingency Table    ")
    print("---------------------------------------------")
    print(mtab)

    cat("\n")

    print("---------------------------------------------")
    print("     Probabilities for Contingency Table     ")
    print("---------------------------------------------")
    print(round(ptab,4))

    cat("\n")

    print("---------------------------------------------")
    print("      Contingency Table W/ Expectations      ")
    print("---------------------------------------------")

    cat("\n")

    print(round(extab,4))

    cat("\n")

    # Testing primary assumption (Expectation > 5?)

    STOP = FALSE

    for (i in 1:length(ex)){
      if (ex[i] < 5){
        print("The data is not suitable for an asymptotic X-squared test.")
        print("The p_value calculated below can be used to check against the Null.")
        cat("\n")
        STOP = TRUE
        break
      }
      if (STOP){break}
    }
  }

  # Displaying plots and visualizations based on user's selection

  if (visual){

    y_upper = max(ex)*1.25
    y_label = max(ex)*1.1

    ex_plot <- barplot(ex,
                       col = c("cyan", "violet", "cyan", "violet"),
                       xlab = "Cells",
                       ylab = "Counts",
                       main = "Expected Cell Counts",
                       ylim = c(0,y_upper))

    abline(h=5, lwd = 2.75)

    text(ex_plot,
         y_label,
         round(ex,2),
         col = 'black',
         cex=2)

    mosaicplot(table,
               main = "Mosaic plot",
               color = c("violet","cyan"))
  }

  # Storing table values in easy to reference objects

  a <- as.numeric(mtab[1,1])
  b <- as.numeric(mtab[1,2])
  c <- as.numeric(mtab[1,3])
  d <- as.numeric(mtab[2,1])
  e <- as.numeric(mtab[2,2])
  f <- as.numeric(mtab[2,3])
  g <- as.numeric(mtab[3,1])
  h <- as.numeric(mtab[3,2])
  i <- as.numeric(mtab[3,3])

  # Printing Summary

  print("---------------------------------------------")
  print("                   Summary                   ")
  print("---------------------------------------------")

  cat("\n")

  print("Null Hypothesis:")
  print("The two directions of data classification are dependent.")

  cat("\n")

  print("Alternative Hypothesis:")
  print("The two directions of data classification are independent.")

  cat("\n")

  # Calculating Sample Estimate of Odds Ratio

  OR <- (a*e)/(b*d)

  # Calculating p-value base on user selection of alternative hypothesis

  p_value <- 0

  if (alt == 'less'){
    print("Type of Alternate Hypothesis: Less")
    p_value <- sum(dhyper(0:a,g,h,c))
  }

  if (alt == 'greater'){
    print("Type of Alternate Hypothesis: Greater")
    p_value <- sum(dhyper(3:g,g,h,c))
  }

  if (alt == 'two.sided'){
    probs <- dhyper(0:g,g,h,c)
    print("Type of Alternate Hypothesis: Two-Sided")
    p_value <- sum(probs[probs <= dhyper(a,g,h,c)])
  }

  cat("\n")

  print(paste("Sample Estimate of Odds Ratio: ", round(OR,4), sep =""))

  cat("\n")

  print(paste("Fisher's Exact P-Value: ", round(p_value,4), sep=""))
  print("(Probability of being at least as extreme as observed table in direction of alternative)")

  cat("\n")

  if (p_value < alpha){
    print(paste("Given the confidence interval provided(", alpha, "), we cannot accept the null hypothesis", sep=""))
    print("It is plausible that the two classifications are dependent.")
  }

  if (p_value > alpha){
    print(paste("Given the confidence interval provided(", alpha, "), we cannot reject the null hypothesis", sep=""))
    print("It is plausible that the two classifications are independent.")
  }

  tf <- list("Data.Table" = table,
             "Contingency.Table" = mtab,
             "Probability.Table" = ptab,
             "Expecation.Table" = extab,
             "p.value" = p_value,
             "Sample.OR" = OR)

  invisible(tf)

}
