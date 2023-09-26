#' mycontingency - Chi-squared test for contingency tables
#'
#' @param table an rXc contingency table
#' @param alpha A confidence coefficient (default = .05)
#' @param yates A Boolean for whether or not to use Yates Correction (default = TRUE)
#' @param verbose A Boolean for whether or not to print extra information (default = TRUE)
#' @param visual A Boolean for whether or not to construct plots (default = TRUE)
#'
#' @return an invisible list called 'mc' that contains various tables,
#' an exhaustive list of all the values related to the chi-squared test
#' @export
#'
#' @examples
#' mycontingency(table = table(c(1,2),c(1,2)))
mycontingency <- function(table,
                          alpha = .05,
                          yates = TRUE,
                          verbose = TRUE,
                          visual = TRUE){

  cat("\n")

  if (verbose){
    print("---------------------------------------------")
    cat("\n")
    print("                DATA SUMMARY                 ")
    cat("\n")
    print("---------------------------------------------")
    cat("\n")
  }

  tab <- table

  if (verbose){
    print("XTABS TABLE FROM GIVEN FORMULA")
    print(tab)
    cat("\n")
  }

  if(visual){
    mosaicplot(tab,
               main = "Mosaic plot",
               color = cm.colors(length(colnames(table))))
  }

  mtab <- addmargins(tab)
  if (verbose){
    print("OBSERVED COUNTS FOR CONTINGENCY TABLE")
    print(mtab)
    cat("\n")
  }

  ptab <- addmargins(prop.table(tab))
  if (verbose){
    print("PROBABILITIES FOR CONTINGENCY TABLE")
    print(round(ptab,4))
    cat("\n")
  }

  numcols <- length(colnames(tab))
  numrows <- length(rownames(tab))

  total <- as.numeric(mtab[numrows+1,numcols+1])

  # Creating a vector of row sums
  r_sums <- c()
  for (i in 1:numrows)
    r_sums <- c(r_sums, mtab[i,numcols+1])

  # Creating a vector of column sums
  c_sums <- c()
  for (i in 1:numcols)
    c_sums <- c(c_sums, mtab[numrows+1,i])

  # Computing Chi_calc
  chi_calc <- 0
  expected_val <- c()
  assumption = TRUE
  for (i in 1:numrows){
    for (j in 1:numcols){
      val <- tab[i,j]
      expectation <- (r_sums[i] * c_sums[j]) / total
      expected_val <- c(expected_val, expectation)
      if (expectation < 5)
        assumption = FALSE
      new <- (val - expectation)^2 / expectation
      chi_calc <- chi_calc + new
    }
  }

  # Computing Yates Corrected Chi-sqaure
  chi.yates <- 0
  expected_val <- c()
  for (i in 1:numrows){
    for (j in 1:numcols){
      val <- tab[i,j]
      expectation <- (r_sums[i] * c_sums[j]) / total
      expected_val <- c(expected_val, expectation)
      new <- (abs(val - expectation)-.5)^2 / expectation
      chi.yates <- chi.yates + new
    }
  }


  # Computing Chi_crit

  df <- (numrows - 1) * (numcols - 1)
  chi_test <- qchisq(p = 1-alpha, df = df)

  if(yates){
    p_value <- pchisq(q = chi.yates, df = df,lower.tail = FALSE)
  } else{
    p_value <- pchisq(q = chi_calc, df = df,lower.tail = FALSE)
  }

  ######################## PRINTING REPORT ###########################

  cat("\n")
  print("---------------------------------------------")
  cat("\n")
  print("                 TEST SUMMARY                ")
  cat("\n")
  print("---------------------------------------------")
  cat("\n")
  if(verbose){
    # Checking primary assumption
    if (assumption == TRUE){
      print("Remarks:")
      print('The estimated expected value of each cell is at least five')
      print('The critical assumption for this test holds.')
      print("An asymptotic X-squared test is valid for this data.")
    }
    if (assumption == FALSE){
      print("Remarks:")
      print('One or more of the estimated expected values of a cell is less than five')
      print('The critical assumption for this test does not hold.')
      print("An asymptotic X-squared test is not suggested for this data.")
      print("Consider using Fisher's exact test instead.")
    }

    cat("\n")
  }

  print("Test Values:")

  if(!yates)
    print(paste('Calculated X-squared Value: ', sprintf("%.5f", round(chi_calc, 5)), sep = ""))

  if (yates)
    print(paste('Yates Corrected X-squared Value: ', sprintf("%.5f", round(chi.yates, 5)), sep = ""))

  print(paste('Critical X-squared Value: ', sprintf("%.5f", round(chi_test, 5)), sep = ""))

  print(paste("Degree(s) of freedom: ", df, sep = ""))

  print(paste("Confidence Coefficient: ", alpha, sep = ""))

  print(paste("P-value: ", p_value, sep = ""))

  if(yates){
    values <- list("chi.yates.corrected" = chi.yates,
                   "chi.crit" = chi_test,
                   "df" = df,
                   "confidence.coefficient" = alpha,
                   "p.value" = p_value)
  }else{
    values <- list("chi.calc" = chi_calc,
                   "chi.crit" = chi_test,
                   "df" = df,
                   "confidence.coefficient" = alpha,
                   "p.value" = p_value)
  }

  cat("\n")
  if(verbose){
    print("---------------------------------------------")
    cat("\n")
    print("                  HYPOTHESIS                 ")
    cat("\n")
    print("---------------------------------------------")
    cat("\n")

    print("Null: The two classifications are independent")
    print("Alt: The two classifications are dependent")

    cat("\n")
    print("---------------------------------------------")
    cat("\n")
    print("                  CONCLUSION                 ")
    cat("\n")
    print("---------------------------------------------")
    cat("\n")

    if(p_value < alpha){
      print("There is sufficient evidence to reject the null hypothesis.")
      print("This means its plausible that the two classifications are dependent.")
    }else{
      print("There is not enough evidence to reject the null hypothesis.")
      print("This means its plausible that the two classifications are independent.")
    }
  }

  mc <- list("Data.Table" = table,
             "Contingency.Table" = mtab,
             "Probability.Table" = ptab,
             "test.vals" = values)

  invisible(mc)
}
