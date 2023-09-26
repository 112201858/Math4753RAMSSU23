#' bigfish
#'
#' @param table an rXc contingency table (for 2X2 use twofish)
#' @param alpha A confidence coefficient (default = .05)
#' @param sim_lower A lower bound for number of iterations used to simulate p-value (default = 1000)
#' @param sim_upper An upper bound for number of iterations used to simulate p-value (default = 2000)
#' @param verbose A Boolean for whether or not to print extra information (default = TRUE)
#' @param visual A Boolean for whether or not to construct plots (default = TRUE)
#'
#' @return an invisible list called 'bf' that contains various tables,
#' a vector of simulated p-values,
#' the mean p-value from simulation
#' @export
#' @importFrom stats fisher.test
#' @importFrom grDevices cm.colors
#' @importFrom graphics mosaicplot
#'
#' @examples
#' bigfish(table = table(c(1,2,3),c(1,2,3)))
bigfish <- function(table,
                    alpha = .05,
                    sim_lower = 1000,
                    sim_upper = 5000,
                    verbose = TRUE,
                    visual = TRUE){

  # Constructing relevant tables

  mtab <- addmargins(table)

  ptab <- addmargins(prop.table(table))

  # Simulating P-Values

  p_values <- c()

  for (i in sim_lower:sim_upper){

    result <- fisher.test(table, simulate.p.value = TRUE, B = i)$p.value
    p_values <- c(p_values, result)

  }

  # Computing the mean of the simulated P-values

  mean <- mean(p_values)

  # Printing Tables when verbose == TRUE

  if(verbose){

    print("-----------------------------------------------------------------")
    print("               Observed Counts for Contingency Table             ")
    print("-----------------------------------------------------------------")
    print(mtab)

    cat("\n")

    print("-----------------------------------------------------------------")
    print("               Probabilities for Contingency Table               ")
    print("-----------------------------------------------------------------")
    print(ptab)

    cat("\n")

  }


  # Plots

  if (visual){

    # Mosaic Plot of Data
    mosaicplot(table,
               main = "Mosaic plot",
               color = cm.colors(length(colnames(table))))

    # Plotting simulated values when visual == TRUE
    plot(sim_lower:sim_upper,
         p_values,
         type = 'l',
         pch=19,
         xlab ='Iterations',
         ylab = 'Simulated Value',
         col='cyan',
         main = paste("Simulated P Values(mean = ", round(mean,4), ")",sep=""))

    abline(h = mean,
           col='blue',
           lw = 3)

  }

  print("-----------------------------------------------------------------")
  print("                             Summary                             ")
  print("-----------------------------------------------------------------")

  cat("\n")

  print("Null Hypothesis:")
  print("The two directions of data classification are dependent.")

  cat("\n")

  print("Alternative Hypothesis:")
  print("The two directions of data classification are independent.")

  cat("\n")

  print("Simulated Fisher's Exact P-Value:")
  print(round(mean,4))

  cat("\n")

  if (mean < alpha){
    print(paste("Given the confidence coefficient provided(",alpha,"), we cannot accept the null hypothesis.", sep=""))
    print("It is plausible that the two directions of data classification are dependent.")
  }else{
    print(paste("Given the confidence coefficient provided(",alpha,"), we cannot reject the null hypothesis.", sep=""))
    print("It is plausible that the two directions of data classification are independent.")
  }

  cat("\n")


  bf <- list("Data.Table" = table,
             "Contingency.Table" = mtab,
             "Probability.Table" = ptab,
             "p.values" = p_values,
             "p.mean" = mean)

  invisible(bf)

}
