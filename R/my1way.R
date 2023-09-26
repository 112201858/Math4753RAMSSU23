#' my1way
#'
#' @param table a one way table
#' @param alpha significance value
#' @param adjusted Boolean for whether or not to Bonferonni correct
#' @param test Boolean for whether or not to perform Chi-squared test
#' @param verbose Boolean for whether or not to print extra information
#' @param visual Boolean for whether or not to show plots
#' @param screed Boolean for extra text dump of covariance matrices
#'
#' @return a table
#' @export
#'
#' @examples
#' my1way(table = table(c(1,2,1,2)))
my1way <- function(table,
                   alpha = .05,
                   adjusted = FALSE,
                   test = TRUE,
                   verbose = TRUE,
                   visual = TRUE,
                   screed = FALSE){

  ############### Making A TABLE ###############

  cat("\n")

  tab <- table

  # Instantiating a second table with margins
  mtab <- addmargins(tab)

  # Instantiating a third table with proportions
  ptab <- prop.table(tab)

  # Printing table to console
  if (verbose == TRUE){


    print('-----------------------------------------------')
    print("                      Table                    ")
    print('-----------------------------------------------')
    print(tab)
    cat("\n")

    print('-----------------------------------------------')
    print("           Category Counts W/ Margins          ")
    print('-----------------------------------------------')
    print(mtab)
    cat("\n")

  }

  # Displaying a bar plot of category counts
  if (visual == TRUE){

    bar_cat <- barplot(height = as.numeric(tab),
                       names = rownames(tab),
                       col = rgb(as.numeric(tab)/max(as.numeric(tab)),.3,0.3),
                       xlab = "Categories",
                       ylab = "Count",
                       main = "Category Counts",
                       ylim = c(0,max(as.numeric(tab))))

    text(bar_cat,
         4,
         as.numeric(tab),
         col = 'white',
         cex=2)

  }

  ############### DECLARING KEY VARIABLES ###############

  # A vector containing numeric data from each category
  category_counts <- as.numeric(tab)

  # The row sum from the table created above
  n <- as.numeric(mtab[length(tab)+1])

  # Vector containing the proportion estimator for each category
  phats <- category_counts/n

  # Z-score for the given alpha/2 value
  z <- qnorm(1-(alpha/2),lower.tail = TRUE)

  # Bonferroni Corrected Z Score
  Bz <- qnorm(1-alpha/(2*length(category_counts)),lower.tail = TRUE)

  # Confidence as percentage from given alpha
  confidence <- (1-alpha)

  ############### EXPECTATION & VARIANCE ###############

  # Variance of individual category count proportion
  var_phats <- phats*(1-phats) / n

  # Instantiating an NA matrix of dimension length(n_i) x length(n_i)
  expected_diff = matrix(data = NA,
                         nrow = length(category_counts),
                         ncol = length(category_counts))

  # Naming rows and columns
  rownames( expected_diff) = rownames(tab)
  colnames( expected_diff) = rownames(tab)

  # Filling matrix with expectation of difference in proportions
  for (i in 1:length(category_counts)){
    for (j in 1: length(category_counts)){
      expected_diff[i,j] = phats[i]-phats[j]
    }
  }

  ############### CONFIDENCE INTERVALS ###############

  if (verbose == TRUE){

    print('-----------------------------------------------')
    print("              Category Proportions             ")
    print('-----------------------------------------------')
    print(round(ptab, 4))
    cat("\n")

  }

  # Computing Intervals (Adjusted == FALSE)

  if (adjusted == FALSE) {

    lower = phats - z * sqrt(var_phats)
    upper = phats + z * sqrt(var_phats)

    CI_pi_mat = matrix(data = NA,
                       nrow = length(category_counts),
                       ncol = 2)

    rownames(CI_pi_mat) = rownames(tab)
    colnames(CI_pi_mat) = c('Lower', 'Upper')

    for (i in 0:length(category_counts)){

      CI_pi_mat[i,1] = lower[i]
      CI_pi_mat[i,2] = upper[i]

    }

    if (verbose == TRUE){

      print('-----------------------------------------------')
      print(paste("           ", confidence*100, '% CI For Cell Proportions         ', sep = ""))
      print('-----------------------------------------------')
      print(round(CI_pi_mat,4))
      cat("\n")

    }

    if (visual == TRUE){

      # Creating a bar plot for proportions
      l <- barplot(height = phats,
                   names = rownames(tab),
                   ylim = c(0,max(CI_pi_mat)+.05),
                   col = rgb(.3,0.3,phats/max(phats)),
                   ylab = "Probability",
                   xlab = "cell category",
                   main = paste("Category Proportions Showing ", (1-alpha)*100, "% CI", sep = ""))

      # Fixing a length in inches equivalent to the bar width given by window size
      length = .5 * par("pin")[1L]/diff(par("usr")[1:2])

      # Adding a line to indicate mean

      # Adding arrows to the plot
      arrows(l,upper,
             l,lower,
             angle=90,
             code=3,
             length = length,
             lwd = 2,
             col = 'red')

    }
  }
  # Computing Intervals (Adjusted == TRUE)

  if (adjusted == TRUE) {

    Blower = phats - Bz * sqrt(var_phats)
    Bupper = phats + Bz * sqrt(var_phats)

    CI_pi_mat_B = matrix(data = NA,
                         nrow = length(category_counts),
                         ncol = 2)

    rownames(CI_pi_mat_B) = rownames(tab)
    colnames(CI_pi_mat_B) = c('Lower', 'Upper')

    for (i in 0:length(category_counts)){

      CI_pi_mat_B[i,1] = Blower[i]
      CI_pi_mat_B[i,2] = Bupper[i]

    }

    if (verbose == TRUE){

      print('-----------------------------------------------')
      print(paste("           ", confidence*100, '% CI For Cell Proportions         ', sep = ""))
      print(paste('     Bonferonni Adjusted For ', length(category_counts), ' Hypothesis      ', sep =""))
      print('-----------------------------------------------')
      print(round(CI_pi_mat_B,4))
      cat("\n")

    }

    if (visual == TRUE){

      # Creating a bar plot for proportions
      l <- barplot(height = phats,
                   names = rownames(tab),
                   ylim = c(0,max(CI_pi_mat_B)+.05),
                   col = rgb(.3,0.3,phats/max(phats)),
                   ylab = "Probability",
                   xlab = "cell category",
                   main = paste("Category Proportions Showing Adjusted ", (1-alpha)*100, "% CI", sep = ""))

      # Fixing a length in inches equivalent to the bar width given by window size
      length = .5 * par("pin")[1L]/diff(par("usr")[1:2])

      # Adding a line to indicate mean

      # Adding arrows to the plot
      arrows(l,Bupper,
             l,Blower,
             angle=90,
             code=3,
             length = length,
             lwd = 2,
             col = 'red')

    }
  }

  ############### CHI-SQUARE TEST FOR UNIFORMITY ###############

  if (test == TRUE){

    p_null = 1 / length(category_counts)
    expected = (p_null * n)

    OECC_mat = matrix(data = NA,
                      nrow = 2,
                      ncol = length(category_counts))

    rownames(OECC_mat) = c('Observed', 'Expected')
    colnames(OECC_mat) = rownames(tab)

    for (i in 0:length(category_counts)){

      OECC_mat[1,i] = round(category_counts[i],0)
      OECC_mat[2,i] = expected
    }

    if (verbose == TRUE){

      print('-----------------------------------------------')
      print("      Observed & Expected Category Counts      ")
      print('-----------------------------------------------')
      print(round(OECC_mat,3))
      cat("\n")

    }

    # Performing Chi-Square Test

    print('-----------------------------------------------')
    print("         Chi-Square Test For Uniformity        ")
    print('-----------------------------------------------')

    df = length(category_counts)-1

    chisq_calc <- sum((category_counts-expected)^2/expected)
    print(paste("X-calc:", round(chisq_calc,3)))

    chisq_crit <- qchisq(1-alpha, length(category_counts)-1)
    print(paste("X-critical", round(chisq_crit,3)))

    p_value <- 1 - pchisq(q = chisq_calc,df = length(category_counts)-1,lower.tail = TRUE)
    print(paste("p-value:",round(p_value,6)))

    print(paste("alpha:", alpha))

    print(paste("df:", length(category_counts)-1))

    assumption = FALSE

    if (expected >= 5){
      print("Primary assumption holds.")
    } else{
      print("The data violates the primary assumption.")
    }
    cat("\n")

    if ((chisq_calc > chisq_crit) & (assumption == TRUE)){
      print("Kick The Null Out The Door!")
    } else{
      print("We cannot reject the null hypothesis.")
    }

    test_list <- list("chisq_calc" = chisq_calc,
                      'chisq_crit' = chisq_crit,
                      'p_value' = p_value,
                      'alpha' = alpha,
                      'df'= df)

    if(visual){
      #create density curve
      curve(dchisq(x, df = length(category_counts)-1), from = 0, to = 25,
            main = paste('Chi-Square Distribution (df =', length(category_counts)-1,')'),
            ylab = 'Density',
            xlab = 'It\'s A Rock & Roll Man!',
            lwd = 3)

      #create vector of x values
      x_vector <- seq(0, chisq_crit)
      x_vector2 <- seq(chisq_crit, 25)

      #create vector of chi-square density values
      p_vector <- dchisq(x_vector, length(category_counts)-1)
      p_vector2 <- dchisq(x_vector2, length(category_counts)-1)

      #fill in acceptance region
      polygon(c(x_vector, rev(x_vector)), c(p_vector, rep(0, length(p_vector))),
              col = adjustcolor('blue', alpha.f=0.3), border = NA)
      #fill in acceptance region
      polygon(c(x_vector2, rev(x_vector2)), c(p_vector2, rep(0, length(p_vector2))),
              col = adjustcolor('red', alpha.f=0.3), border = NA)
    }

  }

  ######## Computing Diff Intervals ########

  if (verbose){
    if(adjusted){
      print(cell_diff_proportion_CIs(table = tab,confidence = 1-alpha,adjusted = TRUE))
    }else{
      print(cell_diff_proportion_CIs(table = tab,confidence = 1-alpha,adjusted = FALSE))
    }
  }




  ############### COVARIANCE CATEGORY COUNTS ###############

  # Instantiating an NA matrix of dimension length(n_i) x length(n_i)
  cov_mat1 <- matrix(data = NA,
                     nrow = length(category_counts),
                     ncol = length(category_counts))

  # Naming rows and columns
  rownames(cov_mat1) = rownames(tab)
  colnames(cov_mat1) = rownames(tab)

  for (i in 1:length(category_counts)){
    for (j in 1: length(category_counts)){
      cov_mat1[i,j] = -1*n*phats[i]*phats[j]
    }
  }

  ############### COVARIANCE ESTIMATORS ###############

  # Instantiating an NA matrix of dimension length(n_i) x length(n_i)
  cov_mat2 <- matrix(data = NA,
                     nrow = length(category_counts),
                     ncol = length(category_counts))

  # Naming rows and columns
  rownames(cov_mat2) = rownames(tab)
  colnames(cov_mat2) = rownames(tab)

  for (i in 1:length(category_counts)){
    for (j in 1: length(category_counts)){
      cov_mat2[i,j] = (-1*phats[i]*phats[j])/n
    }
  }

  ############### screed  ###############

  if (screed == TRUE){

    print("*****************************************")
    print(" Pair-wise covariance of category counts ")
    print("*****************************************")
    print(cov_mat1)
    cat("\n")

    print("*****************************************")
    print(" Pair-wise covariance of estimators ")
    print("*****************************************")
    print(round(cov_mat2,6))
    cat("\n")

  }

  if (adjusted == FALSE){
    if (test == TRUE){
      mow <- list("Table" = tab,
                  "Mtable" = mtab,
                  "Ptable" = ptab,
                  "CIs" = CI_pi_mat,
                  "UniformityTest" = test_list)
    }
    if (test == FALSE){
      mow <- list("Table" = tab,
                  "Mtable" = mtab,
                  "Ptable" = ptab,
                  "CIs" = CI_pi_mat)
    }
  }

  if (adjusted == TRUE){
    if (test == TRUE){
      mow <- list("Table" = tab,
                  "Mtable" = mtab,
                  "Ptable" = ptab,
                  "CIs" = CI_pi_mat_B,
                  "UniformityTest" = test_list)
    }
    if (test == FALSE){
      mow <- list("Table" = tab,
                  "Mtable" = mtab,
                  "Ptable" = ptab,
                  "CIs" = CI_pi_mat_B)
    }
  }


  invisible(mow)
}
