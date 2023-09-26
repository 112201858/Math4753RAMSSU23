# Defining a function to simulate hypergeometric distribution
# Function takes four parameters:
# 1) iter - the number of iterations (default initialized to 100)
# 2) N - the size of the population (default initialized to 20)
# 3) r - the number of successes in the population (default initialized to 10)
# 4) n - the number of trials (default initialized to 5)


#' Simulating Hypergeometric Distribution
#'
#' @param iter number of iterations
#' @param N size of population
#' @param r number of successes in population
#' @param n sample size
#'
#' @return a table and a bar plot
#' @export
#'
#' @examples
#' simHyper()
simHyper = function(iter = 100, N = 20, r = 10, n = 5){

  # Constructing a matrix object to hold the samples
  # Initialized with NA in each entry

  simMat = matrix(NA, nrow = n, ncol = iter, byrow = TRUE)

  # Constructing a vector object to hold the number of successes in each trial

  succ = c()

  # For loop
  # Fills each column with a new sample
  # Calculates the column sum

  for( i in 1:iter){
    simMat[,i] = sample(rep(c(1,0), c(r, N - r)), n, replace = FALSE)
    succ[i] = sum(simMat[,i])
  }

  #Making a table of the successes

  recordedSuccesses <- table(factor(succ, levels = 0:n))

  #Making a bar plot of the proportions

  barplot(recordedSuccesses/(iter),
          col = rainbow(n+1),
          main = "HYPERGEOMETRIC simulation",
          xlab = "Number of successes")


  # Storing proportions in a table

  simulatedProportions <- recordedSuccesses / iter
  simulatedProportions
}
