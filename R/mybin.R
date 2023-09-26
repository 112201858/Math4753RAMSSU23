# Defining a function to simulate binomial distribution
# Function takes three parameters:
# 1) iter - the number of iterations (default initialized to 100)
# 2) n - the number of trials (default initialized to 10)
# 3) p - the probability of a success on a single trial (default initialized to .5)

#' Binomial Simulation
#'
#' @param iter number of iterations
#' @param n sample size
#' @param p probability of success
#'
#' @return a table and a bar plot
#' @export
#'
#' @examples
#' mybin()
mybin = function(iter = 100, n = 10, p = 0.5){

  # Constructing a matrix object to hold the samples
  # Initialized with NA in each entry

  simMat = matrix(NA, nrow = n, ncol = iter, byrow = TRUE)

  # Constructing a vector object to hold the number of successes in each trial

  succ = c()

  # For loop
  # Fills each column with a new sample
  # Calculates the column sum

  for( i in 1:iter){
    simMat[,i] = sample(c(1,0), n, replace = TRUE, prob = c(p,1-p))
    succ[i] = sum(simMat[,i])
  }

  #Making a table of the successes

  recordedSuccesses <- table(factor(succ, levels = 0:n))

  #Making a bar plot of the proportions

  barplot(recordedSuccesses/(iter),
          col = rainbow(n+1),
          main = "Binomial simulation",
          xlab = "Number of successes")

  # Storing proportions in a table

  simulatedProportions <- recordedSuccesses / iter

  simulatedProportions
}


