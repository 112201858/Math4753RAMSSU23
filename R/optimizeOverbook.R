#' optimizeOverbook
#'
#' ── R CMD check results ────────────────────────────────────────────── MATH4753RAMSSU23 0.1.0 ────
#' Duration: 15s
#'
#' 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
#'
#' R CMD check succeeded
#'
#' @param N The total number of seats on a flight
#' @param gamma The pain tolerance
#' @param p The probability that a ticket holder shows up
#'
#' @return A named list, a graph of discrete case, a graph of continuous approximation
#' @export
#'
#' @examples
#' optimizeOverbook()
optimizeOverbook <- function(N = 200, gamma = .02, p = .95){

  ############ Upper Bound ############

  for (i in 1:100000) {
    targ <- 1 - .99
    if (pbinom(N,N+i,p) <= targ) {
      index <- i
      upper <- N + index
      break
    }
  }

  ############ Building Layout ############

  layout(matrix(1:2, nrow = 2, ncol = 1))

  ############ DISCRETE ############

  # Creating a vector and storing probability of overbooking on N:upper trials
  dProb <- pbinom(N,N:upper,p)

  # Creating a vector to store the outputs of the objective function (discrete case)
  objective1 <- (1 - gamma - dProb)

  # Finding the element in |objective1| closest to zero
  n1 <-  which.min(abs(objective1))

  # Storing the output of objective1 @ minimum
  val1 <- objective1[n1]

  # Storing the optimal number of tickets sold
  result1 <- n1+N-1

  # Constructing an object to hold the title of the discrete plot
  title1 <- paste("Objective Vs n to find optimal tickets sold",
                  "\n(", result1, ") gamma =", gamma,"N =",N,"discrete")

  # Plotting 'objective1' vector
  plot(N:upper,
       objective1,
       col = ifelse(objective1 == val1, "red", "black"),
       pch = 21,
       bg = ifelse(objective1 == val1, "red", "blue"),
       type = "b",
       ylab = "Objective",
       xlab = "n",
       main = title1)

  # Highlighting the optimal value
  abline(h = 0, v = result1 , col = "red", lwd=2.0)

  # Adding an embellished point to aid the viewer
  points(x=result1, y=0, col = "red", pch = 19, cex = 1.25)


  ############ Continuous ############

  # Creating a sequence to use as the domain
  domain <- seq(N, N+upper, length.out = 100000)

  # Calculating the probabilities using domain and storing them in a vector
  cProb <- pnorm(N+.5, domain*p, sqrt(domain*p*(1-p)))

  # Creating a vector to store the outputs of the objective function (continuous case)
  obj2 <- 1 - gamma - cProb

  # Finding the element in |objective2| closest to zero
  n2 <- which.min(abs(obj2))

  # Storing the number of tickets sold
  result2 <- domain[n2]

  # Creating objective2 function for curve() expr parameter
  objective2 <- function(x){
    1 - gamma - pnorm(q = N+.5, mean = x*p, sd = sqrt(x*p*(1-p)))
  }

  # Constructing an object to hold the title of the continuous curve
  title2 <- paste("Objective Vs n to find optimal tickets sold",
                  "\n(", result2, ") gamma =", gamma,"N =",N,"continuous")

  # Using objective2 (function) to generate a curve
  curve(expr = objective2(x),
        from = N,
        to = upper,
        n = 100000,
        ylab = 'objective',
        xlab = 'n',
        main = title2)

  # Highlighting the optimal value
  abline(h = 0, v = result2 , col = "red")

  # Adding an embellished point to aid the viewer
  points(x = result2, y = 0, col = "red", pch = 19, cex = 1.25)


  ############ Named List ############

  # Creating the named list
  output <- list(nd = result1,
            nc = result2,
            N = N,
            p = p,
            gamma = gamma)

  # Printing List
  print(output)

}

