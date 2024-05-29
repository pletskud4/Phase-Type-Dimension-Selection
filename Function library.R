#SAMLET KODEBIBLIOTEK
#vi skal have følgende funktioner: TVR, Q-matrix, confidence bands, check confidence, first guess dimensions, add/remove dim
#alt skal være funktion af datasæt, dimensioner, konfidens og change kun


#standard ting vi skal have med
library(matrixdist)
set.seed(999)
confidence <- 0.05
change <- 0.01
directory<- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets"
csv_files <- list.files(directory, pattern = "*.csv", full.names = TRUE)
dataset <- as.numeric(read.csv(csv_files[2])[,2])
largest_values <- sort(read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")[,2])

#FIRST GUESS DIMENSIONS
Firstguess <- function(dataset, confidence) {
  empirical <- ecdf(dataset)
  mean <- sum(dataset)/length(dataset)
  rate <- 1/mean
  probs <- empirical(sort(dataset))
  cdfprobs <- pexp(sort(dataset),rate)
  value <- sum((probs-cdfprobs)^2)/sum((probs-mean(probs))^2)
  dim <- ceiling((sqrt(((confidence)^0.1)*(value)^0.75)*50))
  return(dim)
}

#FIRST GUESS DIMENSIONS
Q <- function(dim, dataset) { 
  N <- length(dataset)
  n <- max(ceiling(dim / max(dataset)) - 1, 1)
  limit <- dim 
  
  ordered_data <- dataset[order(dataset)]
  Q_values <- numeric(length = limit)
  
  for (x in 1:limit) {
    pois_dens <- dpois(x - 1, lambda = n * ordered_data)
    pois_dens_times_frac <- pois_dens * log((N + 1 - (1:length(pois_dens) - 1)) / (N - (1:length(pois_dens) - 1)))
    Q_values[x] <- abs(1 - sum(pois_dens_times_frac[-N]))
  }
  
  Q_values[Q_values > 1] <- 1
  
  T_matrix <- matrix(0, nrow = limit, ncol = limit)
  diag(T_matrix) <- -n
  
  for (i in 1:(limit - 1)) {
    T_matrix[i, i + 1] <- n * Q_values[i]
  }
  
  for (i in 1:(limit - 1)) {
    non_diag_indices <- setdiff(1:limit, i)
    T_matrix[i, non_diag_indices] <- (n * Q_values[i]) / (limit - 1)
  }
  
  for (j in 1:limit) {
    if (T_matrix[limit, j] == 0) {
      T_matrix[limit, j] <- (n * Q_values[limit]) / limit
    }
  }
  
  for (i in 1:nrow(T_matrix)) {
    while (sum(T_matrix[i, ]) >= 0) {
      for (j in 1:ncol(T_matrix)) {
        if (i != j) {
          T_matrix[i, j] <- T_matrix[i, j] * 0.99
        }
      }
    }
  }
  
  alpha <- c(1, rep(0, limit - 1))
  ph2 <- ph(structure = "general", dimension = limit)
  ph2@pars$alpha <- alpha
  ph2@pars$S <- T_matrix
  
  return(ph2)
}

#code to remove state
Remove <- function(dist, n) { 
  GM <- solve((-dist@pars$S)) 
  CV <- dist@pars$alpha %*% GM 
  
  rewards <- rep(1, length(dist@pars$alpha))
  sorted_indices <- order(CV)
  smallest_elements <- sorted_indices[1:n]
  rewards[smallest_elements] <- 0
  
  return(TVR(dist, rewards))
}


#calculate confidence
Calc_confidence <- function(dataset, confidence){
  empirical1 <- ecdf(dataset)
  sequence <- seq(min(dataset), max(dataset), length.out = 10000)
  xs <- empirical1(sequence)
  size <- sqrt(1 /  length(dataset))*quantile(largest_values, probs = 1-confidence)
  ks_bound <- matrix(nrow = 10000, ncol = 2)
  ks_bound[,1] <- pmax(xs - size, 0)
  ks_bound[,2] <- pmin(xs + size, 1)
  return(ks_bound)
}

#check confidence
Check_confidence <- function(dist, lowerbound, upperbound){
  sequence <- seq(min(dataset), max(dataset), length.out = 10000)
  cdfs <- cdf(dist,sequence)
  return(sum(cdfs <= upperbound & cdfs >= lowerbound))
}


#Total model
Total_model <- function(dataset, confidence, change){
  Start_time <- Sys.time()
  Confidence_bands <- Calc_confidence(dataset, confidence)
  Initial_dim <- Firstguess(dataset, confidence)
  Result_vector <- numeric(4)
  names(Result_vector) <- c("Initial dimension", "Final dimension", "Time in last dimension", "Total time")
  Initial_dist <- Q(Initial_dim, dataset)
  FirstFit <- fit(Initial_dist, y = dataset, stepsEM = 1)
  Fit1 <- fit(Initial_dist, y = dataset, stepsEM = 1)
  Loglik <- logLik(Fit1)
  status <- 0
  for (i in 1:50){
    Start_time2 <- Sys.time()
    if (length(Fit1@pars$alpha) >= 50){
      break
    }
    for (iter in 1:50000){
      Fit1 <- fit(Fit1, y = dataset, stepsEM = 1)
      if (abs(Loglik - logLik(Fit1))<change){
        break
      }
      Loglik <- logLik(Fit1)
      print(length(Fit1@pars$alpha))
    }
    Evaluation <- Check_confidence(Fit1, Confidence_bands[,1],Confidence_bands[,2])
    print(Evaluation)
    if (Evaluation == 10000 && length(Fit1@pars$alpha) == 1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha)
      break
    } else if (Evaluation == 10000 && status == -1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha)
      break
    } else if (Evaluation == 10000 && (status == 1 || status == 0)){
      Fit1 <- Remove(Fit1, 1)
      Fit1 <- fit(Fit1, y = dataset, stepsEM = 0)
      Loglik <- logLik(Fit1)
      status <- 1
    } else if (Evaluation < 10000 && status == 1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha)
      break
    } else if (Evaluation < 10000 && (status == -1 || status == 0)){
      Fit1 <- Q(length(Fit1@pars$alpha)+1,dataset)
      Fit1 <- fit(Fit1, y = dataset, stepsEM = 0)
      Loglik <- logLik(Fit1)
      status <- -1
    } 
    Fit2 <- Fit1
    
  }
  Result_vector[1] <- Initial_dim
  Result_vector[2] <- Dimensions
  Result_vector[3] <- End_time - Start_time2
  Result_vector[4] <- End_time - Start_time
  return(list(Result_vector = Result_vector, FirstFit = FirstFit, Fit2 = Fit2))
}

#99.9%
PH20 <- Total_model(dataset, confidence, change)
Multimodal <- Total_model(dataset, confidence, change)

#95%
PH202 <- Total_model(dataset, confidence, change)
Multimodal2 <- Total_model(dataset, confidence, change)


