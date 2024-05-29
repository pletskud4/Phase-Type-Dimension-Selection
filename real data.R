#SAMLET KODEBIBLIOTEK
#vi skal have følgende funktioner: TVR, Q-matrix, confidence bands, check confidence, first guess dimensions, add/remove dim
#alt skal være funktion af datasæt, dimensioner, konfidens og change kun

#standard ting vi skal have med
library(matrixdist)
set.seed(999)
confidence <- 0.05
change <- 0.001
directory<- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Real-Life data"
csv_files <- list.files(directory, pattern = "*.csv", full.names = TRUE)
dataset <- read.csv(csv_files[1])[,2]
largest_values <- sort(read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")[,2])
csv_files
#FIRST GUESS DIMENSIONS
Firstguess <- function(dataset, confidence) {
  empirical <- ecdf(dataset)
  mean <- sum(dataset)/length(dataset)
  rate <- 1/mean
  probs <- empirical(sort(dataset))
  cdfprobs <- pexp(sort(dataset),rate)
  value <- sum((probs-cdfprobs)^2)/sum((probs-mean(probs))^2)
  dim <- ceiling((sqrt(((confidence)^0.5)*(value)^(0.5))*50))
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
    for(j in 1:(limit)){
      if(T_matrix[i,j] == 0){
        T_matrix[i, j] <- (n * Q_values[limit]) / (limit-1)
      }
    }
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
  
  alpha <- c(0.99, rep(0.01*(1/(limit-1)), limit - 1))
  ph2 <- ph(structure = "general", dimension = limit)
  ph2@pars$alpha <- alpha
  ph2@pars$S <- T_matrix
  
  return(ph2)
}

Remove <- function(dist, n) { 
  GM <- solve((-dist@pars$S)) 
  CV <- dist@pars$alpha %*% GM 
  
  rewards <- rep(1, length(dist@pars$alpha))
  sorted_indices <- order(CV)
  smallest_elements <- sorted_indices[1:n]
  rewards[smallest_elements] <- 0
  
  return(TVR(dist, rewards))
}


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

Check_confidence <- function(dist, lowerbound, upperbound,dataset){
  sequence <- seq(min(dataset), max(dataset), length.out = 10000)
  cdfs <- cdf(dist,sequence)
  return(sum(cdfs <= upperbound & cdfs >= lowerbound))
}

Total_model <- function(dataset, confidence, change){
  #immediate results
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
    Evaluation <- Check_confidence(Fit1, Confidence_bands[,1],Confidence_bands[,2],dataset)
    print(Evaluation)
    if (Evaluation == 10000 && length(Fit1@pars$alpha) == 1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha)
      Fit2 <- Fit1
      break
    } else if (Evaluation == 10000 && status == -1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha)
      Fit2 <- Fit1
      break
    } else if (Evaluation == 10000 && (status == 1 || status == 0)){
      Fit2 <- Fit1
      Fit1 <- Remove(Fit1, 1)
      Fit1 <- fit(Fit1, y = dataset, stepsEM = 0)
      Loglik <- logLik(Fit1)
      status <- 1
    } else if (Evaluation < 10000 && status == 1){
      End_time <- Sys.time()
      Dimensions <- length(Fit1@pars$alpha+1)
      break
    } else if (Evaluation < 10000 && (status == -1 || status == 0)){
      Fit1 <- Q(length(Fit1@pars$alpha)+1,dataset)
      Fit1 <- fit(Fit1, y = dataset, stepsEM = 0)
      Loglik <- logLik(Fit1)
      status <- -1
    } 
  }
  Result_vector[1] <- Initial_dim
  Result_vector[2] <- Dimensions
  Result_vector[3] <- End_time - Start_time2
  Result_vector[4] <- End_time - Start_time
  return(list(Result_vector = Result_vector, FirstFit = FirstFit, Fit2 = Fit2))
}



#BMI data først 
dataset2 <- (log((na.omit(dataset)+0.001)/min(na.omit(dataset))))
empirical <- ecdf(dataset2)
#HISTOGRAM
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "BMI hist.png", width = 10, height = 8, units = "in", res = 600)
par(mfrow = c(1, 2))
hist(dataset, main = "BMI Data", xlab = "Values", col = "darkgreen", prob = T, breaks = 10)
hist(dataset2, main = "Transformed BMI Data", xlab = "Values", col = "darkgreen", prob = T, breaks = 10)
par(mfrow = c(1, 1))
dev.off()
#sample
sample_indices <- sample(seq_along(dataset2), size = 20)
BMIsample <- dataset2[sample_indices]
BMInonsample <- dataset2[-sample_indices]
empirical2 <- ecdf(BMInonsample)
empiricalnonsample <- empirical2(sort(BMInonsample))
sequence <- seq(min(dataset2), max(dataset2), length.out = 10000)
BMI999fit <- Total_model(BMIsample, 0.001,change)
BMI95fit <- Total_model(BMIsample, 0.05,change) 
BMI999fitted <- cdf(BMI999fit$Fit2, sequence)
BMI95fitted <- cdf(BMI95fit$Fit2, sequence)
size999 <- sqrt(1 /  length(BMIsample))*quantile(largest_values, probs = 1-0.001)
size95 <- sqrt(1 /  length(BMIsample))*quantile(largest_values, probs = 1-0.05)
conf99up <- pmin(BMI999fitted + size999,1)
conf99down <- pmax(BMI999fitted - size999,0)
conf95up <- pmin(BMI95fitted + size95,1)
conf95down <- pmax(BMI95fitted - size95,0)
#sampleplots
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "BMI95sample fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), max(dataset2)), ylim=c(0, 1))
lines(sort(BMInonsample),empiricalnonsample,col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
#full
BMI999fit <- Total_model(dataset2, 0.001,change)
BMI95fit <- Total_model(dataset2, 0.05,change) 
BMI999fitted <- cdf(BMI999fit$Fit2, sequence)
BMI95fitted <- cdf(BMI95fit$Fit2, sequence)
size999 <- sqrt(1 /  length(dataset2))*quantile(largest_values, probs = 1-0.001)
size95 <- sqrt(1 /  length(dataset2))*quantile(largest_values, probs = 1-0.05)
conf99up <- pmin(empirical(sequence) + size999,1)
conf99down <- pmax(empirical(sequence) - size999,0)
conf95up <- pmin(empirical(sequence) + size95,1)
conf95down <- pmax(empirical(sequence) - size95,0)
#full plots
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "BMI full fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), max(dataset2)), ylim=c(0, 1))
lines(sort(dataset2),empirical(sort(dataset2)),col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()























Firstguess(dataset3, confidence)
Initial_dist <- Q(30, dataset3)
Fit1 <- fit(Fit1, y = dataset3, stepsEM = 500, every =1)
sequence <- seq(min(dataset3), max(dataset3), length.out = 10000)
fitted <- cdf(Fit1, sequence)
bounds <- Calc_confidence(dataset3, 0.05)
Check_confidence(Fit1,bounds[,1], bounds[,2],dataset3)

plot(sort(dataset3),empirical(sort(dataset3)), col = "yellow")
points(sequence,fitted,cex=0.2)
lines(sequence,bounds[,1],cex=0.2)
lines(sequence,bounds[,2],cex=0.2)
0.5^0.5
0.25^0.5
mean(log((na.omit(datasetclean)+0.001)/min(na.omit(datasetclean))))










#dataset cleaned
datasetclean <- numeric(length = length(dataset)-1)
for (i in 2:(length(dataset))){
  datasetclean[i-1] <- dataset[i]/dataset[i-1]
}
#hernæst currency data først 
dataset2 <- (log((na.omit(datasetclean)+0.001)/min(na.omit(datasetclean))))^5
empirical <- ecdf(dataset2)
#HISTOGRAM
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Currency Data histogram.png", width = 10, height = 8, units = "in", res = 600)
par(mfrow = c(1, 2))
hist(dataset, main = "Currency Data", xlab = "Values", col = "darkgreen", prob = T, breaks = 50)
hist(dataset2, main = "Transformed Currency Data", xlab = "Values", col = "darkgreen", prob = T, breaks = 50)
par(mfrow = c(1, 1))
dev.off()
#sample
sample_indices <- sample(seq_along(dataset2), size = 200)
BMIsample <- dataset2[sample_indices]
BMInonsample <- dataset2[-sample_indices]
empirical2 <- ecdf(BMInonsample)
empiricalnonsample <- empirical2(sort(BMInonsample))
sequence <- seq(min(dataset2), max(dataset2), length.out = 10000)
BMI999fit <- Total_model(BMIsample, 0.001,change)
BMI95fit <- Total_model(BMIsample, 0.05,change) 
BMI999fitted <- cdf(BMI999fit$Fit2, sequence)
BMI95fitted <- cdf(BMI95fit$Fit2, sequence)
size999 <- sqrt(1 /  length(BMIsample))*quantile(largest_values, probs = 1-0.001)
size95 <- sqrt(1 /  length(BMIsample))*quantile(largest_values, probs = 1-0.05)
conf99up <- pmin(BMI999fitted + size999,1)
conf99down <- pmax(BMI999fitted - size999,0)
conf95up <- pmin(BMI95fitted + size95,1)
conf95down <- pmax(BMI95fitted - size95,0)
#sampleplots
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "currency95sample fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), max(dataset2)), ylim=c(0, 1))
lines(sort(BMInonsample),empiricalnonsample,col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
#zoomer
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "currency95sample fit zoom.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), 0.00004), ylim=c(0, 1))
lines(sort(BMInonsample),empiricalnonsample,col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
#full
BMI999fit <- Total_model(dataset2, 0.001,change)
BMI95fit <- Total_model(dataset2, 0.05,change)  #denne eer i gang
BMI999fitted <- cdf(BMI999fit$Fit2, sequence)
BMI95fitted <- cdf(BMI95fit$Fit2, sequence)
size999 <- sqrt(1 /  length(dataset2))*quantile(largest_values, probs = 1-0.001)
size95 <- sqrt(1 /  length(dataset2))*quantile(largest_values, probs = 1-0.05)
conf99up <- pmin(empirical(sequence) + size999,1)
conf99down <- pmax(empirical(sequence) - size999,0)
conf95up <- pmin(empirical(sequence) + size95,1)
conf95down <- pmax(empirical(sequence) - size95,0)
#full plots
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Currency full fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), max(dataset2)), ylim=c(0, 1))
lines(sort(dataset2),empirical(sort(dataset2)),col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
#zoom
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Currency full fit zoom.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset2), 0.00005), ylim=c(0, 1))
lines(sort(dataset2),empirical(sort(dataset2)),col = "black", type = "s",lwd = 2)
lines(sequence,  BMI999fitted, col = "darkgreen", type = "s",lwd = 2)
lines(sequence, conf99up, col = "darkgreen",type = "s", lty = 2)
lines(sequence, conf99down, col = "darkgreen",type = "s", lty = 2)
lines(sequence,  BMI95fitted, col = "lightgreen", type = "s",lwd = 2)
lines(sequence, conf95up, col = "lightgreen",type = "s", lty = 2)
lines(sequence, conf95down, col = "lightgreen",type = "s", lty = 2)
legend("bottomright", legend = c("Empirical of non-sampled data", "Fit and confidence 99.9%", "Fit and confidence 95%"),
       col = c("black", "darkgreen", "lightgreen"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()

