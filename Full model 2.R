#SAMLET KODEBIBLIOTEK
#vi skal have følgende funktioner: TVR, Q-matrix, confidence bands, check confidence, first guess dimensions, add/remove dim
#alt skal være funktion af datasæt, dimensioner, konfidens og change kun

#standard ting vi skal have med
library(matrixdist)
set.seed(999)
confidence <- 0.05
change <- 0.001
directory<- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets"
csv_files <- list.files(directory, pattern = "*.csv", full.names = TRUE)
dataset <- as.numeric(read.csv(csv_files[10])[,2])
largest_values <- sort(read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")[,2])

#FIRST GUESS DIMENSIONS
Firstguess <- function(dataset, confidence) {
  empirical <- ecdf(dataset)
  mean <- sum(dataset)/length(dataset)
  rate <- 1/mean
  probs <- empirical(sort(dataset))
  cdfprobs <- pexp(sort(dataset),rate)
  value <- sum((probs-cdfprobs)^2)/sum((probs-mean(probs))^2)
  dim <- ceiling((sqrt(((confidence)^0.15)*(value)^(1/3))*50))
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
Check_confidence <- function(dist, lowerbound, upperbound,dataset){
  sequence <- seq(min(dataset), max(dataset), length.out = 10000)
  cdfs <- cdf(dist,sequence)
  return(sum(cdfs <= upperbound & cdfs >= lowerbound))
}


#Total model
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
  #loop until succesfull
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
PH20 <- Total_model(dataset, 0.001, change)
Multimodal <- Total_model(dataset, 0.001, change)

#95%
PH202 <- Total_model(dataset, 0.05, change)
Multimodal2 <- Total_model(dataset, 0.05, change)

#plots for PH
Bounds95 <- Calc_confidence(dataset, 0.05)
Bounds999 <- Calc_confidence(dataset, 0.001)
empirical <- ecdf(dataset)
sequence <- seq(min(dataset), max(dataset),length.out = 10000)
ph2095first <- cdf(PH202$FirstFit, sequence)
ph2099first <- cdf(PH20$FirstFit, sequence)
ph2095final <- cdf(PH202$Fit2, sequence)
ph2099final <- cdf(PH20$Fit2, sequence)
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "PH20 full automatic fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0, max(dataset)), ylim=c(0, 1))
lines(sequence,  empirical(sequence), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, Bounds95[,1], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds95[,2], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds999[,1], col = "green",type = "s", lty = 2)
lines(sequence, Bounds999[,2], col = "green",type = "s", lty = 2)
lines(sequence, ph2095final, col = "blue")
lines(sequence, ph2099final, col = "red")
lines(sequence, ph2099first, col = "purple")
lines(sequence, ph2095first, col = "orange")
legend("bottomright", legend = c("Final fit 95%", "Final fit 99.9%", "Initial guess 95%", "Initial guess 99.9%"),
       col = c("blue", "red", "orange","purple"), lwd = 2, 
       bty = "n", cex = 1.5)
legend("topleft", legend = c("95% confidence", "99.9% confidence"),
       col = c("lightgreen", "green"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
png(filename = "PH20 full automatic fit zoom.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0.8, 1), ylim=c(0.3, 0.4))
lines(sort(dataset),  empirical(sort(dataset)), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, Bounds95[,1], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds95[,2], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds999[,1], col = "green",type = "s", lty = 2)
lines(sequence, Bounds999[,2], col = "green",type = "s", lty = 2)
lines(sequence, ph2095final, col = "blue")
lines(sequence, ph2099final, col = "red")
lines(sequence, ph2099first, col = "purple")
lines(sequence, ph2095first, col = "orange")
legend("bottomright", legend = c("Final fit 95%", "Final fit 99.9%", "Initial guess 95%", "Initial guess 99.9%"),
       col = c("blue", "red", "orange","purple"), lwd = 2, 
       bty = "n", cex = 1.5)
legend("topleft", legend = c("95% confidence", "99.9% confidence"),
       col = c("lightgreen", "green"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()

#plots for Multimodal
Bounds95 <- Calc_confidence(dataset, 0.05)
Bounds999 <- Calc_confidence(dataset, 0.001)
empirical <- ecdf(dataset)
sequence <- seq(min(dataset), max(dataset),length.out = 10000)
Multimodal99first <- cdf(Multimodal$FirstFit, sequence)
Multimodal95final <- cdf(Multimodal2$Fit2, sequence)
Multimodal99final <- cdf(Multimodal$Fit2, sequence)
Multimodal95first <- cdf(Multimodal2$FirstFit, sequence)
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Multimodal full automatic fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0, max(dataset)), ylim=c(0, 1))
lines(sort(dataset),  empirical(sort(dataset)), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, Bounds95[,1], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds95[,2], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds999[,1], col = "green",type = "s", lty = 2)
lines(sequence, Bounds999[,2], col = "green",type = "s", lty = 2)
lines(sequence, Multimodal95final, col = "blue")
lines(sequence, Multimodal99final, col = "red")
lines(sequence, Multimodal99first, col = "purple")
lines(sequence, Multimodal95first, col = "orange")
legend("bottomright", legend = c("Final fit 95%", "Final fit 99.9%", "Initial guess 95%", "Initial guess 99.9%"),
       col = c("blue", "red", "orange","purple"), lwd = 2, 
       bty = "n", cex = 1.5)
legend("topleft", legend = c("95% confidence", "99.9% confidence"),
       col = c("lightgreen", "green"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()




















#endelig test
empirical <- ecdf(dataset)
dataset2 <- log((dataset+0.001)/min(dataset))
dataset3 <- log(dataset+1)
bounds <- Calc_confidence(dataset,0.001)
bounds2 <- Calc_confidence(dataset2,0.001)
bounds3 <- Calc_confidence(dataset3,0.001)
#fitter
Firstguess(dataset, confidence)
Firstguess(dataset2, confidence)
Firstguess(dataset3, confidence)
sequence <- seq(min(dataset), max(dataset), length.out = 10000)
Initial_dist <- Q(50,dataset)
Fit <- fit(Initial_dist, y = dataset, stepsEM = 800, every = 1)
datasetdone <- cdf(Fit,sequence)
sequence2 <- seq(min(dataset2), max(dataset2), length.out = 10000)
Initial_dist2 <- Q(36,dataset2)
Fit2 <- fit(Initial_dist2, y = dataset2, stepsEM = 800, every = 1)
dataset2done <- cdf(Fit2,sequence2)
sequence3 <- seq(min(dataset3), max(dataset3), length.out = 10000)
Initial_dist3 <- Q(41,dataset3)
Fit3 <- fit(Initial_dist3, y = dataset3, stepsEM = 800, every = 1)
dataset3done <- cdf(Fit3,sequence3)
#checker resultat
sum(dataset2done >= bounds2[,1] & dataset2done <= bounds2[,2])
sum(dataset3done >= bounds3[,1] & dataset3done <= bounds3[,2])
#plotter
plot(dataset2, empirical(dataset2), col ="yellow")
points(sequence2, dataset2done, cex = 0.2)
lines(sequence2,bounds2[,1])
lines(sequence2,bounds2[,2])
plot(dataset3, empirical(dataset3), col ="yellow")
points(sequence3, dataset3done, cex = 0.2)
lines(sequence3,bounds3[,1])
lines(sequence3,bounds3[,2])





#plot simpel lognormal (uden transformations)
Bounds95 <- Calc_confidence(dataset, 0.05)
Bounds999 <- Calc_confidence(dataset, 0.001)
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Log-Normal 50d.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(min(dataset), max(dataset)), ylim=c(0, 1))
lines(sort(dataset),  empirical(sort(dataset)), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, Bounds95[,1], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds95[,2], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds999[,1], col = "green",type = "s", lty = 2)
lines(sequence, Bounds999[,2], col = "green",type = "s", lty = 2)
lines(sequence, datasetdone, col = "red")
legend("bottomright", legend = c("Log-Normal final fit 50d"),
       col = c("red"), lwd = 2, 
       bty = "n", cex = 1.5)
legend("topleft", legend = c("95% confidence", "99.9% confidence"),
       col = c("lightgreen", "green"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()


#plots for transformeret lognormal
#99.9%
LogTrans <- Total_model(dataset2, 0.001, change)
#95%
LogTrans2 <- Total_model(dataset2, 0.05, change)
#plots for Multimodal
Bounds95 <- Calc_confidence(dataset2, 0.05)
Bounds999 <- Calc_confidence(dataset2, 0.001)
empirical <- ecdf(dataset2)
sequence <- seq(min(dataset2), max(dataset2),length.out = 10000)
LogTrans99first <- cdf(LogTrans$FirstFit, sequence)
LogTrans95final <- cdf(LogTrans2$Fit2, sequence)
LogTrans99final <- cdf(LogTrans$Fit2, sequence)
LogTrans95first <- cdf(LogTrans$FirstFit, sequence)
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Log transformed full automatic fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0, max(dataset2)), ylim=c(0, 1))
lines(sort(dataset2),  empirical(sort(dataset2)), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, Bounds95[,1], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds95[,2], col = "lightgreen",type = "s", lty = 2)
lines(sequence, Bounds999[,1], col = "green",type = "s", lty = 2)
lines(sequence, Bounds999[,2], col = "green",type = "s", lty = 2)
lines(sequence, LogTrans95final, col = "blue")
lines(sequence, LogTrans99final, col = "red")
lines(sequence, LogTrans99first, col = "purple")
lines(sequence, LogTrans95first, col = "orange")
legend("bottomright", legend = c("Final fit 95%", "Final fit 99.9%", "Initial guess 95%", "Initial guess 99.9%"),
       col = c("blue", "red", "orange","purple"), lwd = 2, 
       bty = "n", cex = 1.5)
legend("topleft", legend = c("95% confidence", "99.9% confidence"),
       col = c("lightgreen", "green"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
