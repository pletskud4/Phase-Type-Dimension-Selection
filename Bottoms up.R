library(matrixdist)


set.seed(999)
confidence <- 0.25
change <- 0.01
dimensions <- c(2, 5, 10, 20)
DADfunction <- function(dist, n) { 
  GM <- solve((-dist@pars$S)) 
  CV <- dist@pars$alpha %*% GM 
  oldT <- dist@pars$S
  oldalpha <- dist@pars$alpha
  copystatenumber <- which.min(CV) 
  oldT[copystatenumber,] <- dist@pars$S[length(dist@pars$alpha),]
  oldT[length(dist@pars$alpha),] <- dist@pars$S[copystatenumber,]
  copystate <- oldT[,copystatenumber]
  oldT[,copystatenumber] <- oldT[,length(dist@pars$alpha)]
  oldT[,length(dist@pars$alpha)] <- copystate
  
  oldalpha[copystatenumber] <- dist@pars$alpha[length(dist@pars$alpha)]
  oldalpha[length(dist@pars$alpha)] <- dist@pars$alpha[copystatenumber]

  Tmatrix <- matrix(NA, nrow = (length(dist@pars$alpha)+n), ncol = length(dist@pars$alpha)+n)
  Alphavector <- vector("numeric", length = length(dist@pars$alpha)+n)
  epsilon <- min(oldT[length(dist@pars$alpha), oldT[length(dist@pars$alpha),] > 0],abs(oldT[length(dist@pars$alpha),length(dist@pars$alpha)]))
  for (i in seq_len(length(Alphavector))) {
    for (j in seq_len(length(Alphavector))) {
      if (i<length(dist@pars$alpha)){
        if (j<length(dist@pars$alpha)){
          if (i==j){
            Tmatrix[i,j]<- oldT[i,j]
          } else {
            Tmatrix[i,j]<- oldT[i,j]
          }
        } else {
          Tmatrix[i,j] <- oldT[i, length(dist@pars$alpha)]/(n+1)
        }
      } else {
        if (j<length(dist@pars$alpha)){
          Tmatrix[i,j] <- oldT[length(dist@pars$alpha),j]
        } else {
          if (i==j) {
            Tmatrix[i,j] <- (oldT[length(dist@pars$alpha),length(dist@pars$alpha)]-(n+1)*epsilon)
          } else {
            Tmatrix[i,j] <- epsilon
          }
        }
      }
    }
    if (i<length(dist@pars$alpha)){
      Alphavector[i] <- oldalpha[i]
      } else {
      Alphavector[i] <- oldalpha[length(dist@pars$alpha)]/(n+1)
    }
  }
  return(ph(alpha = Alphavector, S = Tmatrix))
}

#nu skal vi simulere en kolmogorov fordeling
largest_values <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")[,2]

#laver nu matrix med plads til alt 4 fits pr datasæt, så 64*4 rækker
#8 søjler, en til beskrivelse af fit, en til beskrivelse af data
#tre til beskrivelse af hvornå konvergens er sket
#tre til log-likelihood værdi når det sker
results <- matrix(nrow = 256, ncol = 20)
col_names <- c("Fit", "Data", "Iteration DKW", "Iteration KS", "Iteration Convergence", "loglik DKW", "Loglik ks", "loglik convergence", "initial loglik","Average point dist DKW final","Average point dist DKW final","Average point value","Time -1 state", "Iterations -1 state", "loglik initial -1 state","loglik conv -1 state","in ks initial -1 state","in ks conv -1 state", "in ks initial","in ks conv")
colnames(results) <- col_names
for (i in seq(4, nrow(results), by = 4)) {
  results[i-3, 1] <- "2d" 
  results[i-2, 1] <- "5d"
  results[i-1, 1] <- "10d" 
  results[i, 1] <- "20d" 
}



#nu er vi klar til at fylde den ud
directory<- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets"
csv_files <- list.files(directory, pattern = "*.csv", full.names = TRUE)
for (i in seq_along(csv_files)) {
  #i yderste loop behøves blot det som kan ses fra data
  if (((i-1) %% 4) == 0) {
    next
  }
  data_frame <- read.csv(csv_files[i])
  data1 <- data_frame[, 2]
  empirical1 <- ecdf(data1)
  quantiles1 <- empirical1(sort(data1))
  data2 <- data_frame[, 3]
  empirical2 <- ecdf(data2)
  quantiles2 <- empirical2(sort(data2))
  epsilon <- sqrt((1 / (2 * length(data1))) * log(2 / (confidence)))
  dkw_lower_bound1 <- pmax(quantiles1 - epsilon, 0)
  dkw_upper_bound1 <- pmin(quantiles1 + epsilon, 1)
  dkw_lower_bound2 <- pmax(quantiles2 - epsilon, 0)
  dkw_upper_bound2 <- pmin(quantiles2 + epsilon, 1)
  size <- sqrt(1 /  length(data1))*quantile(largest_values, probs = 1-confidence)
  ks_lower_bound1 <- pmax(quantiles1 - size, 0)
  ks_upper_bound1 <- pmin(quantiles1 + size, 1)
  ks_lower_bound2 <- pmax(quantiles2 - size, 0)
  ks_upper_bound2 <- pmin(quantiles2 + size, 1)
  for (j in 1:4) {
    if (i == 1) {
      if (j == 1){
        next
      }
    }
    old_dkw1 <- FALSE
    old_ks1 <- FALSE
    old_dkw2 <- FALSE
    old_ks2 <- FALSE
    fitteddkw1 <- numeric(length(data1))
    fitteddkw2 <- numeric(length(data1))
    fittedks1 <- numeric(length(data1))
    fittedks2 <- numeric(length(data1))
    print(i)
    print(j)
    results[(i-1)*4+j,2]=paste0(basename(csv_files[i]), "&1")
    results[(i-1)*4+j+128,2]=paste0(basename(csv_files[i]), "&2")
    #så er alt det nemme fyldt, og nu skal vi her have resultater
    #data1
    ph_model1 <- ph(structure = "general", dimension = dimensions[j]-1)
    fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
    
    sequence1 <- seq(min(data1), max(data1), length.out = 10000)
    results[(i-1)*4+j,17] <- sum(cdf(fit1, sequence1) <= pmin(empirical1(sequence1) + size, 1) & cdf(fit1, sequence1) >= pmax(empirical1(sequence1) - size, 0))/10000
    results[(i-1)*4+j,15] <- logLik(fit1)
    start_time1 <- Sys.time()
    logLik1 <- logLik(fit1)
    for (iter in seq(1, 50000, by = 1)) {
      fit1 <- fit(fit1, y = data1, stepsEM = 1, every = 1)
      if (abs(logLik1 - logLik(fit1)) < change) {
        results[(i-1)*4+j,13] <- Sys.time()-start_time1
        results[(i-1)*4+j,14] <- iter
        results[(i-1)*4+j,18] <- sum(cdf(fit1, sequence1) <= pmin(empirical1(sequence1) + size, 1) & cdf(fit1, sequence1) >= pmax(empirical1(sequence1) - size, 0))/10000
        results[(i-1)*4+j,16] <- logLik(fit1)
        ph_model1fin <- DADfunction(fit1,1)
        print("HALL3O")
        fit1fin <- fit(ph_model1fin, y = data1, stepsEM = 1)
        print("HALL2O")
        results[(i-1)*4+j,19] <- sum(cdf(fit1fin, sequence1) <= pmin(empirical1(sequence1) + size, 1) & cdf(fit1fin, sequence1) >= pmax(empirical1(sequence1) - size, 0))/10000
        results[(i-1)*4+j,9]=logLik(fit1fin)
        logLik1fin <- logLik(fit1fin)
        for (iter2 in seq(1, 50000, by = 1)) {
          fit1fin <- fit(fit1fin, y = data1, stepsEM = 1, every = 1)
          fitted1 <- cdf(fit1fin, sort(data1))
          
          if (!old_dkw1) {
            within_dkw1 <- all(fitted1 >= dkw_lower_bound1 & fitted1 <= dkw_upper_bound1)
            if (within_dkw1) {
              results[(i-1)*4+j,3]=iter2
              results[(i-1)*4+j,6]=logLik(fit1fin)
              fitteddkw1=fitted1
              old_dkw1 <- within_dkw1
            }
          }
          #kode checker om indenfor ks bounds
          if (!old_ks1) {
            within_ks1 <- all(fitted1 >= ks_lower_bound1 & fitted1 <= ks_upper_bound1)
            if (within_ks1) {
              results[(i-1)*4+j,4]=iter2
              results[(i-1)*4+j,7]=logLik(fit1fin)
              fittedks1=fitted1
              old_ks1 <- within_ks1
            }
          }
          #kode nedenfor checker konvergens
          if (abs(logLik1fin - logLik(fit1fin)) < change) {
            results[(i-1)*4+j,5]=iter2
            results[(i-1)*4+j,8]=logLik(fit1fin)
            results[(i-1)*4+j,10]=sum(abs(fitteddkw1-fitted1))/length(data1)
            results[(i-1)*4+j,11]=sum(abs(fittedks1-fitted1))/length(data1)
            results[(i-1)*4+j,12]=mean(fitted1)
            results[(i-1)*4+j,20] <- sum(cdf(fit1fin, sequence1) <= pmin(empirical1(sequence1) + size, 1) & cdf(fit1fin, sequence1) >= pmax(empirical1(sequence1) - size, 0))/10000
            break  
          }
          logLik1fin <- logLik(fit1fin)
        }
        break  
      }
      logLik1 <- logLik(fit1)
    }
    #data2
    ph_model2 <- ph(structure = "general", dimension = dimensions[j]-1)
    fit2 <- fit(ph_model2, y = data1, stepsEM = 1)
    sequence2 <- seq(min(data2), max(data2), length.out = 10000)
    results[(i-1)*4+j+128,17] <- sum(cdf(fit2, sequence2) <= pmin(empirical2(sequence2) + size, 1) & cdf(fit2, sequence2) >= pmax(empirical2(sequence2) - size, 0))/10000
    results[(i-1)*4+j+128,15] <- logLik(fit2)
    start_time2 <- Sys.time()
    logLik2 <- logLik(fit2)
    for (iter in seq(1, 50000, by = 1)) {
      fit2 <- fit(fit2, y = data2, stepsEM = 1, every = 1)
      if (abs(logLik2 - logLik(fit2)) < change) {
        results[(i-1)*4+j+128,13] <- Sys.time()-start_time2
        results[(i-1)*4+j+128,14] <- iter
        results[(i-1)*4+j+128,18] <- sum(cdf(fit2, sequence2) <= pmin(empirical2(sequence2) + size, 1) & cdf(fit2, sequence2) >= pmax(empirical2(sequence2) - size, 0))/10000
        results[(i-1)*4+j+128,16] <- logLik(fit2)
        ph_model2fin <- DADfunction(fit2,1)
        fit2fin <- fit(ph_model2fin, y = data1, stepsEM = 1)
        results[(i-1)*4+j+128,19] <- sum(cdf(fit2fin, sequence2) <= pmin(empirical2(sequence2) + size, 1) & cdf(fit2fin, sequence2) >= pmax(empirical2(sequence2) - size, 0))/10000
        results[(i-1)*4+j+128,9]=logLik(fit2fin)
        logLik2fin <- logLik(fit2fin)
        for (iter2 in seq(1, 50000, by = 1)) {
          fit2fin <- fit(fit2fin, y = data2, stepsEM = 1, every = 1)
          fitted2 <- cdf(fit2fin, sort(data1))
          
          if (!old_dkw2) {
            within_dkw2 <- all(fitted2 >= dkw_lower_bound2 & fitted2 <= dkw_upper_bound2)
            if (within_dkw2) {
              results[(i-1)*4+j+128,3]=iter2
              results[(i-1)*4+j+128,6]=logLik(fit2fin)
              fitteddkw2=fitted2
              old_dkw2 <- within_dkw2
            }
          }
          #kode checker om indenfor ks bounds
          if (!old_ks2) {
            within_ks2 <- all(fitted2 >= ks_lower_bound2 & fitted2 <= ks_upper_bound2)
            if (within_ks2) {
              results[(i-1)*4+j+128,4]=iter2
              results[(i-1)*4+j+128,7]=logLik(fit2fin)
              fittedks2=fitted2
              old_ks2 <- within_ks2
            }
          }
          #kode nedenfor checker konvergens
          if (abs(logLik2fin - logLik(fit2fin)) < change) {
            results[(i-1)*4+j+128,5]=iter2
            results[(i-1)*4+j+128,8]=logLik(fit2fin)
            results[(i-1)*4+j+128,10]=sum(abs(fitteddkw2-fitted2))/length(data1)
            results[(i-1)*4+j+128,11]=sum(abs(fittedks2-fitted2))/length(data1)
            results[(i-1)*4+j+128,12]=mean(fitted2)
            results[(i-1)*4+j+128,20] <- sum(cdf(fit2fin, sequence2) <= pmin(empirical2(sequence2) + size, 1) & cdf(fit2fin, sequence2) >= pmax(empirical2(sequence2) - size, 0))/10000
            break  
          }
          logLik2fin <- logLik(fit2fin)
        }
        break  
      }
      logLik2 <- logLik(fit2)
    }
  }
}
#gemmer resultatet
folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/"
file_name <- paste0(folder_path, "resultsbottomup.csv")
write.csv(results, file = file_name)
results <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/resultsbottomup.csv")

















#første plot i sektionen
data_frame <- read.csv(csv_files[1])
data1 <- data_frame[, 2]
empirical1 <- ecdf(data1)
quantiles1 <- empirical1(sort(data1))
size <- sqrt(1 /  length(data1))*quantile(largest_values, probs = 1-confidence)
ks_lower_bound1 <- pmax(quantiles1 - size, 0)
ks_upper_bound1 <- pmin(quantiles1 + size, 1)
sequence <- seq(min(data1), max(data1), length.out = 10000)
#finder færdigt fit i 19 dim
ph_model1 <- ph(structure = "general", dimension = 19)
fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
logLik1 <- logLik(fit1)
for (iter in seq(1, 50000, by = 1)) {
  fit1 <- fit(fit1, y = data1, stepsEM = 1, every = 1)
  if (abs(logLik1 - logLik(fit1)) < change) {
    break
  }
  logLik1 <- logLik(fit1)
}
converge19 <- cdf(fit1,sequence)
#laver random første fit
ph_model2 <- ph(structure = "general", dimension = 20)
fit2 <- fit(ph_model2, y = data1, stepsEM = 1)
random20 <- cdf(fit2,sequence)
#laver planlagt første fit
ph_model3 <- DADfunction(fit1,1)
fit3 <- fit(ph_model3, y = data1, stepsEM = 1)
first20 <- cdf(fit3,sequence)
#laver færdigt fit i 20 dim
fit4 <- fit(ph_model3, y = data1, stepsEM = 1)
logLik4 <- logLik(fit4)
for (iter in seq(1, 50000, by = 1)) {
  fit4 <- fit(fit4, y = data1, stepsEM = 1, every = 1)
  if (abs(logLik4 - logLik(fit4)) < change) {
    break
  }
  logLik4 <- logLik(fit4)
}
converge20 <- cdf(fit4,sequence)
#og nu kan vi endelig lave plottet
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Bottom up example.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0, max(data1)), ylim=c(0, 1))
lines(sort(data1),  empirical1(sort(data1)), col = "darkgreen", type = "s",lwd = 2)
lines(sort(data1), ks_lower_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), ks_upper_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sequence, converge19, col = "blue")
lines(sequence, converge20, col = "red")
lines(sequence, first20, col = "grey")
lines(sequence, random20, col = "purple")
legend("bottomright", legend = c("Fit at convergence 19d", "Fit at convergence 20d", "First guess 20d", "Random 20d"),
       col = c("blue", "red", "grey", "purple"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()
#Hernæst laves plot hvor vi zoomer ind
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Bottom up example zoom.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0.79, 0.8), ylim=c(0.78, 0.8))
lines(sort(data1),  empirical1(sort(data1)), col = "darkgreen", type = "s",lwd = 2)
lines(sort(data1), ks_lower_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), ks_upper_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sequence, converge19, col = "blue")
lines(sequence, converge20, col = "red")
lines(sequence, first20, col = "grey")
lines(sequence, random20, col = "purple")
legend("topleft", legend = c("Fit at convergence 19d", "Fit at convergence 20d", "First guess 20d", "Random 20d"),
       col = c("blue", "red", "grey", "purple"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()








#her samler vi de tal som bruges til vores resultater
non_na_count <- 0
pointsinitial <- 0
pointsnineteen <- 0
for (i in 1:nrow(results)) {
  if (!anyNA(results[i, 3])) {
    non_na_count <- non_na_count + 1
    pointsinitial <- pointsinitial + results[i, 20] / results[i, 21]
  }
}
pointsinitial/non_na_count
#skal have det samme for et tilfældigt første gæt
randoms <- matrix(nrow = 256, ncol = 3)
for (i in seq(4, nrow(randoms), by = 4)) {
  randoms[i-3, 1] <- "2d" 
  randoms[i-2, 1] <- "5d"
  randoms[i-1, 1] <- "10d" 
  randoms[i, 1] <- "20d" 
}
for (i in seq_along(csv_files)) {
  #i yderste loop behøves blot det som kan ses fra data
  if (((i-1) %% 4) == 0) {
    next
  }
  data_frame <- read.csv(csv_files[i])
  data1 <- data_frame[, 2]
  empirical1 <- ecdf(data1)
  quantiles1 <- empirical1(sort(data1))
  data2 <- data_frame[, 3]
  empirical2 <- ecdf(data2)
  quantiles2 <- empirical2(sort(data2))
  size <- sqrt(1 /  length(data1))*quantile(largest_values, probs = 1-confidence)
  ks_lower_bound1 <- pmax(quantiles1 - size, 0)
  ks_upper_bound1 <- pmin(quantiles1 + size, 1)
  ks_lower_bound2 <- pmax(quantiles2 - size, 0)
  ks_upper_bound2 <- pmin(quantiles2 + size, 1)
  fitteddkw1 <- numeric(length(data1))
  fitteddkw2 <- numeric(length(data1))
  fittedks1 <- numeric(length(data1))
  fittedks2 <- numeric(length(data1))
  for (j in 1:4) {
    if (i == 1) {
      if (j == 1){
        next
      }
    }
    randoms[(i-1)*4+j,2]=paste0(basename(csv_files[i]), "&1")
    randoms[(i-1)*4+j+128,2]=paste0(basename(csv_files[i]), "&2")
    #så er alt det nemme fyldt, og nu skal vi her have resultater
    #data1
    ph_model1 <- ph(structure = "general", dimension = dimensions[j])
    fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
    sequence1 <- seq(min(data1), max(data1), length.out = 10000)
    randoms[(i-1)*4+j,3] <- as.numeric(sum(cdf(fit1, sequence1) <= pmin(empirical1(sequence1) + size, 1) & cdf(fit1, sequence1) >= pmax(empirical1(sequence1) - size, 0))/10000)
    #data2
    ph_model2 <- ph(structure = "general", dimension = dimensions[j]-1)
    fit2 <- fit(ph_model2, y = data1, stepsEM = 1)
    sequence2 <- seq(min(data2), max(data2), length.out = 10000)
    randoms[(i-1)*4+j+128,3] <- as.numeric(sum(cdf(fit2, sequence2) <= pmin(empirical2(sequence2) + size, 1) & cdf(fit2, sequence2) >= pmax(empirical2(sequence2) - size, 0))/10000)
  }
}
non_na_count <- 0
pointsinitial <- 0
for (i in 1:nrow(randoms)) {
  if (!anyNA(randoms[i, 3])) {
    non_na_count <- non_na_count + 1
    pointsinitial <- pointsinitial + as.numeric(randoms[i, 3]) / results[i, 21]
  }
}
pointsinitial/non_na_count





#her gøres det så for hver dimension
non_na_count <- 0
pointsinitial <- 0
for (i in 1:nrow(results)) {
  if (!is.na(results[i, 3]) && results[i, 2] == "20d" ) {
    non_na_count <- non_na_count + 1
    pointsinitial <- pointsinitial + results[i, 20] / results[i, 21]
  }
}
pointsinitial/non_na_count
#og så også for et random gæt
non_na_count <- 0
pointsinitial <- 0
for (i in 1:nrow(randoms)) {
  if (!is.na(randoms[i, 3]) && randoms[i, 1] == "20d" ) {
    non_na_count <- non_na_count + 1
    pointsinitial <- pointsinitial + as.numeric(randoms[i, 3]) / results[i, 21]
  }
}
pointsinitial/non_na_count









#hernæst skal vi have tal på iterationer
results2 <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/results.csv")
non_na_count <- 0
convcurrent <- 0
total <- 0
total2 <- 0
for (i in 1:nrow(results)) {
  if (!is.na(results[i, 3]) ) {
    non_na_count <- non_na_count + 1
    convcurrent <- convcurrent + results[i, 6] / results2[i, 6]
    total <- total + results2[i, 6]
    total2 <- total2 + results[i,6]
  }
}
convcurrent/non_na_count
total/non_na_count
#og gør så det samme for forskellige mængder dimensioner
non_na_count <- 0
convcurrent <- 0
total <- 0
total2 <- 0
for (i in 1:nrow(results)) {
  if (!is.na(results[i, 3]) && results[i, 2] == "2d" ) {
    non_na_count <- non_na_count + 1
    convcurrent <- convcurrent + results[i, 6] / results2[i, 6]
    total <- total + results2[i, 6]
    total2 <- total2 + results[i,6]
  }
}
convcurrent/non_na_count
total/non_na_count
convcurrent
non_na_count
total
total2
sum(results2[,6], na.rm = TRUE)






#sidste plot. Her looper vi indtil vi kommer ind i konfidensintervallet
data_frame <- read.csv(csv_files[2])
data1 <- data_frame[, 2]
empirical1 <- ecdf(data1)
quantiles1 <- empirical1(sort(data1))
size <- sqrt(1 /  length(data1))*quantile(largest_values, probs = 1-confidence)
ks_lower_bound1 <- pmax(quantiles1 - size, 0)
ks_upper_bound1 <- pmin(quantiles1 + size, 1)
sequence <- seq(min(data1), max(data1), length.out = 10000)
ph_model1 <- ph(structure = "general", dimension = 2)
fit1 <- fit(ph_model1, y = data1, stepsEM = 5000)
cdfs <- matrix(NA, nrow = 10000, ncol = 50)
cdfs2 <- matrix(NA, nrow = 10000, ncol = 50)
iterations <- rep(NA, 50)
iterations2 <- rep(NA, 50)
cdfs[,2] <- cdf(fit1,sequence)
cdfs2[,2] <- cdf(fit1,sequence)

for (i in 3:20) {
  ph_model1 <- DADfunction(fit1,1)
  fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
  logLik1 <- logLik(fit1)
  logLik2 <- logLik(fit2)
  counter <- 0
  counter2 <- 0
  størst1 <- 0
  størst2 <- 0
  samlet1 <- 0
  samlet2 <- 0
  for (iter in 1:50000){
    if (counter == 0){
      fit1 <- fit(fit1, y = data1, stepsEM = 1)
        if (abs((logLik(fit1)-logLik1)) < change) {
            iterations[i] <- iter
            counter <- 1
            cdfs[,i] <- cdf(fit1,sequence)
      }
      størst1 <- max(abs(logLik(fit1)-logLik1), størst1)
      samlet1 <- samlet1 + abs(logLik(fit1)-logLik1)
      logLik1 <- logLik(fit1)
    }
    if (counter == 1){
      break
    }
  }
}
folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/"
write.csv(cdfs, file = paste0(folder_path, "bottomupcdfs.csv"))
write.csv(cdfs2, file = paste0(folder_path, "bottomupcdfs2.csv"))
write.csv(iterations, file = paste0(folder_path, "bottomupiterations.csv"))
write.csv(iterations2, file = paste0(folder_path, "bottomupiterations2.csv"))
cdfs <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/bottomupcdfs.csv")

setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
#finally we can create our plots
png(filename = "Bottom-up full.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability by random start", col = "black", 
     xlim = c(0, max(data1)), ylim=c(0, 1))
lines(sort(data1),  empirical1(sort(data1)), col = "darkgreen", type = "s",lwd = 2)
lines(sequence, pmax(empirical1(sequence) - size, 0), col = "lightgreen",type = "s", lty = 2)
lines(sequence, pmin(empirical1(sequence) + size, 1), col = "lightgreen",type = "s", lty = 2)
lines(sequence, cdfs[,2], col = "black")
lines(sequence, cdfs2[,2], col = rainbow(4)[1])
lines(sequence, cdfs2[,5], col = rainbow(4)[2])
lines(sequence, cdfs2[,10], col = rainbow(4)[3])
lines(sequence, cdfs2[,20], col = rainbow(4)[4])
legend_labels <- character(5)
legend_labels[1] <- "DAD-method fits"
for (i in 2:5) {
  legend_labels[i] <- paste0(dimensions[i-1], "d random initial parameters")
}
legend("bottomright", legend = legend_labels, col = rainbow(4), lty = 1, cex = 1)
dev.off()
iterations
sum(cdfs[,20] <= pmin(empirical1(sequence) + size, 1) & cdfs[,20] >= pmax(empirical1(sequence) - size, 0))
cdfs
