library(matrixdist)

#starter med at lave reproducerbart
set.seed(999)
confidence <- 0.01
change <- 0.01
dimensions <- c(2, 5, 10, 20)

#nu skal vi simulere en kolmogorov fordeling
largest_values <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")[,2]
drift <- 0
variance <- 1
largest_values <- numeric(1000000)
for (i in 1:1000000) {
  increments <- rnorm(100000, mean = drift, sd = sqrt(variance/100000))
  baby_bm <- cumsum(increments)
  bridge_bm <- baby_bm - (seq_along(baby_bm)/100000) * tail(baby_bm, n = 1)
  largest_values[i] <- max(abs(bridge_bm))
}
write.csv(largest_values, file = "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/kolmogorov.csv")

#laver nu matrix med plads til alt 4 fits pr datasæt, så 64*4 rækker
#8 søjler, en til beskrivelse af fit, en til beskrivelse af data
#tre til beskrivelse af hvornå konvergens er sket
#tre til log-likelihood værdi når det sker
results <- matrix(nrow = 256, ncol = 12)
col_names <- c("Fit", "Data", "Iteration DKW", "Iteration KS", "Iteration Convergence", "loglik DKW", "Loglik ks", "loglik convergence", "initial loglik","Average point dist DKW final","Average point dist DKW final","Average point value")
colnames(results) <- col_names
for (i in seq(4, nrow(results), by = 4)) {
  results[i-3, 1] <- "1d" 
  results[i-2, 1] <- "5d"
  results[i-1, 1] <- "10d" 
  results[i, 1] <- "20d" 
}


#nu er vi klar til at fylde den ud
directory <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets"
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
    ph_model1 <- ph(structure = "general", dimension = dimensions[j])
    fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
    results[(i-1)*4+j,9]=logLik(fit1)
    logLik1 <- logLik(fit1)
    for (iter in seq(1, 50000, by = 1)) {
      fit1 <- fit(fit1, y = data1, stepsEM = 1, every = 1)
      fitted1 <- cdf(fit1, sort(data1))
      #kode checker om indenfor dkw bounds
      if (!old_dkw1) {
        within_dkw1 <- all(fitted1 >= dkw_lower_bound1 & fitted1 <= dkw_upper_bound1)
        if (within_dkw1) {
          results[(i-1)*4+j,3]=iter
          results[(i-1)*4+j,6]=logLik(fit1)
          fitteddkw1=fitted1
          old_dkw1 <- within_dkw1
        }
      }
      #kode checker om indenfor ks bounds
      if (!old_ks1) {
        within_ks1 <- all(fitted1 >= ks_lower_bound1 & fitted1 <= ks_upper_bound1)
        if (within_ks1) {
          results[(i-1)*4+j,4]=iter
          results[(i-1)*4+j,7]=logLik(fit1)
          fittedks1=fitted1
          old_ks1 <- within_ks1
        }
      }
      #kode nedenfor checker konvergens
      if (abs(logLik1 - logLik(fit1)) < change) {
        results[(i-1)*4+j,5]=iter
        results[(i-1)*4+j,8]=logLik(fit1)
        results[(i-1)*4+j,10]=sum(abs(fitteddkw1-fitted1))/length(data1)
        results[(i-1)*4+j,11]=sum(abs(fittedks1-fitted1))/length(data1)
        results[(i-1)*4+j,12]=mean(fitted1)
        break  
      }
      logLik1 <- logLik(fit1)
    }
    #data2
    ph_model2 <- ph(structure = "general", dimension = dimensions[j])
    fit2 <- fit(ph_model2, y = data2, stepsEM = 1)
    results[(i-1)*4+j+128,9]=logLik(fit2)
    logLik2 <- logLik(fit2)
    for (iter in seq(1, 50000, by = 1)) {
      fit2 <- fit(fit2, y = data2, stepsEM = 1, every = 1)
      fitted2 <- cdf(fit2, sort(data2))
      #kode checker om indenfor dkw bounds
      if (!old_dkw2) {
        within_dkw2 <- all(fitted2 >= dkw_lower_bound2 & fitted2 <= dkw_upper_bound2)
        if (within_dkw2) {
          results[(i-1)*4+j+128,3]=iter
          results[(i-1)*4+j+128,6]=logLik(fit2)
          fitteddkw2=fitted2
          old_dkw2 <- within_dkw2
        }
      }
      #kode checker om indenfor ks bounds
      if (!old_ks2) {
        within_ks2 <- all(fitted2 >= ks_lower_bound2 & fitted2 <= ks_upper_bound2)
        if (within_ks2) {
          results[(i-1)*4+j+128,4]=iter
          results[(i-1)*4+j+128,7]=logLik(fit2)
          fittedks2=fitted2
          old_ks2 <- within_ks2
        }
      }
      #kode nedenfor checker konvergens
      if (abs(logLik2 - logLik(fit2)) < change) {
        results[(i-1)*4+j+128,5]=iter
        results[(i-1)*4+j+128,8]=logLik(fit2)
        results[(i-1)*4+j+128,10]=sum(abs(fitteddkw2-fitted2))/length(data1)
        results[(i-1)*4+j+128,11]=sum(abs(fittedks2-fitted2))/length(data1)
        results[(i-1)*4+j+128,12]=mean(fitted2)
        break  
      }
      logLik2 <- logLik(fit2)
    }
  }
}
#gemmer resultatet
folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/"
file_name <- paste0(folder_path, "results.csv")
write.csv(results, file = file_name)
results <- read.csv("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/results.csv")


#vil gerne regne bedste konfidensniveau
niveauer <- c(0.01, 0.05, 0.1, 0.25)
obser <- c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)
forskelle <- matrix(NA, nrow = length(niveauer), ncol = length(obser))
rownames(forskelle) <- as.character(niveauer)
colnames(forskelle) <- as.character(obser)
for (i in 1:length(niveauer)) {
  for (j in 1:length(obser)) {
    forskelle[i, j] <- sqrt((1 / (2 * obser[j])) * log(2 / (niveauer[i]))) - sqrt(1 /  obser[j])*quantile(largest_values, probs = 1-niveauer[i])
  }
}
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Confidence levels.png", width = 10, height = 8, units = "in", res = 600)
barplot(as.matrix(forskelle), beside = TRUE, 
        col = c("red", "green", "blue", "orange"),
        xlab = "Observations", ylab = "DKW-KS",
        names.arg = as.character(obser),
        col.axis = "darkgreen")
legend("topright", legend = rownames(forskelle), fill = c("red", "green", "blue", "orange"), title = "Confidence levels",cex = 2)
dev.off()



#undersøger forskel i størrelse nu
# Set up the data
obser <- c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)
confidence <- 0.01
# Calculate the differences
forskelle <- sapply(obser, function(obs) {
  sqrt((1 / (2 * obs)) * log(2 / confidence)) / (sqrt(1 / obs) * quantile(largest_values, probs = 1 - confidence))
})
# Plotting
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Relative Confidence levels.png", width = 10, height = 8, units = "in", res = 600)
plot(obser, forskelle, type = "h",
     col = "darkgreen", lwd = 10,
     xlab = "Observations", ylab = "DKW/KS", ylim=c(1,1.002))
dev.off()




#plot af fitting af 1000 observationer
data_frame <- read.csv(csv_files[3])
data1 <- data_frame[, 2]
empirical1 <- ecdf(data1)
quantiles1 <- empirical1(sort(data1))
epsilon <- sqrt((1 / (2 * length(data1))) * log(2 / (confidence)))
dkw_lower_bound1 <- pmax(quantiles1 - epsilon, 0)
dkw_upper_bound1 <- pmin(quantiles1 + epsilon, 1)
size <- sqrt(1 /  length(data1))*quantile(largest_values, probs = 1-confidence)
ks_lower_bound1 <- pmax(quantiles1 - size, 0)
ks_upper_bound1 <- pmin(quantiles1 + size, 1)
points <- seq(0, max(data1), length.out = 10000)
for (j in 1:4) {
  #så er alt det nemme fyldt, og nu skal vi her have resultater
  #data1
  ph_model1 <- ph(structure = "general", dimension = dimensions[j])
  fit1 <- fit(ph_model1, y = data1, stepsEM = 1)
  logLik1 <- logLik(fit1)
  for (iter in seq(1, 50000, by = 1)) {
    fit1 <- fit(fit1, y = data1, stepsEM = 1, every = 1)
    #kode nedenfor checker konvergens
    if (abs(logLik1 - logLik(fit1)) < change) {
      results <- cdf(fit1, points)
      vector_name <- paste0("Dimensions", dimensions[j])
      assign(vector_name, results)
      break  
    }
    logLik1 <- logLik(fit1)
  }
}
#klar til at plotte nu
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Multimodal fit.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0, max(data1)), ylim=c(0, 1))
lines(sort(data1),  empirical1(sort(data1)), col = "black", type = "s",lwd = 2)
lines(sort(data1), dkw_lower_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), dkw_upper_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), ks_lower_bound1, col = "darkgreen",type = "s", lty = 2)
lines(sort(data1), ks_upper_bound1, col = "darkgreen",type = "s", lty = 2)
lines(points, Dimensions2, col = "blue")
lines(points, Dimensions5, col = "red")
lines(points, Dimensions10, col = "grey")
lines(points, Dimensions20, col = "purple")
legend("topleft", legend = c("DKW", "KS"), col = c("lightgreen", "darkgreen"), lty = 2, 
       bty = "n", cex = 1.5)  
legend("bottomright", legend = c("Dimensions 2", "Dimensions 5", "Dimensions 10", "Dimensions 20"),
       col = c("blue", "red", "grey", "purple"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()

#laver også et zoomet plot
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Multimodal fit zoom.png", width = 10, height = 8, units = "in", res = 600)
plot(0, type = "n",
     xlab = "Data", ylab = "Cumulative Probability", col = "black", 
     xlim = c(0.6, 0.61), ylim=c(0.7025, 0.703))
lines(sort(data1),  empirical1(sort(data1)), col = "black", type = "s",lwd = 2)
lines(sort(data1), dkw_lower_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), dkw_upper_bound1, col = "lightgreen",type = "s", lty = 2)
lines(sort(data1), ks_lower_bound1, col = "darkgreen",type = "s", lty = 2)
lines(sort(data1), ks_upper_bound1, col = "darkgreen",type = "s", lty = 2)
lines(points, Dimensions2, col = "blue")
lines(points, Dimensions5, col = "red")
lines(points, Dimensions10, col = "grey")
lines(points, Dimensions20, col = "purple")
legend("topleft", legend = c("DKW", "KS"), col = c("lightgreen", "darkgreen"), lty = 2, 
       bty = "n", cex = 1.5)  
legend("bottomright", legend = c("Dimensions 2", "Dimensions 5", "Dimensions 10", "Dimensions 20"),
       col = c("blue", "red", "grey", "purple"), lwd = 2, 
       bty = "n", cex = 1.5)
dev.off()


#Fremskaffer tal som er anvendt i tabellen
count <- sum(!is.na(results[, 3]))
countdkw <- sum(!is.na(results[, 4]))
countks <- sum(!is.na(results[, 5]))
mean <- mean(results[!is.na(results[, 6]), 6])
dkw_mean <- mean(results[!is.na(results[, 4]), 4])
ks_mean <- mean(results[!is.na(results[, 5]), 5])
sum <- sum(results[!is.na(results[, 9]), 9])/sum(results[!is.na(results[, 10]), 10])
dkw_sum <- sum(results[!is.na(results[, 7]), 7])/sum(results[!is.na(results[, 7]), 10])
ks_sum <- sum(results[!is.na(results[, 8]), 8])/sum(results[!is.na(results[, 8]), 10])
sum_end <- sum(results[!is.na(results[, 9]), 9])/sum(results[!is.na(results[, 10]), 10])
dkw_sum_end <- sum(results[!is.na(results[, 7]), 9])/sum(results[!is.na(results[, 7]), 7])
ks_sum_end <- sum(results[!is.na(results[, 8]), 9])/sum(results[!is.na(results[, 8]), 8])

relativedkw <- numeric(196)
relativeks <- numeric(196)
iterdkw <- numeric(196)
iterks <- numeric(196)
for (i in 1:196){
  if (!is.na(results[i, 7])){
    relativedkw[i] <- (results[i,7]-results[i, 10])/(results[i,9]-results[i, 10])
    iterdkw[i] <- (results[i,4])
  }
  if (!is.na(results[i, 8])){
    relativeks[i] <- (results[i,8]-results[i, 10])/(results[i,9]-results[i, 10])
    iterks[i] <- (results[i,5])
  }
}
mean(relativedkw[!is.na(relativedkw) & relativedkw != 0])
mean(relativeks[!is.na(relativeks) & relativeks != 0])
sum(!is.na(relativedkw) & relativedkw != 0)
sum(!is.na(relativeks) & relativeks != 0)
mean(iterdkw[!is.na(iterdkw) & iterdkw != 0])
mean(iterks[!is.na(iterks) & iterks != 0])






#fremskaffer iteration tabel
countmatrixdkw <- matrix(NA, nrow = 3, ncol = 4)
countmatrixks <- matrix(NA, nrow = 3, ncol = 4)
for (j in 13:16) {
  non_blank_values <- 0
  count <- 0
  for (i in seq(j, 256, by = 16)) {
    if (!is.na(results[i, 4])) {  
      non_blank_values = non_blank_values + (results[i,8]-results[i,10])/(results[i,9]-results[i,10])
      count = count +1
    }
  }
  countmatrixks[1,j-12] = count
}


#kolmogorov statistic
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = "Kolmogorov statistic.png", width = 10, height = 8, units = "in", res = 600)
hist(largest_values, col = "darkgreen", breaks = 200, freq = F, xlab = "Kolmogorov-statistic")
dev.off()

