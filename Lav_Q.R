library(readr)
library(matrixdist)
Q <- function(n, data, main_hist) {
  N=length(data)
  limit = ceiling(n*max(data)) 
  ordered_data <- data[order(data, decreasing = F)]
  Q_values <- numeric(length = limit)
  
  for (x in 1:limit) {
    pois_dens <- dpois(x = x-1, lambda = n * ordered_data)
    pois_dens_times_frac <- pois_dens * log((N+1-(1:length(pois_dens) - 1))/(N-(1:length(pois_dens) - 1)))
    Q_values[x] <- 1 - sum(pois_dens_times_frac[-N])
  }
  T_matrix <- matrix(0, nrow = limit, ncol = limit)
  diag(T_matrix) <- -n
  for (i in 1:(limit - 1)) {
    for (j in 1:(ncol(T_matrix) - 1)) {
      if (T_matrix[i, j] == -n) {
        T_matrix[i, j + 1] = n*Q_values[i]
      }
    }
  }
  alpha <- c(1, rep(0, length = limit-1))
  ph2<-ph(structure="general",dimension=limit)
  ph2@pars$alpha = alpha
  ph2@pars$S = T_matrix
  
  calculate_expression <- function(y, n, p, Q_tilde) {
    result <- numeric(length(y))
    
    for (i in seq_along(y)) {
      term <- 0
      for (l in 1:p) {
        inner_product <- prod(Q_tilde[1:(l - 1)])
        poisson_term <- n*dpois(l - 1, lambda =  n * y[i])
        term <- term + (1 - Q_tilde[l]) * inner_product * poisson_term
      }
      
      result[i] <- term
    }
    
    return(result)
  }
  work <- calculate_expression(seq(0, 1.25*max(data), length.out = 100), n, limit, Q_values)
  hist(data, breaks = 100
       , probability = T, col = "darkgreen", ylim = c(0, 1.25*max(1, density(data)$y)), xlim = c(0, 1.25*max(data)), main =main_hist)
  lines(seq(0, 1.25*max(data), length.out = 100), work, col = "blue", type = "p")
  lines(seq(0, 1.25*max(data), length.out = 100), work, col = "blue", type = "l")
  
  quantiles <- seq(0, 1.25*max(data), length.out = 100)  
  cdf_values <- cumsum(work)  
  cdf_values <- cdf_values / sum(work)
  lines(quantiles, cdf_values, type = "step", col= "red")
  lines(ecdf(data))
  
  return(invisible(list(ph2 =ph2, density = work, dim = limit)))
}
########
folder_path <- "/Users/Thesis/ph5_10obs (1)_kopi"

files <- list.files(folder_path, full.names = TRUE)

file_names <- tools::file_path_sans_ext(basename(files))

#Q_150 <- function(data){
#  Q(500, data)
#}
#######
### get datasets
all_datasets <- list()

for (file in files) {
  data <- read.csv(file)

  dataset1 <- data[, 2]  
  dataset2 <- data[, 3]  
  

  datasets <- list(dataset1 = dataset1, dataset2 = dataset2)
  

  all_datasets[[basename(file)]] <- datasets
}
#setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Bounds/Plots")
png(filename = " 1000_Q_1000.png", width = 10, height = 8, units = "in", res = 600)
Q(1000, all_datasets$` multimodal_1000obs.csv`$dataset1, "1000 observations, n=1000")
dev.off()
k <- Q(10, all_datasets$` multi_10000.csv`$dataset1, "10000 observations, n=1000")
k$dim


for (file_name in names(all_datasets)) {
  datasets <- all_datasets[[file_name]]
  

  for (i in 1:length(datasets)) {
    dataset <- datasets[[i]]
    Q_150(dataset)

  }
}

for (file_name in names(all_datasets)) {
  datasets <- all_datasets[[file_name]]
  result <- Q(100, datasets$dataset1, file_name) 
}








