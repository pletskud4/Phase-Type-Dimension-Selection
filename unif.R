#matrix til opbevaring af data
library(matrixdist)
library(actuar)
rows <- 10000
matrix_data <- matrix(NA, nrow = rows, ncol = 2)

set.seed(999)  
matrix_data[,1] <- runif(rows, min = 0, max = 1)
matrix_data[,2] <- runif(rows, min = 0, max = 1)

folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste0(folder_path, "unif_10000obs.csv")
write.csv(matrix_data, file = file_name)

