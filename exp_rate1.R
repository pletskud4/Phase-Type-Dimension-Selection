
library(matrixdist)
rows <- 10
dim <- 5
matrix_data <- matrix(NA, nrow = rows, ncol = 2)


set.seed(999)  
matrix_data[,1] <- rexp(10, rate=1)  
matrix_data[,2] <- rexp(10, rate=1)  

folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste0(folder_path, "exp_rate1_10obs.csv")
write.csv(matrix_data, file = file_name)
