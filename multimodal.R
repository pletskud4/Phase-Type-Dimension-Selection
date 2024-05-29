rows <- 10000
matrix_data <- matrix(NA, nrow = rows, ncol = 2)


set.seed(999)  
data <- c(rgamma(20000,shape=30,rate=40),rgamma(10000,shape=10,rate=40),rgamma(5000,shape=2,rate=40))
matrix_data[,1] <- sample(data, size = rows, replace = FALSE)
matrix_data[,2] <- sample(data, size = rows, replace = FALSE)

folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste(folder_path, "multimodal_10000obs.csv")
write.csv(matrix_data, file = file_name)

