
rows <- 10
matrix_data <- matrix(NA, nrow = rows, ncol = 2)


set.seed(999)  
matrix_data[,1] <- rlnorm(rows, meanlog = 0, sdlog = 0.1)
matrix_data[,2] <- rlnorm(rows, meanlog = 0, sdlog = 0.1)


folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste0(folder_path, "ln_10obs.csv")
write.csv(matrix_data, file = file_name)
