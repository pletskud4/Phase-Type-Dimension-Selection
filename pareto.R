
install.packages("VGAM")
library(VGAM)
rows <- 10000
matrix_data <- matrix(NA, nrow = rows, ncol = 2)


set.seed(999)  
matrix_data[,1] <- rpareto(rows, shape = 1, scale = 1)
matrix_data[,2] <- rpareto(rows, shape = 1, scale = 1)

folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste0(folder_path, "pareto_10000obs.csv")
write.csv(matrix_data, file = file_name)
