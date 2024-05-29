
library(matrixdist)
library(actuar)
rows <- 10
dim <- 10
matrix_data <- matrix(NA, nrow = rows, ncol = 2)


set.seed(999)  
ph_model1 <- ph(structure = "general", dimension = dim)
ph_model2 <- ph(structure = "general", dimension = dim)
matrix_data[,1] <- rphtype(rows, ph_model1@pars$alpha, ph_model1@pars$S)
matrix_data[,2] <- rphtype(rows, ph_model2@pars$alpha, ph_model2@pars$S)


folder_path <- "C:/Users/bjorn/OneDrive/Desktop/speciale/Final files/Simulated datasets/Datasets/"
file_name <- paste0(folder_path, "ph5_10000obs.csv")
write.csv(matrix_data, file = file_name)
