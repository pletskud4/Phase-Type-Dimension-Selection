library(matrixdist)
options(digits = 22)
# Load data
firedata <- scan(file = "C:/Users/bjorn/OneDrive/Desktop/speciale/Determination of n/data/firedata.dat")

# View histogram og kør standard
h <- hist(firedata, breaks = 30, probability = TRUE, col = "darkgreen")

ph_model <- ph(structure = "general", dimension = 8)
fit <- fit(ph_model, y = firedata, stepsEM = 1)
fit <- fit(fit, y = firedata, stepsEM = 1000, every = 5)
lines(seq(0,10, 0.01), dens(fit, seq(0,10, 0.01)), col = "blue")

#prøver også lidt forskelligt med diverse dimensioner
exp_data <- rexp(1, rate = 0.5)
ph_model <- ph(structure = "general", dimension = 2)
fit <- fit(ph_model, y = exp_data, stepsEM = 1)
fit <- fit(fit, y = exp_data, stepsEM = 500, every = 1)
lines(seq(0,10, 0.01), dens(fit, seq(0,10, 0.01)), col = "blue")

# Store log-likelihood values in a matrix
#denne tager 5 timer at køre
logLik_matrix <- matrix(NA, nrow = 1000, ncol = 10)  # 200 rows (from 5 to 1000 by 5) and 20 columns (dimensions 1 to 20)

for (dim in 1:10) {
  ph_model <- ph(structure = "general", dimension = dim)
  fit <- fit(ph_model, y = firedata, stepsEM = 1)
  
  for (iter in seq(1, 1000, by = 1)) {
    fit <- fit(fit, y = firedata, stepsEM = iter, every = 1)
    logLik <- logLik(fit)
    logLik_matrix[iter, dim] <- logLik
  }
}
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Determination of n/")
write.table(logLik_matrix, file = "firedatatendimensions.txt", sep = "\t", row.names = FALSE)
 
#denne vil jeg nu kunne læse når som helst og lave et plot til
file_path <- "firedatatendimensions.txt"

# Read structured data from the text file (e.g., tab-separated or comma-separated)
logLik_matrix <- read.table(file_path, header = TRUE)
#tilføjer titler
row_headers <- paste(seq(from = 1, to = 1000, by = 1))
col_headers <- paste("dim", seq(from = 1, to = 10))
rownames(logLik_matrix) <- row_headers
colnames(logLik_matrix) <- col_headers
matplot(logLik_matrix, type = "p", pch = 19, col = rainbow(20),cex = 0.3,
        xlab = "Iterations", ylab = "Log-likelihood",
        main = "Scatter plot of achieved log likelihood", xaxt = "n")
axis(1, at = seq(1, nrow(logLik_matrix)), labels = seq(1, 1000, by = 1))
legend("topright", legend = colnames(logLik_matrix), col = rainbow(20), pch = 19, cex = 0.38)

#development of matrix:
changematrix <- logLik_matrix[2:60,] - logLik_matrix[1:59,]
matplot(changematrix, type = "l", pch = 19, col = rainbow(20, alpha = 0.7), lty = 1, lwd = 2,
        xlab = "Iterations", ylab = "Log-likelihood change", xaxt = "n")
axis(1, at = seq(1, nrow(changematrix)), labels = seq(2, 60, by = 1))
legend("topright", legend = colnames(logLik_matrix), col = rainbow(20, alpha = 0.7), pch = 19, cex = 0.8)
max_values <- apply(changematrix, 2, max)
max_values
apply(changematrix, 2, sum)
# Initialize a new matrix with NAs
new_matrix <- matrix(NA, nrow = nrow(changematrix), ncol = ncol(changematrix))
k <- 1
new_row <- c(NA, NA, NA)
maxima <- matrix(NA, nrow = 1, ncol = 3)
# Loop through each cell and check the condition
for (i in 2:(nrow(changematrix) - 1)) {
  for (j in 1:ncol(changematrix)) {
    if (changematrix[i - 1, j] < changematrix[i, j] && changematrix[i, j] > changematrix[i + 1, j]) {
      new_matrix[i, j] <- changematrix[i, j]
      if (k == 1) {
        maxima <- c(i, j, new_matrix[i, j])
      } else {
        maxima <- rbind(maxima, c(i, j, new_matrix[i, j]))
      }
      k = k+1
    }
  }
}
plot(maxima[, 1], maxima[, 3], pch = 19, col = rainbow(length(unique(maxima[, 2])))[as.numeric(factor(maxima[, 2]))],
     xlab = "X Axis", ylab = "Y Axis", main = "Scatterplot with Colored Points")
legend("topright", legend = levels(factor(maxima[, 2])),
       fill = rainbow(length(unique(maxima[, 2]))), title = "Values",cex=0.8)




#hernæst vil vi lave histogram plot med alle dimensionerne
# View histogram
h <- hist(firedata, breaks = 30, probability = TRUE, col = "darkgreen", xlab = "firedata", ylab = "density")

# Loop through dimensions from 1 to 20 and add lines to the existing plot
for (dim in 1:20) {
  ph_model <- ph(structure = "general", dimension = dim)
  fit <- fit(ph_model, y = firedata, stepsEM = 1)
  fit <- fit(fit, y = firedata, stepsEM = 1000, every = 5)
  
  lines(seq(0, 12, 0.01), dens(fit, seq(0, 12, 0.01)), col = rainbow(20)[dim])
}

# Add legend with smaller text and modified labels
legend("topright", legend = paste0("dim ", 1:20), col = rainbow(20), lty = 1, cex = 0.38)




#find første gang log-likelihood værdier er præcis ens
check_convergence <- function(matrix_data) {
  convergence_vector <- numeric(ncol(matrix_data))
  for (i in 1:ncol(matrix_data)) {
    for (j in 2:nrow(matrix_data)) {
      if (matrix_data[j, i] == matrix_data[j-1, i]) {
        convergence_vector[i] <- (j - 1) * 5  # Noting iteration
        break
      } else if (j == nrow(matrix_data)){
        convergence_vector[i] <- NA  # Not converged, mark as N/A
        break  # Move to the next column
      } else {
        next
      }
    }
  }
  return(convergence_vector)
}
# Checking convergence for each column in the matrix
convergence_result <- check_convergence(logLik_matrix)
table_values <- matrix(NA, nrow = 2, ncol = length(convergence_result) + 1)
table_values[1, 1] <- "Dimensions"
table_values[2, 1] <- "Iterations until convergence"
table_values[1, 1:length(convergence_result) + 1] <- paste(1:length(convergence_result))
table_values[2, 1:length(convergence_result) + 1] <- convergence_result
# printing the table nicely
table_values





#vil gerne se udviklingen i log-likelihood
difference_logLik_matrix <- matrix(NA, nrow = nrow(logLik_matrix) - 1, ncol = ncol(logLik_matrix))

# Calculate differences and store them in the new matrix
for (j in 1:(ncol(logLik_matrix))) {
  for (i in 1:(nrow(logLik_matrix) - 1)) {
    difference_logLik_matrix[i, j] <- logLik_matrix[i + 1, j] - logLik_matrix[i, j]
  }
}
#kun interessert i første 200 iterationer
difference_matrix <- difference_logLik_matrix[1:20, ]
#endelig vil vi gerne lave denne lidt pænere
row_headers <- paste(seq(from = 5, to = 100, by = 5))
col_headers <- paste("dim", seq(from = 1, to = 20))
rownames(difference_matrix) <- row_headers
colnames(difference_matrix) <- col_headers
matplot(difference_matrix, type = "p", pch = 19, col = rainbow(20),cex = 0.5,
        xlab = "Iterations", ylab = "Log-likelihood difference",
        main = "Scatter plot of achieved log likelihood", xaxt = "n")
axis(1, at = seq(1, nrow(difference_matrix)), labels = seq(5, 100, by = 5))
legend("topright", legend = colnames(logLik_matrix), col = rainbow(20), pch = 19, cex = 0.38)







#kører lige de første fem iterationer for alle dim
logLik_matrix <- matrix(NA, nrow = 5, ncol = 20)  # 200 rows (from 5 to 1000 by 5) and 20 columns (dimensions 1 to 20)
for (dim in 1:20) {
  ph_model <- ph(structure = "general", dimension = dim)
  fit <- fit(ph_model, y = firedata, stepsEM = 1)
  
  for (iter in seq(1, 5, by = 1)) {
    fit <- fit(fit, y = firedata, stepsEM = iter, every = 1)
    logLik <- logLik(fit)
    logLik_matrix[iter/1, dim] <- logLik
  }
}
#og vil printe dem 
row_headers <- paste(seq(from = 1, to = 5, by = 1))
col_headers <- paste("dim", seq(from = 1, to = 20))
rownames(logLik_matrix) <- row_headers
colnames(logLik_matrix) <- col_headers
matplot(logLik_matrix, type = "p", pch = 19, col = rainbow(20),cex = 0.6,
        xlab = "Iterations", ylab = "Log-likelihood",
        main = "Scatter plot of achieved log likelihood", xaxt = "n")
axis(1, at = seq(1, nrow(logLik_matrix)), labels = seq(1, 5, by = 1))
legend("topleft", legend = colnames(logLik_matrix), col = rainbow(20), pch = 19, cex = 0.3)


# Store log-likelihood values in a matrix
#samme som før men kører kun til hundrede iterationer (til gengæld alle trækkes ud)
logLik_matrix <- matrix(NA, nrow = 100, ncol = 20)  # 200 rows (from 5 to 1000 by 5) and 20 columns (dimensions 1 to 20)
for (dim in 1:20) {
  ph_model <- ph(structure = "general", dimension = dim)
  fit <- fit(ph_model, y = firedata, stepsEM = 1)
  
  for (iter in seq(1, 100, by = 1)) {
    fit <- fit(fit, y = firedata, stepsEM = iter, every = 1)
    logLik <- logLik(fit)
    logLik_matrix[iter, dim] <- logLik
  }
}
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Determination of n/")
write.table(logLik_matrix, file = "firedatadimensions2.txt", sep = "\t", row.names = FALSE)
difference_logLik_matrix <- matrix(NA, nrow = nrow(logLik_matrix) - 1, ncol = ncol(logLik_matrix))
# Calculate differences and store them in the new matrix
for (j in 1:(ncol(logLik_matrix))) {
  for (i in 1:(nrow(logLik_matrix) - 1)) {
    difference_logLik_matrix[i, j] <- logLik_matrix[i + 1, j] - logLik_matrix[i, j]
  }
}
row_headers <- paste(seq(from = 1, to = 99, by = 1))
col_headers <- paste("dim", seq(from = 1, to = 20))
rownames(difference_logLik_matrix) <- row_headers
colnames(difference_logLik_matrix) <- col_headers
matplot(difference_logLik_matrix, type = "p", pch = 19, col = rainbow(20),cex = 0.6,
        xlab = "Iterations", ylab = "Log-likelihood difference",
        main = "Scatter plot of achieved log likelihood", xaxt = "n")
axis(1, at = seq(1, nrow(difference_logLik_matrix)), labels = seq(1, 99, by = 1))
legend("topleft", legend = colnames(difference_logLik_matrix), col = rainbow(20), pch = 19, cex = 0.3)

#vil gerne lave en tabel med start log-lik, konvergens log-lik, osv
# Create an empty matrix to store the results
new_matrix <- matrix(NA, nrow = 3, ncol = ncol(difference_logLik_matrix))
colnames(new_matrix) <- colnames(difference_logLik_matrix)
# Loop through each column
for (col in 1:ncol(difference_logLik_matrix)) {
  # Get the first cell value in the column
  first_value <- difference_logLik_matrix[1, col]
  # Initialize variables to store row number and cell value
  row_number <- NA
  cell_value <- NA
  # Loop through rows starting from the second row
  for (row in 2:nrow(difference_logLik_matrix)) {
    # Check if the cell value is less than or equal to the first cell value
    if (difference_logLik_matrix[row, col] <= first_value) {
      # Note down row number and cell value
      row_number <- row
      cell_value <- difference_logLik_matrix[row, col]
      break  # Exit loop when condition is met
    }
  }
  # Store the information in the new matrix
  new_matrix[, col] <- c(row_number, cell_value, first_value)
}
#tilføjer data til denne matrix
additional_rows <- matrix(NA, nrow = 2, ncol = ncol(new_matrix))
colnames(additional_rows) <- colnames(new_matrix)
# Fill the first additional row
additional_rows[1, ] <- logLik_matrix[1, ]
# Fill the second additional row based on the first row of 'new_matrix'
for (i in 1:ncol(new_matrix)) {
  additional_rows[2, i] <- logLik_matrix[new_matrix[1, i], i]
}
# Combine 'new_matrix' and the additional rows
extended_matrix <- rbind(new_matrix, additional_rows)
# tilføjer også hvornår konvergens opnås
new_row <- table_values[2, 2:ncol(table_values)]
# Add the new row to the extended matrix
final_matrix <- rbind(extended_matrix, new_row)
# tilføjer til sidst konvergensværdi
additional_rows2 <- matrix(NA, nrow = 1, ncol = ncol(final_matrix))
for (i in 1:ncol(final_matrix)) {
  if (is.na(final_matrix[6, i])) {
    additional_rows2[1, i] <- NA
  } else {
    additional_rows2[1, i] <- logLik_matrix[as.numeric(final_matrix[6, i]) / 5, i]
  }
}
final_matrix2 <- rbind(final_matrix, additional_rows2)
#til sidst regner vi forskel på loglik, og gennemsnitlig udvikling pr iteration
#både før og efter vores forslag
calculated_rows <- matrix(NA, nrow = 4, ncol = ncol(final_matrix2))
calculated_rows[1, ] <- as.numeric(final_matrix2[5,])-as.numeric(final_matrix2[4,])
calculated_rows[2, ] <- as.numeric(final_matrix2[7,])-as.numeric(final_matrix2[5,])
calculated_rows[3, ] <- as.numeric(calculated_rows[1,])/as.numeric(final_matrix2[1,])
calculated_rows[4, ] <- as.numeric(calculated_rows[2,])/(as.numeric(final_matrix2[6,])-as.numeric(final_matrix2[1,]))
final_matrix3 <- rbind(final_matrix2, calculated_rows)
row_names <- c("Iterations until stopping criteria", "Development", "Development first iteration", "Initial loglik", "Stopping criteria loglik", "Convergence iteration", "Convergence value","Development until stopping criteria", "Development from stopping to convergence", "Average development (5 ite) until stopping", "Average development (5 ite) stopping until convergence")  # Replace these with your desired row names
# Assign row names to the matrix
rownames(final_matrix3) <- row_names
#vil se tabellen





#vil lige sikre at der kun er et toppunkt til hver
maxfinder <- matrix(NA, nrow = nrow(difference_logLik_matrix)-2, ncol = ncol(difference_logLik_matrix))
for (i in 1:ncol(maxfinder)) {
  for (j in 1:nrow(maxfinder)){
    if (as.numeric(difference_logLik_matrix[j+1, i])>as.numeric(difference_logLik_matrix[j, i]) && as.numeric(difference_logLik_matrix[j+1, i])>as.numeric(difference_logLik_matrix[j+2, i])) {
      maxfinder[j, i] <- 1 
    } else {
      maxfinder[j,i] <- 0
    }
  }
}
colSums(maxfinder)






#op til ti dimensioner tabel
#denne tager 5 timer at køre
logLik_matrix <- matrix(NA, nrow = 1000, ncol = 10)  # 200 rows (from 5 to 1000 by 5) and 20 columns (dimensions 1 to 20)

for (dim in 1:10) {
  ph_model <- ph(structure = "general", dimension = dim)
  fit <- fit(ph_model, y = firedata, stepsEM = 1)
  
  for (iter in seq(1, 1000, by = 1)) {
    fit <- fit(fit, y = firedata, stepsEM = iter, every = 1)
    logLik <- logLik(fit)
    logLik_matrix[iter, dim] <- logLik
  }
}
setwd("C:/Users/bjorn/OneDrive/Desktop/speciale/Determination of n/")
write.table(logLik_matrix, file = "firedatadimensions3.txt", sep = "\t", row.names = FALSE)
#denne vil jeg nu kunne læse når som helst og lave et plot til
file_path <- "firedatadimensions2.txt"
#laver nu en matrix med forskellene
#vil gerne se udviklingen i log-likelihood
tolerance <- 1e-12  # Define a tolerance threshold (you can adjust this value based on your requirement)
difference_logLik_matrix <- matrix(NA, nrow = nrow(logLik_matrix) - 1, ncol = ncol(logLik_matrix))
convergence_ite <- matrix(NA, nrow = 1, ncol = ncol(logLik_matrix))
# Calculate differences and store them in the new matrix
for (j in 1:(ncol(logLik_matrix))) {
  for (i in 1:(nrow(logLik_matrix) - 1)) {
    difference_logLik_matrix[i, j] <- logLik_matrix[i + 1, j] - logLik_matrix[i, j]
    if (abs(logLik_matrix[i + 1, j] - logLik_matrix[i, j]) <= tolerance && 
        abs(logLik_matrix[i, j] - logLik_matrix[max(i - 1, 1), j]) > tolerance) {
      convergence_ite[1, j] <- i
    }
  }
}

iterations <- numeric(length = ncol(difference_logLik_matrix))  # To store the number of iterations for each column
# Loop through each column of the matrix
for (col_index in 2:ncol(difference_logLik_matrix)) {
  # Loop backward starting from the vector value for the column
  for (i in convergence_ite[1,col_index]:2) {
    if (difference_logLik_matrix[i, col_index] > difference_logLik_matrix[i - 1, col_index]) {
      iterations[col_index] <- convergence_ite[1,col_index] - i + 1  # Note the number of iterations
      break  # Exit the loop once condition is met
    }
  }
}
column_to_plot <- difference_logLik_matrix[54:60, 10]  # Selecting the first column
# Plotting the values in the selected column
plot(column_to_plot, type = "p", xlab = "Index", ylab = "Values", main = "Plot of a Matrix Column")
file_path <- "firedatadimensions3.txt"

# Read structured data from the text file (e.g., tab-separated or comma-separated)
logLik_matrix <- read.table(file_path, header = TRUE)


