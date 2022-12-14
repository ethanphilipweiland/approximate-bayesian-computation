library(dplyr)
set.seed(911)

#Inserting observed data tables

  #Table 2: Influenza A (H3N2) infection in 1977-78 (middle column) and 1980-81 (right column) epidemics, Tecumseh, 
  #Michigan 
  row_0 <- c(66, 87, 25, 22, 4, 44, 62, 47, 38, 9)
  row_1 <- c(13, 14, 15, 9, 4, 10, 13, 8, 11, 5)
  row_2 <- c(NA, 4, 4, 9, 1, NA, 9, 2, 7, 3)
  row_3 <- c(NA, NA, 4, 3, 1, NA, NA, 3, 5, 1)
  row_4 <- c(NA, NA, NA, 1, 1, NA, NA, NA, 1, 0)
  row_5 <- c(NA, NA, NA, NA, 0, NA, NA, NA, NA, 1)
  Table_2 <- rbind(row_0, row_1, row_2, row_3, row_4, row_5)
  row.names(Table_2) <- c("0", "1", "2", "3", "4", "5")
  colnames(Table_2) <- c("1", "2", "3", "4", "5", "1", "2", "3", "4", "5")

  #Table 3: Influenza B infection in 1975-76 epidemic (middle column) and influenza A (H1N1) infection in 1978-79 
  #epidemic (right column), Seattle, Washington
  row_0 <- c(9, 12, 18, 9, 4, 15, 12, 4)
  row_1 <- c(1, 6, 6, 4, 3, 11, 17, 4)
  row_2 <- c(NA, 2, 3, 4, 0, NA, 21, 4)
  row_3 <- c(NA, NA, 1, 3, 2, NA, NA, 5)
  row_4 <- c(NA, NA, NA, 0, 0, NA, NA, NA)
  row_5 <- c(NA, NA, NA, NA, 0, NA, NA, NA)
  Table_3 <- rbind(row_0, row_1, row_2, row_3, row_4, row_5)
  row.names(Table_3) <- c("0", "1", "2", "3", "4", "5")
  colnames(Table_3) <- c("1", "2", "3", "4", "5", "1", "2", "3")
  
  #table_2_outbreak_1
  table_2_outbreak_1 <- Table_2[,1:5]
  
  #table_2_outbreak_2
  table_2_outbreak_2 <- Table_2[,6:10]
  
  #table_3_outbreak_1
  table_3_outbreak_1 <- Table_3[,1:5]
  
  #table_3_outbreak_2
  table_3_outbreak_2 <- Table_3[1:4,6:8]
  
#Probability matrix
generate_probability_matrix <- function(q_c, q_h, observed_data) {
  probability_matrix <- matrix(nrow=nrow(observed_data), ncol=ncol(observed_data)) #Probability matrix should 
  #be the same size as the observed data we are trying to simulate, varies by outbreak
  #The first row of the probability matrix (the probability that 0 out of the s susceptible individuals in a household
  #become infected) is given by w_0_s = (q_c)^s
  for (s in 1:ncol(probability_matrix)) {
    probability_matrix[1,s] <- (q_c)^s
  }
  #For the rest of the rows, the probability that j out of the s susceptibles in a household become infected is given
  #by the formula on pg. 107
  for (j in 1:(nrow(probability_matrix)-1)) { #j=1 occurs at row 2
    for (s in j:ncol(probability_matrix)) { #Number of susceptible individuals can't be less than number infected (logically)
      #w_j_j is equal to 1 minus the sum of the elements of column j from the element in the first row to the element
      #in the jth-1 row (will be from the first row to the jth row in the probability matrix because 
      #this matrix starts at row 1 rather than row 0 like in the observed data matrix)
      w_j_j <- 1 - sum(probability_matrix[1:j,j])
      probability_matrix[j+1,s] <- choose(s,j) * w_j_j * (q_c * (q_h)^j)^(s-j) #Formula from pg. 107
    }
  }
  #Informal and basic test to ensure the probability matrix is constructed correctly: Does the diagonal contain only real 
  #numbers, like in the observed data matrix?
  for (i in diag(probability_matrix)) {
    if (is.na(i)) {
      print("Warning! The probability matrix contains NAs along the diagonal")
      break
    }
  }
  return(probability_matrix)
}

#Data generating mechanism
generate_data <- function(probability_matrix, observed_data) {
  data <- matrix(nrow=nrow(observed_data), ncol=ncol(observed_data)) #Generated data should be the same size as the observed data
  #we are comparing it to, varies by outbreak
  for (i in 1:ncol(data)) {
    n <- 1 #We want one draw
    size <- sum(observed_data[,i], na.rm=T) #We obtain the total number of cases from the observed data
    prob <- probability_matrix[1:(i+1),i] #Obtaining the corresponding probabilities from the probability matrix
    draws <- rmultinom(n=n, size=size, prob=prob) #Simulating using a multinomial distribution
    data[1:(i+1), i] <- draws
  }
  return(data)
}

#Figure 3a.) and 3c.)
generate_figure <- function(observed_data_1, observed_data_2, N, tolerance) {
  #Creating empty vectors to store results
  q_c_1_vector <- c()
  q_c_2_vector <- c()
  q_h_1_vector <- c()
  q_h_2_vector <- c()
  model_vector <- c()
  for (i in 1:N) {
    #Simulating the first outbreak. Will be the same for both the two parameter and four parameter models,
    #as these models only vary in terms of q_c_2 and q_h_2
      #Drawing parameters from a uniform prior distribution
      q_c_1 <- runif(n=1)
      q_h_1 <- runif(n=1)
      #Generating the probability matrix based off these parameters
      probability_matrix_1 <- generate_probability_matrix(q_c=q_c_1, q_h=q_h_1, observed_data=observed_data_1)
      #Simulating data based off this probability matrix
      simulated_data_1 <- generate_data(probability_matrix=probability_matrix_1, observed_data=observed_data_1)
    #Simulating the second outbreak. In the two parameter model q_c and q_h will remain constant, while in the four parameter model they will 
    #be redrawn
      pick_model <- rbinom(n=1, size=1, prob=.5)
      #Four parameter model
        if (pick_model == 0) { 
          #Drawing parameters from a uniform prior distribution
          q_c_2 <- runif(n=1)
          q_h_2 <- runif(n=1)
          #Generating the probability matrix based off these parameters
          probability_matrix_2 <- generate_probability_matrix(q_c=q_c_2, q_h=q_h_2, observed_data=observed_data_2)
          #Simulating data based off this probability matrix
          simulated_data_2 <- generate_data(probability_matrix=probability_matrix_2, observed_data=observed_data_2)
          #Calculating distance and deciding whether to keep based off chosen tolerance level
          distance <- SUMMARY_STATISTIC_FUNCTION(observed_data_1, simulated_data_1, observed_data_2, simulated_data_2)
          if (distance <= tolerance) {
            q_c_2_vector[i] <- q_c_2
            q_h_2_vector[i] <- q_h_2
            model_vector[i] <- "Four Parameter"
          }
        }
      #Two parameter model
        if (pick_model == 1) {
          #Parameters remain constant
          q_c_2 <- q_c_1
          q_h_2 <- q_h_1
          #Generating the probability matrix based off these parameters
          probability_matrix_2 <- generate_probability_matrix(q_c=q_c_2, q_h=q_h_2, observed_data=observed_data_2)
          #Simulating data based off this probability matrix
          simulated_data_2 <- generate_data(probability_matrix=probability_matrix_2, observed_data=observed_data_2)
          #Calculating distance and deciding whether to keep based off chosen tolerance level
          distance <- SUMMARY_STATISTIC_FUNCTION(observed_data_1, simulated_data_1, observed_data_2, simulated_data_2)
          if (distance <= tolerance) {
            q_c_2_vector[i] <- q_c_2
            q_h_2_vector[i] <- q_h_2
            model_vector[i] <- "Two Parameter"
          }
        }
    }
  #NAs are in the storage vectors for iterations when the sample was rejected, fixing this:
    q_c_2_vector <- na.omit(q_c_2)
    q_h_2_vector <- na.omit(q_h_2)
    model_vector <- na.omit(model_vector)
  #Repeating q_c_1 and q_h_1 the number of times the sample was accepted as they are constant. Everything will be matched by position
    q_c_1_vector <- rep(q_c_1, length(q_c_2_vector))
    q_h_1_vector <- rep(q_h_1, length(q_h_2_vector))
  #Combining the vectors into a data frame (they are all the same length due to the last two steps)
    accepted_samples <- data.frame(q_c_1_vector, q_h_1_vector, q_c_2_vector, q_h_2_vector, model_vector)
  #Creating Figure 3
    plot <- ggplot(accepted_samples, aes(x=q_h_1_vector, y=q_c_1_vector, color=model_vector)) +
              geom_point(x=q_h_2_vector, y=q_c_2_vector, color=model_vector) +
              xlim(0,1) +
              ylim(0,1) +
              xlab("q_h") +
              ylab("q_c") +
              ggtitle("Figure 3")
    return(plot)
}
