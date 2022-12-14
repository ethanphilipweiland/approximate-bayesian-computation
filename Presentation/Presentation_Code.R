library(dplyr)

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
generate_probability_matrix <- function(q_c, q_h, observed_data_matrix) {
  probability_matrix <- matrix(nrow=nrow(observed_data_matrix), ncol=ncol(observed_data_matrix)) #Probability matrix should 
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
generate_data <- function(probability_matrix, observed_data_matrix) {
  data <- matrix(nrow=nrow(observed_data_matrix), ncol=ncol(observed_data_matrix)) #Generated data matrix should be the same
  #size as the observed data we are comparing it to, varies by outbreak
  for (i in 1:nrow(probability_matrix)) {
    for (j in 1:ncol(probability_matrix)) {
      if (!is.na(probability_matrix[i,j])) { #For every non-NA element in the probability matrix
        if (sum(!is.na(probability_matrix[,j])) == 2) { #If there are two non-NA values in the column containing
        #the element [i,j]
          size <- sum(observed_data_matrix[,j], na.rm=T)
          data[i,j] <- rbinom(n=1, size=size, p=probability_matrix[i,j]) #Since there are two non-NA elements in the column, we 
          #sample from a binomial distribution
        } else { #Otherwise we use a multinomial distribution
          size <- sum(observed_data_matrix[,j], na.rm=T)
          data[i,j] <- rmultinom(n=1, size=size, p=probability_matrix[i,j])
        }
      }
    }
  }
  return(data)
}



# primary goal of this function: compute a similarity measure between the simulated data and an observed data
# and keep or discard parameter samples based on the simularity

# similarity measure: the Frobenius norm (also called Euclidean norm)
# 
# C = A - B
# || C || = sqrt ( sum ( (C_ij)^2 ) )


# distance function (summary function)
distance_fun <- function(observed_data, simulated_data){
  observed_data[is.na(observed_data)] <- 0
  simulated_data[is.na(simulated_data)] <- 0
  new_matrix <- as.matrix(observed_data - simulated_data)
  #distance <- sqrt(sum(new_matrix^2))
  distance <- norm(new_matrix, type = "F")
  return(distance)
}

# informal testing
data1 <- matrix(rep(3, 9), nrow = 3, ncol = 3)
data1
# distance between a dataset and itself should be 0 (suggests high similarity)
distance_fun(data1, data1)

data2 <- matrix(rep(0, 9), nrow = 3, ncol = 3)
data2
# distance between data1 and data2 should be 9
distance_fun(data1, data2)


data3 <- matrix(rep(NA, 9), nrow = 3, ncol = 3)
data3
# distance between data1 and data3 should be 9
distance_fun(data1, data3)

data4 <- matrix(c(1:9), nrow = 3, ncol = 3)
data4[2, 1] <- NA
data4[3, 1] <- NA
data4[3, 2] <- NA
data4
# distance between data1 and data4 should be about 10.63
distance_fun(data1, data4)


# ABC algorithm
# observed_data = used for simulation and the comparison between simulated one and observed one
# summary_statistic = a function to give us approximate samples from the posterior instead of exact samples
# data_generating_function = function used for data simulation from p(y_obs | parameters)
# epsilon  = (tolerance); the  desired level of agreement between observed_data and simulated_data
generate_abc_sample <- function(observed_data,
                                summary_statistic,
                                prior_distribution,
                                data_generating_function,
                                epsilon) {
  while(TRUE) {
    theta <- prior_distribution()
    simulated_data <- data_generating_function(theta)
    if(distance_fun(observed_data, simulated_data) <= epsilon) {
      return(theta)
    }
  }
}

posterior_samples <- replicate(n = 1000, 
                               generate_abc_sample(observed_data, summary_statistic,
                                                   prior_distribution, data_generating_function,
                                                   epsilon))
hist(posterior_samples)