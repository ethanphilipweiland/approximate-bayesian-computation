library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyverse)
set.seed(2022)



# ABC algorithm
# observed_data = used for simulation and the comparison between simulated one and observed one
# summary_statistic = a function to give us approximate samples from the posterior instead of exact samples
# data_generating_function = function used for data simulation from p(y_obs | parameters)
# epsilon  = (tolerance); the  desired level of agreement between observed_data and simulated_data
generate_abc_sample <- function(observed_data1, 
                                observed_data2,
                                summary_statistic,
                                prior_distribution,
                                data_generating_function,
                                epsilon) {
  while(TRUE) {
    q_c_1 <- prior_distribution()
    q_h_1 <- prior_distribution()
    q_c_2 <- prior_distribution()
    q_h_2 <- prior_distribution()
    
    prob_matrix_1 <- generate_probability_matrix(q_c_1, q_h_1, observed_data1)
    prob_matrix_2 <- generate_probability_matrix(q_c_2, q_h_2, observed_data2)
    
    simulated_data1 <- data_generating_function(prob_matrix_1, observed_data1)
    simulated_data2 <- data_generating_function(prob_matrix_2, observed_data2)
    
    if(summary_statistic(observed_data1, simulated_data1, observed_data2, simulated_data2) <= epsilon) {
      theta = c(q_c_1, q_h_1, q_c_2, q_h_2)
      return(theta)
    }
  }
}


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


#Distance function
generate_distance <- function(observed_data_1, simulated_data_1, observed_data_2, simulated_data_2) {
  #First, a function to calculate the Frobenious norm 
  #The Frobenious norm is matrix norm of an matrix defined as the square root of the sum of the absolute squares of its elements
  frobenious <- function(matrix) {
    frobenious_norm <- sum(abs(matrix)^2, na.rm=T) %>% sqrt()
    return(frobenious_norm)
  }
  #Formula for the distance function is given on pg. 107
  distance <- (1/2) * (frobenious(observed_data_1 - simulated_data_1) + frobenious(observed_data_2 - simulated_data_2))
  return(distance)
}

prior_distribution <- function() runif(n=1)
summary_statistic <- generate_distance
data_generating_function <- generate_data






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


#Table 2
# sample from posterior ditribution for the different outbreaks of the same strain
posterior_samples <- replicate(n = 1000, 
                               generate_abc_sample(table_2_outbreak_1,
                                                   table_2_outbreak_2,
                                                   summary_statistic,
                                                   prior_distribution, 
                                                   data_generating_function,
                                                   epsilon = 25))

# convert the parameter particles into appropriate data frame
m = as.data.frame(posterior_samples)
row.names(m) = c("q_c_1", "q_h_1", "q_c_2", "q_h_2")


z = data.frame(q_c_1 = t(m[1,]),
               q_h_1 = t(m[2,]),
               q_c_2 = t(m[3,]),
               q_h_2 = t(m[4,]))


#Table3
# sample from posterior ditribution for different outbreaks of different strains
posterior_samples2 <- replicate(n = 1000, 
                                generate_abc_sample(table_3_outbreak_1,
                                                    table_3_outbreak_2,
                                                    summary_statistic,
                                                    prior_distribution, 
                                                    data_generating_function,
                                                    epsilon = 10))


# convert the parameter particles into appropriate data frame
m = as.data.frame(posterior_samples)
row.names(m) = c("q_c_1", "q_h_1", "q_c_2", "q_h_2")

m2 = as.data.frame(posterior_samples2)
row.names(m2) = c("q_c_1", "q_h_1", "q_c_2", "q_h_2")


z2 = data.frame(q_c_1 = t(m2[1,]),
                q_h_1 = t(m2[2,]),
                q_c_2 = t(m2[3,]),
                q_h_2 = t(m2[4,]))


# generate figures
# figure 3a
ggplot(z)+
  geom_point(aes(x=q_h_1, y=q_c_1), color = 'red3') +
  geom_point(aes(x=q_h_2, y=q_c_2), color = 'blue3') +
  xlim(0,1) +
  ylim(0,1) +
  xlab("q_h") +
  ylab("q_c") +
  ggtitle("Figure 3a") +
  theme_few()


#Figure 3c)
ggplot(z2)+
  geom_point(aes(x=q_h_1, y=q_c_1), colour = 'red3') +
  geom_point(aes(x=q_h_2, y=q_c_2), colour = 'blue3') +
  xlim(0,1) +
  ylim(0,1) +
  xlab("q_h") +
  ylab("q_c") +
  ggtitle("Figure 3c") +
  theme_few()
