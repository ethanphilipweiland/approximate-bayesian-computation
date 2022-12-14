library(testthat)

##
#Test 1
##
matrix1<-generate_probability_matrix(0.1,0.1,table_2_outbreak_1)

matrix2<-generate_probability_matrix(0.1,0.1,table_2_outbreak_2)

matrix3<-generate_probability_matrix(0.1,0.1,table_3_outbreak_1)

matrix4<-generate_probability_matrix(0.1,0.1,table_3_outbreak_2)

for (rows in 1:nrow(matrix1)) {
  for (columns in 1:ncol(matrix1)) {
    print(matrix1[rows, columns])
  }
}

length(matrix1)
is.na(matrix1[[7]])

#function 1.1
matrix_diag_count<-function(matrix) {
  count<-sum(is.na(diag(matrix)))
  return(count)
}

#function 1.2
matrix_values<- function(matrix) {
  c<-0
  for (i in 1:length(matrix)){
    if (is.na(matrix[[i]])==0){
      if (matrix[[i]] >=0 & matrix[[i]] <=1) {
        c<-c}
      else {
        c<-c+1
      }
    }
  }
  return(c)
}


#test 1.1
test_that("matrix diagonal has no NA values", {
  expect_equal(matrix_diag_count(matrix1), 0)
  expect_equal(matrix_diag_count(matrix2), 0)
  expect_equal(matrix_diag_count(matrix3), 0)
  expect_equal(matrix_diag_count(matrix4), 0)
})


#test 1.2
test_that("upper-triangle matrix values are in [0,1]", {
  expect_equal(matrix_values(matrix1), 0)
  expect_equal(matrix_values(matrix2), 0)
  expect_equal(matrix_values(matrix3), 0)
  expect_equal(matrix_values(matrix4), 0)
})



##
#Test 2
##
data1 <- matrix(rep(3, 9), nrow = 3, ncol = 3)

data2 <- matrix(rep(0, 9), nrow = 3, ncol = 3)

data3 <- matrix(rep(NA, 9), nrow = 3, ncol = 3)

#function
distance_fun <- function(observed_data, simulated_data){
  observed_data[is.na(observed_data)] <- 0
  simulated_data[is.na(simulated_data)] <- 0
  new_matrix <- as.matrix(observed_data - simulated_data)
  #distance <- sqrt(sum(new_matrix^2))
  distance <- norm(new_matrix, type = "F")
  return(distance)
}

#test
test_that("distance matches", {
  expect_equal(distance_fun(data1, data1), 0)
  expect_equal(distance_fun(data1, data2), 9)
  expect_equal(distance_fun(data1, data3), 9)
})


