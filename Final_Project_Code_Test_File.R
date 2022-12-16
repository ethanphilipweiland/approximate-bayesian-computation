library(testthat)
library(devtools)

#Note to John or Professor Fukuyama : 
# You may try to run the source_url function (line 9) to read the file
# If it doesn't work, you might just need to download the file and use source function to read it locally (line 8)

source('Final_Project_Code.R')
#source_url('https://github.com/ethanphilipweiland/STAT-S-610-Final-Project/blob/main/Final_Project_Code.R?raw=TRUE')

##
#Test 1: testing the probability matrix
##

#generate a few sample matrices
matrix1<-generate_probability_matrix(0.1,0.1,table_2_outbreak_1)

matrix2<-generate_probability_matrix(0.2,0.2,table_2_outbreak_2)

matrix3<-generate_probability_matrix(0.3,0.3,table_3_outbreak_1)

matrix4<-generate_probability_matrix(0.4,0.4,table_3_outbreak_2)


#function 1.1
#counting the number of NAs on the diagonal
matrix_diag_count<-function(matrix) {
  count<-sum(is.na(diag(matrix)))
  return(count)
}

#function 1.2
#counting the number of values that fall outside of [0,1]
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
#testing if the number of NAs on the diagonal is equal to 0
#passed
test_that("matrix diagonal has no NA values", {
  expect_equal(matrix_diag_count(matrix1), 0)
  expect_equal(matrix_diag_count(matrix2), 0)
  expect_equal(matrix_diag_count(matrix3), 0)
  expect_equal(matrix_diag_count(matrix4), 0)
})


#test 1.2
#testing if the number of values that fall outside of [0,1] is equal to 0
#passed
test_that("upper-triangle matrix values are in [0,1]", {
  expect_equal(matrix_values(matrix1), 0)
  expect_equal(matrix_values(matrix2), 0)
  expect_equal(matrix_values(matrix3), 0)
  expect_equal(matrix_values(matrix4), 0)
})



##
#Test 2: testing the distance function
##
#generate a few sample matrices

matrix5 <- matrix(rep(3, 9), nrow = 3, ncol = 3)

matrix6 <- matrix(rep(0, 9), nrow = 3, ncol = 3)

matrix7 <- matrix(rep(1, 9), nrow = 3, ncol = 3)

matrix8 <-matrix(rep(2, 9), nrow = 3, ncol = 3)


matrix9 <- matrix(rep(4, 9), nrow = 3, ncol = 3)

matrix10 <- matrix(rep(5, 9), nrow = 3, ncol = 3)

matrix11 <- matrix(rep(6, 9), nrow = 3, ncol = 3)

matrix12 <-matrix(rep(7, 9), nrow = 3, ncol = 3)


#test
#testing if the distance function works as designed 
#passed 
test_that("distance matches", {
  expect_equal(generate_distance(matrix5 , matrix6, matrix7, matrix8), 6)
  expect_equal(generate_distance(matrix9 , matrix10, matrix11, matrix12), 3)
})
