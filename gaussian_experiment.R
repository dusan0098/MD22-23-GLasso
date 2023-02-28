#synthetic data experiments for validating
#Friedmans glasso package 

#1. Methods for constructing Covariance matrix
#2. Sampling methods
#3. Attempt at recovery of the full sparsity structure

install.packages("glasso")
install.packages("gRbase")
install.packages("mvtnorm")
library(pcalg)
library(glasso)
library(igraph)
library(dplyr)
library(gRbase)
library("mvtnorm")
#set directory so loading data works
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

sample_size <- 150
dimension <- 400

#Precision matrix for sparse model
#Diagonal = 1, Lower,upper diagonal =0.5, others = 0
create_sparse_matrix <- function(dimension){
  mat <- diag(dimension)
  for (i in 1:nrow(mat)){
    if(i-1>=1){
      mat[i,i-1]=0.5
    }
    if(i-1>=1){
      mat[i-1,i]=0.5
    }
  }
  #print(mat)
  return (solve(mat))
}

#Precision matrix for dense model
#Diagonal=2, all others =1
create_dense_matrix <- function(dimension){
  mat <- matrix(1, ncol=dimension, nrow=dimension)
  mat <- mat + diag(dimension)
  #print(mat)
  return(solve(mat))
}

#Creates a sample for the dense/sparse scenario
create_gaussian_sample <- function (sample_size,dimension, sparse=TRUE){
  if(sparse){
    mat <-create_sparse_matrix(dimension)
  }else{
    mat <- create_dense_matrix(dimension)
  }
  sample <- rmvnorm(sample_size, sigma = mat)
  return(sample)
}

#Sparse scenario
sample_sparse <- create_gaussian_sample(sample_size, dimension)
sample_dense <- create_gaussian_sample(sample_size, dimension, sparse = FALSE)
sparse_estimate <-glasso(s = var(sample_sparse), 
                          rho = 4,
                          nobs = sample_size,
                         maxit = 30)
sparse_estimate['wi']
colSums(matrix(unlist(sparse_estimate['wi']),nrow=dimension, byrow=TRUE)!=0)

#Dense scenario
dense_estimate <-glasso(s = var(sample_dense), 
                         rho = 0.02,
                         nobs = sample_size,
                        maxit = 30)
dense_estimate['wi']
colSums(matrix(unlist(dense_estimate['wi']),nrow=dimension, byrow=TRUE)!=0)

