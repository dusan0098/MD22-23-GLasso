install.packages("mgm")
#install.packages("htmltools")
#install.packages("xfun")
library("mgm")

#Sampling a Pair-wise mixed model
#DOCUMENTATION https://rdrr.io/cran/mgm/man/mgmsampler.html
# --------- Example 3: p = 5 Mixed Graphical Model with two 2-way interaction ---------
# ----- 1) Specify Model -----

# a) General Graph Info - 5 gaussian, 5 categorical
type = c("g","g","g","g","g","c","c","c","c","c")
level = c(1, 1, 1, 1, 1,      3,  5,  6,  2,  4)
# b) Define Interaction
factors <- list()
#Pairwise interactions
factors[[1]] <- matrix(c(1,6, #Edge C-D 3 levels
                         1,10,#Edge C-D 4 levels
                         2,3,#Edge C-C
                         2,6,#Edge C-D 3 levels
                         3,7,#Edge C-D 5 levels
                         4,6,#Edge C-D 3 levels
                         4,7,#Edge C-D 5 levels
                         4,9,#Edge C-D 2 levels
                         4,10,#Edge C-D 4 levels
                         5,10,#Edge C-D 4 levels
                         7,8), #Edge D 5 levels-D 6 levels
                        ncol=2, byrow = T)
#Defining weights of pair-wise interactions
#First element of interactions has to be weights of all pair-wise interactions
interactions <- list()
interactions[[1]] <- vector("list", length = length(factors[[1]])/2)

#For each edge we create a scalar/vector/matrix of weights
#For now set all to "Betas,Rhos,Phis" to 0.5
for (i in 1:length(interactions[[1]])){
  edge <- factors[[1]][i,]
  if (i <= 6) {
    #We make edges 1-5 stronger interactions in order to 
    #Observe a difference in the conditional distributions after sampling
    interactions[[1]][[i]] <- array(0.1, dim=c(level[edge[1]],level[edge[2]]))
  } else {
    interactions[[1]][[i]] <- array(0.2, dim=c(level[edge[1]],level[edge[2]]))
  }

}
interactions[[1]][[6]] = interactions[[1]][[6]] + c(0.1,0.3,1)
interactions[[1]][[10]] = interactions[[1]][[10]] + c(0,0,0,0.5)
#Define Thresholds - for determining if interaction exists
thresholds <- list()
for (i in 1:length(type)){
  thresholds[[i]] <- rep(0, level[i])
}

#Define Variances
sds <- rep(.1, length(type))


# ----- 2) Sample cases -----

set.seed(1)
data <- mgmsampler(factors = factors,
                   interactions = interactions,
                   thresholds = thresholds,
                   sds = sds,
                   type = type,
                   level = level,
                   N = 500,
                   nIter = 100,
                   pbar = TRUE)


# ----- 3) Check: Conditional Means -----

# We condition on the categorical variables and check whether
# the conditional means match what we expect from the model:

dat <- data$data
sample=as.data.frame(dat)
sample

# Check interactions for edge 1-6
mean(dat[, 4])
#Should be 0 difference from overall mean
mean(dat[dat[,8] == 1, 4])
mean(dat[dat[,8] == 2, 4])
mean(dat[dat[,8] == 3, 4])
mean(dat[dat[,8] == 4, 4])
mean(dat[dat[,8] == 5, 4])
mean(dat[dat[,8] == 6, 4])

#Should be a difference from overall mean
mean(dat[dat[,6] == 1, 4]) 
mean(dat[dat[,6] == 2, 4])
mean(dat[dat[,6] == 3, 4])

#Edge 5-10
#Should be a difference from overall mean
mean(dat[, 5]) 
mean(dat[dat[,10] == 1, 5]) 
mean(dat[dat[,10] == 2, 5])
mean(dat[dat[,10] == 3, 5])
mean(dat[dat[,10] == 4, 5])

mgm_obj <- mgm(data = data$data,
               type = type,
               level = level,
               k = 2,
               lambdaSel = "CV",
               ruleReg = "AND",
               pbar = TRUE, 
               overparameterize = TRUE, 
               signInfo = FALSE)

library(qgraph)
qgraph(mgm_obj$pairwise$wadj,
       edge.color = mgm_obj$pairwise$edgecolor,
       layout = "spring",
       labels =  data$colnames)

mgm_obj$nodemodels[[2]]$model
############################################

