#Contains data and visualisations used in the implementation section
#of the Final Report

#1.Graphical Lasso
#2.PC algorithm
#3.Visualisations and comparisons
#4.Data - protein cell signaling


#packages for reading the data and visualisation
install.packages("readxl")
install.packages("glasso")
install.packages("gRbase")
library("readxl")
library(pcalg)
library(glasso)
library(igraph)
library(dplyr)
library(gRbase)
#set directory so loading data works
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#Protein data - used in Friedman 2008
files = c("./dataGLasso/1. cd3cd28.xls",
          "./dataGLasso/2. cd3cd28icam2.xls",
          "./dataGLasso/3. cd3cd28+aktinhib.xls",
          "./dataGLasso/4. cd3cd28+g0076.xls",
          "./dataGLasso/5. cd3cd28+psitect.xls",
          "./dataGLasso/6. cd3cd28+u0126.xls",
          "./dataGLasso/7. cd3cd28+ly.xls",
          "./dataGLasso/8. pma.xls",
          "./dataGLasso/9. b2camp.xls")

#Reading the files into a single dataframe
data <- read_excel("./dataGLasso/1. cd3cd28.xls",sheet = 1)
for(f in files[-1]){
  new_data <- read_excel(f,sheet = 1)
  print(new_data)
  data <- rbind(data, new_data)
  print(c("success:",f))
}
#For first 9 files should be 7466
print(nrow(data))

#Change column order and names to those from Friedman 2008
data <- data[, c(9,10,11,1,2,3,4,5,6,7,8)]
data
data <- data %>% 
  rename(Raf = praf,
         Jnk = pjnk ,
         Akt = pakts473,
         Erk = 'p44/42' ,
         Plcg = plcg,
         Mek = pmek )
#data

#Building graph for ground truth
create_ground_truth <- function(symmetric = FALSE){
  gr_truth <- matrix(0,11,11)
  colnames(gr_truth) <- colnames(data)
  rownames(gr_truth) <- colnames(data)
  #PKC edges
  gr_truth[1,2:5] = 1
  gr_truth[6,1] = 1
  gr_truth[7,1] = 1
  #P38 edges
  gr_truth[11,2] = 1
  #
  gr_truth[11,3] = 1
  #
  gr_truth[11,4] = 1
  gr_truth[4,5] = 1 
  #
  gr_truth[5,9] = 1
  gr_truth[11,5] = 1
  #
  gr_truth[6,7:8] = 1
  #
  gr_truth[8,7] = 1
  #
  gr_truth[8,10] = 1
  #
  gr_truth[9,10] = 1
  gr_truth[11,9] = 1
  #
  gr_truth[11,10] = 1
  #Forcing symmetric matrix
  if(symmetric){
    gr_truth[lower.tri(gr_truth)] = t(gr_truth)[lower.tri(gr_truth)]
  }
  return(gr_truth)
}
true_graph <- graph_from_adjacency_matrix(create_ground_truth(symmetric = FALSE),
                                          #Switched
                                          mode="directed")
true_graph
#Plot Friendmans ground truth graph
true_graph <- moralize(true_graph)
plot(true_graph , 
     edge.arrow.size=.5, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White", vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 1.6,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 1.5,
     main = "Ground truth")




#PC algorithm - returns PC_skeleton graphs for a list of alphas
pc_skeletons <- function(data, alphas){
  sufficient_statistics<- list(C = cor(data), n = nrow(data))
  column_names <- colnames(data) 
  skeleton_surveys <- sapply(alphas, pcalg::skeleton,
                      suffStat=sufficient_statistics,
                      indepTest= gaussCItest,
                      labels = column_names)
  #print(skeleton_surveys)
  skel_graphs <- list()
  for(s in skeleton_surveys){
    new_graph <- graph_from_graphnel(s@graph, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
    print(length(skel_graphs))
    skel_graphs[[length(skel_graphs)+1]] <- new_graph
  }
  
  return(skel_graphs)
}


#Difference of two graphs
difference_graph <- function(first, second, both=T){
  G1 <- first %m% second
  G2 <- second %m% first
  G3 <- G1 %u% G2
  if(both){
    return(G3)
  }
  return(G1)
}


#Ploting Graphical model requires matrix of 0s and 1s
threshold_wi <- function(wi){
  limit <- 0
  #wi[wi = limit] = 0
  wi[wi > limit] = 1
  wi[wi < limit] = 1
  diag(wi) <- 0
  return(wi)
}

#Performs GLasso for a list of lambdas and returns networks as list
networks_GLasso <- function(data, lambdas){
  networks <- list()
  S = var(data)
  
  for (lambda in lambdas){
    #Friedmans implementation of GLasso - documentation at:
    #https://cran.r-project.org/web/packages/glasso/glasso.pdf
    #S - empirical covariance
    #rho = "lambda" - penalisation, can also be a PxP matrix if we want different penalties per variable
    #nobs - number of observations
    #zero - indices to be fixed to 0 (known independence)
    #thr (10e-4) - minimal average change of the W matrix in a step to continue optimisation
    #maxit (10e4) - maximal iterations
    #approx - if true performs PAIR-WISE regression instead
    #start (cold) - if set to warm we can provide intial estimates for W and Wi
    inverse_estimate <-glasso(s = S, 
                              rho = lambda,
                              nobs = nrow(data))
    
    #for the structure Wi matrix is sufficient
    inverse_estimate <- inverse_estimate['wi']
    print(inverse_estimate)
    #Conversion to DataFrame, we keep the variablenames for visualisation
    inverse_df <- do.call(rbind.data.frame, inverse_estimate)
    colnames(inverse_df) <- colnames(data)
    rownames(inverse_df) <- colnames(data)
    
    #Determine graph structure - if Wi[i,j] < limit => No edge i-j
    graph_struct <- threshold_wi(inverse_df)
    
    #Visualising graph structure
    graph_struct <- data.matrix(graph_struct)
    network <- graph_from_adjacency_matrix(graph_struct,
                                           mode = "undirected")
    networks[[length(networks)+1]] <- network 
  }
  return(networks)
}

options(scipen = 999)
#Example of function calls 
lambdas <- list(0, 250, 750, 1250, 4200, 10000, 15000, 20000, 50000, 90000)
networks <- networks_GLasso(data = data, 
                            lambdas = lambdas)

par(mfrow=c(2,5))
i<-1
for (network in networks){
  plot(network , 
       edge.arrow.size=.0005, 
       directed = F, 
       layout=layout.circle,
       vertex.color = "White",
       vertex.size = 4,  
       vertex.label.color = "Black",
       vertex.label.dist = 2,
       vertex.label.cex = 1.6,
       vertex.label.font = 1,
       vertex.label.family = "sans",
       edge.color = "Blue",
       edge.width = 1.5,
       main = as.character(lambdas[i]))
  i <- i+1
}

#Ploting difference from true graph - network
differences_true_net <- sapply(networks,difference_graph,
                               first = true_graph, 
                               both = FALSE)
par(mfrow=c(2,5))
i<-1
for (network in differences_true_net){
  
  plot(network , 
       edge.arrow.size=.0005, 
       directed = F, 
       layout=layout.circle,
       vertex.color = "White",
       vertex.size = 4,  
       vertex.label.color = "Black",
       vertex.label.dist = 2,
       vertex.label.cex = 1.6,
       vertex.label.font = 1,
       vertex.label.family = "sans",
       edge.color = "red",
       edge.width = 1.5,
       main = as.character(lambdas[i]))
  i <- i+1
}

#Plotting diference network - true graph
differences_net_true <- sapply(networks,difference_graph,
                               second = true_graph, 
                               both = FALSE)
par(mfrow=c(2,5))
i<-1
for (network in differences_net_true){
  
  plot(network , 
       edge.arrow.size=.0005, 
       directed = F, 
       layout=layout.circle,
       vertex.color = "White",
       vertex.size = 4,  
       vertex.label.color = "Black",
       vertex.label.dist = 2,
       vertex.label.cex = 1.6,
       vertex.label.font = 1,
       vertex.label.family = "sans",
       edge.color = "Green",
       edge.width = 1.5,
       main = as.character(lambdas[i]))
  i <- i+1
}

#Finding closest graphs for the PC algorithm and GLasso to the true graph
#11 vertices (proteins)
gorder(true_graph)
#19 edges
gsize(true_graph)

#Closest in density
sizes_lasso <- sapply(networks, gsize)
index_density <-which.min(abs(sizes_lasso - gsize(true_graph)))

#Closest by size of Edge difference
differences_both_sides <- sapply(networks,difference_graph,
                               second = true_graph, 
                               both = TRUE)
sizes_diff <- sapply(differences_both_sides, gsize)
#Because the structure is sparse the "best" graph is empty
sizes_diff

par(mfrow=c(1,3))
#True graph
plot(true_graph , 
     edge.arrow.size=.0005, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White",
     vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 4.3,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 3,
     main = "True Graph")
#Closest GLasso graph
plot(networks[[index_density]] , 
     edge.arrow.size=.0005, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White",
     vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 4.3,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 3,
     main = "GLasso - 4200")
#Closest PC algo graph
skeletons <- pc_skeletons(data, c(0.005,0.0005, 0.0001, 0.00005))
sizes_pc <- sapply(skeletons,gsize)
index_density_pc <-which.min(abs(sizes_pc - gsize(true_graph)))
plot(skeletons[[index_density_pc]] , 
     edge.arrow.size=.0005, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White",
     vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 4.3,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 3,
     main = "PC - 0.0005")
gsize(true_graph)
gsize(networks[[index_density]])
gsize(skeletons[[index_density_pc]])

#Intersects with true graph
par(mfrow=c(1,2))
#Closest GLasso graph
plot(networks[[index_density]] %s% true_graph, 
     edge.arrow.size=.0005, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White",
     vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 4.3,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 3,
     main = "GLasso - 4200")
#Closest PC algo graph
skeletons <- pc_skeletons(data, c(0.005,0.0005, 0.0001, 0.00005))
sizes_pc <- sapply(skeletons,gsize)
index_density_pc <-which.min(abs(sizes_pc - gsize(true_graph)))
plot(skeletons[[index_density_pc]] %s% true_graph, 
     edge.arrow.size=.0005, 
     directed = F, 
     layout=layout.circle,
     vertex.color = "White",
     vertex.size = 4,  
     vertex.label.color = "Black",
     vertex.label.dist = 2,
     vertex.label.cex = 4.3,
     vertex.label.font = 1,
     vertex.label.family = "sans",
     edge.color = "Blue",
     edge.width = 3,
     main = "PC - 0.0005")

