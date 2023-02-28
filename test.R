#Testing pcalg and visualisation package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("RBGL", version = "3.14")
BiocManager::install("graph", version = "3.14")
BiocManager::install("Rgraphviz", version = "3.14")


#load the pcalg package
library(pcalg)


#load the data
data <- read.csv("./data/Mid-Atlantic_Wage_Data_974_40.csv")
# encode education and health
data
data$education <- as.integer(substr(data$education, 1, 1))
data$health <- as.integer(substr(data$health, 1, 1))
# Drop IDs
drops <- c("X.1","X","logwage")
data=data[ , !(names(data) %in% drops)]

# keep only the numeric or int variables 
data <- data[, sapply(data, is.numeric)]

sufficient_statistics<- list(C = cor(data), n = nrow(data))
column_names <- colnames(data)
skel_survey <- pcalg::skeleton(sufficient_statistics,
                         indepTest= gaussCItest,
                         labels=column_names,
                         alpha=0.0005)

# visualize the skeleton
plot(skel_survey, main="Skeleton of the survey data")

