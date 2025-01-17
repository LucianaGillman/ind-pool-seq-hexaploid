# Load required packages
library(MASS)


# Load your dissimilarity matrix 
diss_matrix<- read.csv("RogersDistances_para_analisis.csv", header=TRUE)
rownames(diss_matrix) <- diss_matrix[,1]

##calculate MDS
MDS3 <- cmdscale(diss_matrix[,-1], k = 3, eig = T, add = FALSE, x.ret = FALSE)

## explained variance
explained_variance <- (MDS3$eig / sum(abs(MDS3$eig))) * 100

write.csv(explained_variance, "explained_variance.csv")




