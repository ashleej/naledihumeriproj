


## Landmark analysis


## Geometric morphometry


## Hypothesis testing
  

# We can see here  that `homo_sapiens_01_21-landmarks.pp` is the problem, so let's drop it, and take another look.     

# keep only specimens with 21 landmarks,
# and ensure landmarks are in order, 
# and convert to array
nested_df_of_pp_data_H_sapiens_21_lmk <- 
nested_df_of_pp_data_H_sapiens %>%
  mutate(n_row = map_int(data, nrow)) %>% 
  filter(n_row == 21) %>% 
  mutate(data = map(data, ~arrange(.x, landmark) %>% 
                           select(-landmark))) %>% 
  mutate(data_array = map(data, simplify2array))

L <-  nested_df_of_pp_data_H_sapiens_21_lmk$data_array
three_d_array <- array(unlist(L), 
                       dim = c(nrow(L[[1]]), 
                               ncol(L[[1]]), 
                               length(L)))
      
# GPA
library(geomorph)

# Obtain Procrustes coordinates
H_sapiens_gpa <-gpagen(three_d_array)
summary(H_sapiens_gpa)
plot(H_sapiens_gpa)

# Convert coordinates to two-dimensional matrix
coords2d <- two.d.array(H_sapiens_gpa$coords)

# Calculate consensus and flatten to single vectors
consensus <- apply(H_sapiens_gpa$coords, c(1,2), mean)
consensusvec <- apply(coords2d, 2, mean)

# Calculate Procrustes residuals (Procrustes coordinates - consensus)
resids <- t(t(coords2d)-consensusvec)

# Calculate covariance matrix
P <- cov(resids)

# Calculate eigenvector and eigenvalues with SVD
pca.stuff <- svd(P)
eigenvalues <- pca.stuff$d
eigenvectors <- pca.stuff$u

# Calculate PCA scores
scores <- resids%*%eigenvectors

plotTangentSpace(three_d_array)
plot(scores[,1:2],asp=1, pch=20,cex=2)




#--------------------------------------
# all!
nested_df_of_pp_data_21_lmk <- 
nested_df_of_pp_data %>%
  mutate(n_row = map_int(data, nrow)) %>% 
  filter(n_row == 21) %>% 
  mutate(data = map(data, ~arrange(.x, landmark) %>% 
                      select(-landmark))) %>% 
  mutate(data_array = map(data, simplify2array))

L <-  nested_df_of_pp_data_21_lmk$data_array
three_d_array <- array(unlist(L), 
                       dim = c(nrow(L[[1]]), 
                               ncol(L[[1]]), 
                               length(L)))

IDs <- nested_df_of_pp_data_21_lmk$specimen

# GPA
library(geomorph)

# Obtain Procrustes coordinates
all_gpa <-gpagen(three_d_array)
summary(all_gpa)
plot(all_gpa)

pca.lands <- plotTangentSpace(all_gpa$coords, label=TRUE)

plot(pca.lands$pc.scores[,1:2],pch=15,xlab="PC1",ylab="PC2")
text(pca.lands$pc.scores[,1:2],IDs,pos=4,cex=.5)

# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12035
# https://www.sciencedirect.com.sci-hub.tw/science/article/pii/S0047248415000512?via%3Dihub

library(scatterplot3d)
scatterplot3d(pca.lands$pc.scores[,1:3])

# The shape differences between group mean can be visualized graphically, by obtaining the average landmark coordinates for each group and the overall mean, and plotting the differences as thin‐plate spline transformation grids:
ref <- mshape(all_gpa$coords)

gp1.mn <- mshape(all_gpa$coords[,,1:20])

plotRefToTarget(ref,gp1.mn,mag=2) # need links?


physignal(plethspecies$phy,Y.gpa$coords,iter=99)

#------------------------------------------------

# with Momocs
library(Momocs)

all_data_ldk <- 
  Ldk(coo=three_d_array %>% a2l,
      fac= data.frame( ID =  IDs))

fgProcrustes(all_data_ldk) # nope





#--------------------------------------

# take a look with the bone
library(geomorph)
hum <- read.ply("data/Humerus.ply")
shade3d(hum, color = "gray80")
with(test_pp, plot3d(x, y, z, type = 's', col = "red", add = T))
with(test_pp ,text3d(x,y,z,row.names(test_pp), add = T, cex =2, adj = c(1,1)))

library(geomorph)
library(rgl)
library(tidyverse)

# put in the specimen ID and the filename
specimen_name <- "modern_human"
file_name <- "Humerus.ply"

# read in PLY file, convert STL to PLY (not binary) with MeshLab
hum <- read.ply(str_glue("data/{file_name}"))

# place points
number_of_landmarks <- 3
hum_points <- digit.fixed(hum, fixed = number_of_landmarks)
# right-click to place, check console and 'y'-enter to confirm

# output points to CSV file to use later, one file per model
write.csv(hum_points, str_glue("{specimen_name}.csv"))

-->