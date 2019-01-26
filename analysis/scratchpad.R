


## Landmark analysis


## Geometric morphometry


## Hypothesis testing
  

# We can see here  that `homo_sapiens_01_21-landmarks.pp` is the problem, so let's drop it, and take another look.     

nested_df_of_pp_data_H_sapiens





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