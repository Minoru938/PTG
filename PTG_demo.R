## ---------------------------------------------------------------------------------------------------------------------------------
## Author: Minoru Kusaba
## Affiliation: School of Multidisciplinary Sciences Department of Statistical Science SOKENDAI
## Contact: kusaba@ism.ac.jp
## File name: PTG_demo.R
## Task: Demonstration of PTG using the elemental data
## Last update: 2019/12/18
## ---------------------------------------------------------------------------------------------------------------------------------

# Set directory
dir_atom = "/Users/kusaba/Downloads/PTG-master/" ## Set directory containing "the_element_data.csv" file (the elemental data used in the paper)
dir_PTG = "/Users/kusaba/Downloads/PTG-master/" ## Set directory containing "PTG.R" and "PTG_length_parameter.cpp" files (functions for PTG)

# Require to load "PTG.R" 
source(file =paste(dir_PTG,"PTG.R",sep = ""))

# Load and scale the element data
element_data = read.csv(paste(dir_atom, "the_element_data.csv", sep = "") ,row.names = 1)
data = scale(element_data)

# Show heatmap of the data
d = data
pheatmap(
  d, cluster_rows = FALSE, cluster_cols = FALSE, 
  legend = FALSE, show_rownames = TRUE, fontsize_col = 10,main = "data")

# Generate color values
library(RColorBrewer) ##for generating color
mycol = brewer.pal(5,name = "Pastel1")
period = element_data[,"period"]
fc = factor(period)
col_p = mycol[fc] ##color values for period

# Create and show nodes (5 × 5 square nodes)
grid_sum = crt_grid(grid_str = "rctg",xDim = 5,yDim = 5)
plot(grid_sum$grid,xlab = "",ylab = "",pch = 19)

# Run PTG (time consuming:it took about 10 minutes in my MacBook Pro)
# Constrain the estimation by lambda:the larger the lambda, the stronger the constraints
# n_sep determines the number of nodes to increase in the nodes expansion process (when n_sep = 1, 5 × 5 nodes are incremented to 9 × 9)
# rate defines a prior distribution of β (gamma distribution)
t_start = proc.time()
set.seed(1122)
ptg = PTG(grid_sum = grid_sum,data = data,lambda = 3,rate = 500,n_sep = 1)
result = get_resultof_PTG(result = ptg,data = data)
calculation_time = proc.time()-t_start

# Show the PTG-created layout of the 54 elements on the 9 × 9 square lattice (each element is color-coded by period)
x = result$map
cex = 4
scatter2D(x = x[,1], y = x[,2], bty = "n",col = col_p
          , pch = 19, cex = cex, main = "",type = "p",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
text2D(x = x[,1], y = x[,2],labels = rownames(data),
       add = TRUE, colkey = FALSE, cex = cex/4,bty = "g",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
text2D(x = x[,1]-0.055, y = x[,2]+0.055,labels = as.character(seq(1:nrow(data))),
       add = TRUE, colkey = FALSE, cex = cex/6,bty = "g",xlab = "",ylab = "",xaxt = "n",yaxt = "n")

# Show the element data and mapped nodes in the feature space (the dimension is reduced to 3 by PCA for visualization)
# Each element is color coded by period and black points indicate mapped nodes
pca = prcomp(x = data, scale = T)
summary(pca)
x = pca$x[, 1:3]
Y_pca = result$Y %*% pca$rotation 
y = Y_pca[,1:3]
plot3d(x,col = col_p,size = 8)
text3d(x+0.1,texts = rownames(data))
points3d(y,size = 8)

# Create and show PTG property landscapes for the learned PTG square table
# n_sep controls the resolution of the landscapes, lambda should be equal to the one used in PTG
# A feature value of the landscape is set by the column number of the data (col_num = 11 corresponds to electron negativity)
krg = kriging_calculation(result = ptg,n_sep = 3,lambda = 3)
show_kriging(krg = krg,data = data,col_num = 11)
