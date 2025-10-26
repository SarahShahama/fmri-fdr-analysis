install.packages(c("fdrtool", "qvalue", "knockoff", "oro.nifti", "ggplot2", "dplyr", "tidyr"))
library(oro.nifti)
library(ggplot2)
# Set the working directory to where the data is stored
setwd("C:/Users/HP/OneDrive/dataset/sub-1/func")
# Load the NIfTI file
fmri_img <- readNIfTI("sub-1_task-objectviewing_run-01_bold.nii.gz", reorient = FALSE)
# Check basic information
print(fmri_img)
# Display the 20th slice at timepoint 1
image(fmri_img[,,20,1], col = gray(0:64/64), main = "Slice 20 - Time 1")
library(readr)
# Read the event (stimulus) file
events <- read_tsv("sub-1_task-objectviewing_run-01_events.tsv")
head(events)
n_vols <- dim(fmri_img)[4]  # number of timepoints
n_vols
#there are 121 volumes
# TR = 2.5 seconds
TR <- 2.5
# Convert onset times to corresponding volume indices
events$volume_index <- round(events$onset / TR) + 1
head(events$volume_index)#these are index of the fmri volume corresponding to each trails
#extracting only these volumes that correspond to actual trials
# Make sure indices don’t exceed total volumes
events <- events[events$volume_index <= dim(fmri_img)[4], ]
# Subset the fMRI image to only those volumes
fmri_task <- fmri_img[,,,events$volume_index]
# Check new dimension (should be 96 volumes)
dim(fmri_task)#now we have only task related data
events_scissors_chair <- events[events$trial_type %in% c("scissors", "chair"), ]
table(events_scissors_chair$trial_type)
# Make sure the number of volumes matches the events
dim(fmri_task)[4] == nrow(events)
# Indices for each condition
idx_chair <- which(events$trial_type == "chair")
idx_scissors  <- which(events$trial_type == "scissors")
length(idx_chair); length(idx_scissors) #how many scans belong to each condition
# Convert the 4D array to a 2D matrix
n_voxels <- prod(dim(fmri_task)[1:3])  # total number of voxels
n_timepoints <- dim(fmri_task)[4]
fmri_matrix <- matrix(fmri_task, nrow = n_voxels, ncol = n_timepoints)
dim(fmri_matrix)
# Initialize vector to store p-values
p_values <- numeric(n_voxels)
# Run voxel-wise t-tests
for (v in seq_len(n_voxels)) {
  voxel_data <- fmri_matrix[v, ]
  group1 <- voxel_data[idx_chair]
  group2 <- voxel_data[idx_scissors]
  
  test_result <- try(t.test(group1, group2), silent = TRUE)
  
  if (inherits(test_result, "htest")) {
    p_values[v] <- test_result$p.value
  } else {
    p_values[v] <- NA
  }
}
# Check summary
summary(p_values)
hist(p_values, breaks = 50, col = "skyblue", main = "Voxel-wise p-value distribution",xlab = "p-value")
# Apply BH correction
p_adj_bh <- p.adjust(p_values, method = "BH")
# Threshold for significance (FDR < 0.05)
alpha <- 0.05
significant_voxels_bh <- which(p_adj_bh < alpha)
# Summary
length(significant_voxels_bh)#This give us how many voxels survived the BH correction
#Visualize the results
# Create an empty 3D mask of zeros
sig_map_bh <- array(0, dim = dim(fmri_task)[1:3])
# Mark significant voxels with 1
# Visualize one slice
slice_index <- 20  
image(sig_map_bh[,,slice_index], col = gray(0:1), main = paste("Significant voxels (BH) - Slice", slice_index))
#Storey’s q-value Method
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
library(qvalue)
# Compute q-values from p-values
qobj <- qvalue(p = p_values)

# Extract q-values
q_values <- qobj$qvalues

# Threshold for significance
alpha <- 0.05
significant_voxels_q <- which(q_values < alpha)

# Summary
length(significant_voxels_q)
#voxels that are significant under Storey’s q-value control is 106796 while under BH was only 88374
# Create a binary map of significant voxels
sig_map_q <- array(0, dim = dim(fmri_task)[1:3])
sig_map_q[significant_voxels_q] <- 1

# Plot one slice
slice_index <- 20
image(sig_map_q[,,slice_index], col = gray(0:1),
      main = paste("Significant voxels (Storey's q-value) - Slice", slice_index))
cat("BH significant voxels:", length(significant_voxels_bh), "\n")
cat("Storey q-value significant voxels:", length(significant_voxels_q), "\n")
#seeing where they overlap
overlap <- intersect(significant_voxels_bh, significant_voxels_q)
cat("Overlap:", length(overlap), "voxels\n")
#visualizing where q-value finds extra voxels not found by BH
extra_q <- setdiff(significant_voxels_q, significant_voxels_bh)
extra_map <- array(0, dim = dim(fmri_task)[1:3])
extra_map[extra_q] <- 1

image(extra_map[,,slice_index], col = gray(0:1),
      main = paste("Extra voxels detected by q-value (Slice", slice_index, ")"))
install.packages("knockoff")
library(knockoff)
#For simplicity, we’ll use a smaller random subset of voxels to make this computationally feasible
#We’ll also create a simple design matrix representing the two conditions (e.g., chair vs scissors)
unique(events$trial_type)
#conditions: chair vs scissors 
cond1 <- "chair"
cond2 <- "scissors"
data_1<-fmri_task[,,,idx_chair]
data_2<-fmri_task[,,,idx_scissors]
# Confirm dimensions
cat("Cond1:", dim(data_1), "\nCond2:", dim(data_2))
# Reshape to 2D (timepoints × voxels)
data_cond1 <- t(apply(data_1, 4, function(x) as.vector(x)))
data_cond2 <- t(apply(data_2, 4, function(x) as.vector(x)))
dim(data_cond1)
dim(data_cond2)
# Combine both conditions
X <- rbind(data_cond1, data_cond2)
# Create response labels
y <- c(rep(1, nrow(data_cond1)), rep(0, nrow(data_cond2)))
# Reduce to a manageable subset (e.g., 1000 voxels)
set.seed(123)
subset_voxels <- sample(ncol(X), 1000)
X_sub <- X[, subset_voxels]
dim(X_sub)
# Remove voxels (columns) that have zero variance or NA values
good_voxels <- apply(X_sub, 2, function(x) all(is.finite(x)) && sd(x) > 0)
X_clean <- X_sub[, good_voxels]
# Check dimensions after cleaning
dim(X_clean)
#standardizing each column
x_scaled<-scale(X_clean,center = TRUE,scale = TRUE)
#computing covariance on standardized data
sigma<-cov(x_scaled)
# Small regularization: add a tiny diagonal term for stability
Sigma <- sigma + diag(1e-6, ncol(sigma))
# Use Lasso-based statistic for feature selection
result <- knockoff.filter(
  X = x_scaled,
  y = y,
  knockoffs = function(X) create.gaussian(X, mu = rep(0, ncol(X)), Sigma = Sigma),
  statistic = stat.lasso_lambdasmax
)
selected_features<-result$selected
length(selected_features)
cat("B-H DISCOVERIES:",length(significant_voxels_bh))
cat("Q VALUE DISCOVERIES:",length(significant_voxels_q))
cat("KNOCKOFF DISCOVERIES:",length(selected_features))
sig_map_bh <- array(0, dim = dim(fmri_task)[1:3])
sig_map_q  <- array(0, dim = dim(fmri_task)[1:3])
# Assign 1s for significant voxels
sig_map_bh[significant_voxels_bh] <- 1
sig_map_q[significant_voxels_q] <- 1
# Select a slice index
slice_index <- 20
# Separate visualizations
par(mfrow = c(1, 2))
image(sig_map_bh[, , slice_index],
      col = gray(0:1),
      main = paste("Significant Voxels (BH) - Slice", slice_index))
image(sig_map_q[, , slice_index],
      col = gray(0:1),
      main = paste("Significant Voxels (Storey q-value) - Slice", slice_index))
# Combine maps into one
combined_map <- array(0, dim = dim(fmri_task)[1:3])
combined_map[significant_voxels_bh] <- 1          # BH significant
combined_map[significant_voxels_q]  <- combined_map[significant_voxels_q] + 2  # q-value significant
# Define colors:
# 0 = white (none)
# 1 = red (BH only)
# 2 = blue (q-value only)
# 3 = purple (both)
overlay_colors <- c("white", "red", "blue", "purple")
# Visualize overlay
image(combined_map[, , slice_index],
      col = overlay_colors,
      main = paste("Overlap: BH vs Storey q-value - Slice", slice_index))
#comparing with knockoff filter
#Create Knockoff map
sig_map_knock <- array(0, dim = dim(fmri_task)[1:3])
sig_map_knock[selected_features] <- 1
# Visualize alongside
par(mfrow = c(1, 2))
image(sig_map_bh[, , slice_index],
      col = gray(0:1),
      main = paste("BH FDR - Slice", slice_index))
image(sig_map_q[, , slice_index],
      col = gray(0:1),
      main = paste("Storey q-value - Slice", slice_index))