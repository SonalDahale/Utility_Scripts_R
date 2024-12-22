# Load required library
library("factoextra")

# Inspect the RNA data
head(rna_raw) # View the first few rows of the raw RNA data

# Step 1: Prepare data for PCA
# Remove the first column (assumed to be gene IDs) for PCA computation
forPCA_rna <- rna_raw[, 2:5]

# Step 2: Perform PCA
# Run Principal Component Analysis (PCA) on the selected data
# 'scale = T' standardizes the variables to have unit variance
rna.pca <- prcomp(forPCA_rna, scale = TRUE)

# Step 3: Visualize Eigenvalues (Scree Plot)
# Display the proportion of variance explained by each principal component
fviz_eig(rna.pca)

# Step 4: Visualize PCA Variables
# Create a PCA variable plot colored by contribution
fviz_pca_var(
  rna.pca,
  col.var = "contrib", # Color variables based on their contribution to the PCs
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Gradient color scheme
  repel = TRUE        # Avoid overlapping text labels
)
