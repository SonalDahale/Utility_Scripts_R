# Set the working directory to where your data files are stored
setwd("C:/Users/sd01165/Downloads/milk")

# Load count data and experimental metadata
cts <- as.matrix(read.csv("subset_milk_cts.txt", sep = "\t", row.names = 1))
coldata <- read.csv("coldata.txt", row.names=1, sep="\t")

# Load the DESeq2 package
library(DESeq2)

###############################
# Simple Model: Condition Only
###############################

# The simple model assesses how the 'Condition' (e.g., 'Good' vs 'Poor') affects gene expression.
# This is a basic differential expression analysis model, where we compare gene expression between different conditions.

simple.model <- as.formula(~ Condition)

# Creating the model matrix for 'Condition', automatically setting the first factor alphabetically as the baseline (e.g., 'Good').
model.matrix(simple.model, data = coldata)

# Applying DESeq2 to create a DESeqDataSet object, specifying the count data, experimental design, and condition variable
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = simple.model)

# Filter out genes with low counts across all samples
keep <- rowSums(counts(dds)) > 5
dds.filt <- dds[keep,]

# Run the DESeq analysis
dds <- DESeq(dds.filt)

# Get results for the 'Condition' variable (e.g., comparing 'Good' vs 'Poor')
results.condition <- results(dds, alpha = 0.05)
results.condition

# Summarize results for significantly upregulated or downregulated genes
sum(results.condition$padj < 0.05)  # How many genes have significant p-values?
sum(is.na(results.condition$padj))  # How many genes have NA p-values?

# Upregulated genes (log2FoldChange > 0) with significant p-values
sum(results.condition$padj < 0.05 & results.condition$log2FoldChange > 0, na.rm = TRUE)

# Downregulated genes (log2FoldChange < 0) with significant p-values
sum(results.condition$padj < 0.05 & results.condition$log2FoldChange < 0, na.rm = TRUE)

####################################
# Additive Model: Speed + Condition
####################################

# The additive model assesses the independent contributions of two variables, 'Speed' and 'Condition', to gene expression.
# Here, we assume that the effects of 'Speed' (e.g., 'Slow' vs 'Fast') and 'Condition' (e.g., 'Good' vs 'Poor') are independent.
# This model is useful when you suspect both variables impact gene expression independently.

additive.model <- as.formula(~ Speed + Condition)

# Creating the model matrix to visualize how 'Speed' and 'Condition' are included in the design
model.matrix(additive.model, data = coldata)

# Applying DESeq2 with the additive model
dds.additive <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = additive.model)

# Filter out genes with low counts
keep <- rowSums(counts(dds.additive)) > 5
dds.additive.filt <- dds.additive[keep,]

# Running DESeq analysis for the additive model
dds.additive <- DESeq(dds.additive.filt)

# Getting results for the 'Speed' and 'Condition' factors
results.additive <- results(dds.additive, alpha = 0.05)
results.additive

# Summarizing results for significant genes
sum(results.additive$padj < 0.05)  # Significant genes
sum(results.additive$padj < 0.05 & results.additive$log2FoldChange > 0, na.rm = TRUE)  # Upregulated genes
sum(results.additive$padj < 0.05 & results.additive$log2FoldChange < 0, na.rm = TRUE)  # Downregulated genes

#######################################################################
# Interactive Model: Speed + Condition + Speed:Condition
#########################################################################

# The interactive model examines the interaction between 'Speed' and 'Condition'.
# This is useful when you suspect that the effect of one variable (e.g., 'Speed') may vary depending on the level of another variable (e.g., 'Condition').
# For example, the effect of 'Speed' on gene expression may be different in 'Good' vs 'Poor' conditions.

interactive.model <- as.formula(~ Speed + Condition + Speed:Condition)

# Creating the model matrix to visualize the interaction term between 'Speed' and 'Condition'
model.matrix(interactive.model, data = coldata)

# Applying DESeq2 with the interactive model
dds.interactive <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = interactive.model)

# Filter out genes with low counts
keep <- rowSums(counts(dds.interactive)) > 5
dds.interactive.filt <- dds.interactive[keep,]

# Running DESeq analysis for the interactive model
dds.interactive <- DESeq(dds.interactive.filt)

# Getting results for the 'Speed', 'Condition', and interaction terms
results.interactive <- results(dds.interactive, alpha = 0.05)
results.interactive

# Summarizing results for significant genes
sum(results.interactive$padj < 0.05)  # Significant genes
sum(results.interactive$padj < 0.05 & results.interactive$log2FoldChange > 0, na.rm = TRUE)  # Upregulated genes
sum(results.interactive$padj < 0.05 & results.interactive$log2FoldChange < 0, na.rm = TRUE)  # Downregulated genes

# You can also extract specific comparisons (e.g., comparing 'Good' vs 'Poor' or 'Slow' vs 'Fast')
resultsNames(dds.interactive)  # Check available comparisons

# Extract specific results
results.interaction.fast <- results(dds.interactive, name = "Condition_Good_vs_Poor", alpha = 0.05)
results.interaction.slow <- results(dds.interactive, contrast = list(c("Condition_Good_vs_Poor", "Speed_Slow_vs_Fast")), alpha = 0.05)

sum(results.interaction.fast$padj < 0.05, na.rm = TRUE)
sum(results.interaction.slow$padj < 0.05, na.rm = TRUE)

# Plotting PCA for visualization (optional)
vstcounts <- vst(dds.interactive, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Speed"))
