library(ggplot2)
library(tidyverse)
library(pheatmap)

# Import RNA-seq data (returns a data frame)
cov <- read.csv("covid19.csv") %>%
        select(-X1)


# Explore your data frame 
print(dim(cov))   # number of rows and columns
print(cov[1:10, 1:10])  # first 10 rows and 10 columns


# Create feature_matrix (input for running dimensionality reduction)
feature_matrix <- cov[, 4:ncol(cov)]
print(feature_matrix[1:5, 1:5])

########################################## PCA ########################################## 


# Run PCA (returns a PCA object)
pca <- prcomp(feature_matrix, 
              center = TRUE,    # Expression counts are zero-centered
              scale = TRUE)     # Expression counts are brought to the same scale 
                                # to avoid dominance by highly expressed genes 
                   

# Check your PCA result 
# Note that when your total variance is set to 1, 
# count how many PCs (dimensions) are needed to explain and cluster your data 
# effectively 
summary(pca)

# Determine minimum number of PCs (dimensions) for your data 

# Extract proportion of variance explained by PC1-10 
# and store as a data frame pca_scree
pve <- pca$sdev^2 / sum(pca$sdev^2) 
cum_pve <- cumsum(pve)
pca_df <- data.frame(prop_var = pve[1:10], 
                     cum_prop_var = cum_pve[1:10], 
                     PC = factor(1:10))

print(pca_df)

# Create a scree (elbow) plot about proportion of variance change 
prop_var_plot <- ggplot(pca_df,
                        aes(x = PC, 
                            y = prop_var,
                            group = 1)) +
        geom_line() + 
        geom_point() +
        labs(title = "Scree (elbow) plot: Proportion of Variance Explained by PC1-10",
             y = "Proportion of Variance Explained")

print(prop_var_plot)



# Create a plot about cumulative proportion of variance change 
cum_prop_var_plot <- ggplot(pca_df,
                            aes(x = PC, 
                                y = cum_prop_var,
                                group = 1)) +
        geom_line() + 
        geom_point() +
        labs(title = "Cumulative Proportion of Variance Explained by PC1-10",
             y = "Cumulative Proportion of Variance Explained")

print(cum_prop_var_plot)


# Extract PCA coordinates (PC1-PC3) and clean data 
# (returns a data frame)
pca_coord <- data.frame(Sample = cov$Sample,
                        X = pca$x[, 1],
                        Y = pca$x[, 2],
                        Z = pca$x[, 3]) %>% 
        inner_join(cov[, 1:3], by = "Sample") %>%
        mutate(Covid19 = ifelse(str_detect(Sample, "N"), "Negative", "Positive"))

print(pca_coord)


# Plot the PCA result in 2D space (total observations)
pca_2D_plot1 <- ggplot(pca_coord,
                      aes(x = X, 
                          y = Y, 
                          color = ICU,
                          shape = gender)) + 
        geom_point(size = 2, alpha = 0.5) +
        labs(title = "PCA",
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(pca_2D_plot1)
# left: nonICU, right: ICU

# Plot the PCA result in 2D space (splitted observations by infected status)
pca_2D_plot2 <- pca_2D_plot1 + 
        facet_grid(~ Covid19)

print(pca_2D_plot2)

########################################## Hierarchical Clustering ##########################################


# Create a coordinate matrix 
coord_matrix <- cov[, 1:3] %>%
        unite(sample_code, gender, ICU, Sample) %>%
        cbind(pca_coord[, 2:4]) %>%
        column_to_rownames(var = "sample_code") %>%
        as.matrix()

print(coord_matrix)


# Calculate distance: (X, Y, Z) coordinates
# (returns a distance matrix)
distance <- dist(coord_matrix, 
                 method = "euclidean")

print(distance)


# Perform hierarchical clustering 
Hierarchical_clustering <- hclust(distance, 
                                  method = "average")

# Create a dendrogram
plot(Hierarchical_clustering)

# Heatmap + dendrogram
meta <- pca_coord[, 5:7] 
rownames(meta) <- rownames(coord_matrix)

pheatmap(coord_matrix,
         clustering_distance_rows = "euclidean",
         annotation_row = meta,
         fontsize_row = 5,
         main = "Heatmap")


# Cluster by cutting the tree 
# (returns a series of cluster number)
cut_by_k <- cutree(Hierarchical_clustering, k = 8)
cut_by_h <- cutree(Hierarchical_clustering, h = 100)

print(cut_by_k)
print(cut_by_h)

# Combine the clustering result with the (x, y, z) coordinates
# (returns a data frame)
pca_coord$hcluster_k <- factor(cut_by_k)
pca_coord$hcluster_h <- factor(cut_by_h)

print(pca_coord)
        
# Clean the data frame before plotting
pca_coord_cleaned <- gather(pca_coord,
                            hclustering_by, 
                            hcluster_number, 
                            c("hcluster_k", "hcluster_h"))

print(pca_coord_cleaned)

# Plot hierarchical clustering results 
hierarchical_plot1 <- ggplot(pca_coord_cleaned,
       aes(x = X,
           y = Y, 
           color = hcluster_number)) +
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(~ hclustering_by) + 
        labs(title = "Hierarchical clustering (k = 8, h = 100)",
             x = "PC1 (28%)",
             y = "PC2 (19%)")

plot(hierarchical_plot1)


hierarchical_plot2 <- ggplot(pca_coord_cleaned,
                             aes(x = X, 
                                 y = Y, 
                                 color = ICU, 
                                 shape = Covid19)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(hcluster_number ~ hclustering_by) + 
        labs(title = "Hierarchical clustering (k = 8, h = 100)",
             x = "PC1 (28%)",
             y = "PC2 (19%)")


plot(hierarchical_plot2)

# hcluster_k 
# cluster1: positive & nonICU 
# cluster2: positive 
# cluster4: ICU 
# cluster 3 & 7: (relatively) negative 






########################################## K-means Clustering ##########################################

# Determine optimal k (= number of clusters) by calculating total within cluster 
# sum of squares set.seed(32)
ttWithinss <- map_dbl(1:10, 
                      function(k) {
                              set.seed(32)
                              km <- kmeans(x = coord_matrix, 
                                           centers = k,
                                           nstart = 25)
                              km$tot.withinss})

print(ttWithinss)

# Create a scree plot 
ttWithinss_plot <- data.frame(tt_within_ss = ttWithinss,
                              k = factor(1:10)) %>%
        ggplot(aes(x = k, 
                   y = tt_within_ss,
                   group = 1)) +
        geom_point() + 
        geom_line() + 
        geom_vline(xintercept = 3, color = "red") + 
        labs(title = "Total Within-Cluster Sum of Squares Change",
             y = "Total Within-Cluster Sum of Squares",
             x = "Number of Clusters (k)")

print(ttWithinss_plot)

# Run k-means clustering with k = 3 
# (returns a kmeans object)
set.seed(32)
kmc <- kmeans(x = coord_matrix, 
              centers = 3,
              nstart = 25)


# Extract cluster results from the kmeans object
# (returns a series of cluster numbers)
kmc_cluster <- factor(kmc$cluster)

print(kmc_cluster)

# Combine the kmeans clustering result with the pca coordinate table 
pca_coord$kmcluster <- kmc_cluster

print(pca_coord)

# Plot the result of k-means clustering
kmeans_plot1 <- ggplot(pca_coord,
                       aes(x = X,
                           y = Y, 
                           color = kmcluster)) + 
        geom_point(size = 2, alpha = 0.5) + 
        labs(title = "K-Means Clustering", 
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(kmeans_plot1)


# Visualize characteristics of each cluster
kmeans_plot2 <- ggplot(pca_coord,
                       aes(x = X,
                           y = Y, 
                           color = ICU, 
                           shape = Covid19)) + 
        geom_point(size = 2, alpha = 0.5) +  
        facet_grid(gender ~ kmcluster) + 
        labs(title = "K-Means Clustering", 
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(kmeans_plot2)




########################################## t-SNE ##########################################


# Load the Rtsne package 
library(Rtsne)


# Set a function for running t-SNE 
tSNE_fn <- function(pp) {
        
        # tSNE
        set.seed(277)
        Rtsne(as.matrix(feature_matrix), 
              PCA = TRUE,              # Preliminary PCA
              perplexity = pp,         # Perplexity
              max_iter = 2000,         # Max iteration number  
              dims = 2)                # Number of output dimensions 
}

# Run t-SNE (returns t-SNE objects)
tsne1 <- tSNE_fn(1)
tsne2 <- tSNE_fn(2)
tsne5 <- tSNE_fn(5)
tsne10 <- tSNE_fn(10)

# Set a function for cleaning data 
data_clean_fn <- function(tsne_object, perplexity) { 
        
        data.frame(X = tsne_object$Y[, 1],
                   Y = tsne_object$Y[, 2],
                   gender = pca_coord$gender, 
                   ICU = pca_coord$ICU,
                   Covid19 = pca_coord$Covid19,
                   Perplexity = factor(perplexity)) 
}

# Clean your data
# (returns a data frame)
tsne_compare_df <- rbind(data_clean_fn(tsne1, 1),
                         data_clean_fn(tsne2, 2),
                         data_clean_fn(tsne5, 5),
                         data_clean_fn(tsne10, 10))

print(tsne_compare_df)

# Check out the relationship btw perplexity and the coordinates 
perplexity_plot <- ggplot(tsne_compare_df,
                          aes(x = X, 
                              y = Y, 
                              color = ICU, 
                              shape = Covid19)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(~ Perplexity) + 
        labs(title = "t-SNE with Perplexity = 1, 2, 5, and 10", 
             x = "Dim1",
             y = "Dim2")

print(perplexity_plot)

# Subset with Perplexity = 2
tsne_pp <- filter(tsne_compare_df, Perplexity == 2)
rownames(tsne_pp) <- rownames(meta)


# Calculate distance: (X, Y, Z) coordinates
# (returns a distance matrix)
distance_tsne <- dist(tsne_pp[, 1:2], 
                 method = "euclidean")

# Perform hierarchical clustering 
Hierarchical_clustering_tsne <- hclust(distance, 
                                  method = "average")

# Create a dendrogram
plot(Hierarchical_clustering_tsne)

# Extract the clustering result (8 clusters)

hcluster_tsne <- cutree(Hierarchical_clustering_tsne,
                        k = 8)

print(hcluster_tsne)

# Clean data
tsne_pp$hcluster <- factor(hcluster_tsne)

print(tsne_pp)

# Plotting tSNE/hierarchical clustering results 
tsne_hclustering1 <- ggplot(tsne_pp,
                           aes(x = X,
                               y = Y, 
                               color = hcluster)) + 
        geom_point(size = 2, alpha = 0.5) + 
        labs(title = "t-SNE and Hierarchical Clustering",
             x = "Dim1",
             y = "Dim2")

print(tsne_hclustering1)

tsne_hclustering2 <- ggplot(tsne_pp,
                            aes(x = X,
                                y = Y, 
                                color = ICU,
                                shape = Covid19)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid( ~ hcluster) + 
        labs(title = "t-SNE and Hierarchical Clustering",
             x = "Dim1",
             y = "Dim2")

print(tsne_hclustering2)

# Cluster1: positive & nonICU 
# Cluster2: positive & ICU 
# Cluster4: ICU 
# Cluster7: nonICU



########################################## Dimensionality Reduction in a scRNA-seq pipleline ##########################################

# Load Seurat package
library(Seurat)

# Data cleaning: row = gene, column = sample
# (returns a count matrix)
rownames(feature_matrix) <- rownames(coord_matrix)
count_matrix <- t(feature_matrix)

# Explore your input count matrix
print(dim(count_matrix))
print(count_matrix[1:10, 1:10])

# Create a seurat object 
seurat <- CreateSeuratObject(counts = count_matrix,
                             meta.data = meta)

# Normalize your data 
seurat <- NormalizeData(seurat,
                        normalization.method = "LogNormalize")
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)


# Run PCA 
seurat <- RunPCA(seurat)

# Cluster 
seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat)

# Run t-SNE (with perplexity 2)
seurat <- RunTSNE(seurat, 
                  perplexity = 2)


# Plot the PCA result 
DimPlot(object = seurat, 
        reduction = "pca") + 
        ggtitle("PCA") 

# Plot the t-SNE result
DimPlot(object = seurat, 
        reduction = "tsne") + 
        ggtitle("t-SNE") 

# Calculate the Number of Significant Features 
VariableFeaturePlot(seurat)

# Plote distribution of eigenvalues of the first dimension in each cluster
VlnPlot(seurat, features = "PC_1")
VlnPlot(seurat, features = "PC_1", split.by = "ICU")
VlnPlot(seurat, features = "PC_1", split.by = "Covid19")
VlnPlot(seurat, features = "PC_1", split.by = "gender")
VlnPlot(seurat, features = "tSNE_1")

# Heatmap (Gene expression profile in each cluster)
DoHeatmap(seurat)


