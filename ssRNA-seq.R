
# Contents:
# Single cells RNA-seq analysis with K-mean clustering (uses all functions created in biomarkersignatures.R Parts 1-4)


# setwd
setwd("~/data")

# Part 5: Single Cell Omics 
#new em with single cell samples
em1 = read.table("EM.csv", header = TRUE, row.names = 1, sep="\t")

# PCA function: accepts a vector of groups, returns a PCA plot for PC1 vs PC2 
make_pc1_pc2 = function(colour_groups, e_data){
  # scale data
  em_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  # run PCA
  xx = prcomp(t(em_scaled)) 
  pca_coordinates = data.frame(xx$x)
  # get % variation
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100 
  prop_y = round(vars["PC2"] / sum(vars),4) * 100 
  x_axis_label = paste("PC1 (", prop_x, "%)", sep="") 
  y_axis_label = paste("PC2 (", prop_y, "%)", sep="")
  # plot
  ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = colour_groups)) +
    geom_point() +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) + 
    theme_bw()
  
  return(ggp) 
  }


# PCA function: accepts a vector of groups, returns a PCA plot for PC3 and PC4
make_pc3_pc4 = function(colour_groups, e_data){
  # scale data
  em_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  # run PCA
  xx = prcomp(t(em_scaled)) 
  pca_coordinates = data.frame(xx$x)
  # get % variation
  vars = apply(xx$x, 2, var)
  print(vars)
  prop_x = round(vars["PC3"] / sum(vars),4) * 100 
  prop_y = round(vars["PC4"] / sum(vars),4) * 100 
  x_axis_label = paste("PC3 (" ,prop_x, "%)", sep="") 
  y_axis_label = paste("PC4 (" ,prop_y, "%)", sep="")
  # plot
  ggp = ggplot(pca_coordinates, aes(x=PC3, y= PC4, colour = colour_groups)) +
    geom_point() +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) + theme_bw()
  
  return(ggp) 
}

#make PCA for comp 1&2 and 3&4
make_pc1_pc2("", em1)
make_pc3_pc4("", em1)

get_component_genes = function(component, e_data) {
    # scale data
    e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
    # run PCA
    xx = prcomp(t(e_data_scaled)) 
    pca_coordinates = data.frame(xx$x)
    # get the samples in the upper and lower quartile
    summary(pca_coordinates[,component])
    q1 = summary(pca_coordinates[,component])[2] 
    q3 = summary(pca_coordinates[,component])[5]
    q1_samples = row.names(subset(pca_coordinates, pca_coordinates[,component] < q1)) 
    q3_samples = row.names(subset(pca_coordinates, pca_coordinates[,component] > q3))
    # get the p and sig genes
    de = get_de(q1_samples, q3_samples, e_data)
    markers = row.names(subset(de, p.adj < 0.01 & abs(log2fold) > 2))
    
    results = list("markers" = markers, "q1_samples" = q1_samples, "q3_samples" = q3_samples) 
    return(results)
  }

#get marker genes for PC1, PC2 and PC3
markers_PC1_em1 = get_component_genes(1, em1)  #### NOT WORKING
markers_PC2_em1 = get_component_genes(2, em1)  
markers_PC3_em1 = get_component_genes(3, em1)

# Get component genes (manually)
# scale data
em1_scaled = na.omit(data.frame(t(scale(t(em1)))))

# run PCA
xx = prcomp(t(em1_scaled)) 
pca_coordinates = data.frame(xx$x)
# get the samples in the upper and lower quartile
summary(pca_coordinates[,1])
q1_PC1_em1 = summary(pca_coordinates[,1])[2] 
q3_PC1_em1 = summary(pca_coordinates[,1])[5]
q1_samples_em1 = row.names(subset(pca_coordinates, pca_coordinates[,1] < q1)) 
q3_samples_em1 = row.names(subset(pca_coordinates, pca_coordinates[,1] > q3))
# get the p and sig genes
de_em1 = get_de(q1_samples_em1, q3_samples_em1, em1)
markers_PC1_em1 = row.names(subset(de_em1, p.adj < 0.05 & abs(log2fold) > 1))

# create vector containing only samples in q1 and q3 to order em
ordered_em1 = c(q1_samples_em1, q3_samples_em1)
#order em
em1_sliced = em1[,ordered_em1] 

# make heatmap using ordered dataset and markers list (PC1)
make_heatmap(em1_sliced, markers_PC1_em1)

# K-means clustering
# get_de function 
get_de = function(group_1,group_2, e_data){
  # slice EM table using the vectors just created 
  group_1_em = e_data[,group_1]
  group_2_em = e_data[,group_2]
  
  #create new DE
  de = data.frame(matrix(nrow=nrow(e_data), ncol=3))
  row.names(de) = row.names(e_data)
  names(de) = c("p", "p.adj", "log2fold")
  
  for (row_number in 1:nrow(e_data)) {
    # gets the expression data for each group, for the current gene
    group_1_candidate = as.numeric(group_1_em[row_number,]) 
    group_2_candidate = as.numeric(group_2_em[row_number,])
    # calculates P
    w_result = wilcox.test(group_1_candidate, group_2_candidate, alternative = "two.sided")
    # stores P
    de[row_number,"p"] = w_result$p.value
    #calculates p.adj
    de[row_number,"p.adj"] = p.adjust(de[row_number,"p"], method = p.adjust.methods, n=3)
    #calculates log2fold change
    log2fold = log(mean(group_1_candidate),2)-log(mean(group_2_candidate),2)
    #stores log2folds
    de[row_number,"log2fold"] = log2fold
  }
  return(de)
}

# 1) Scale the data
em1.s = na.omit(data.frame(t(scale(t(em1)))))
# 2) transpose and perform k-means
km = kmeans(t(em1.s), 3, iter.max = 10, nstart = 1) # 2 is the no. of clusters to find, 10 = no. iterations
# 3) extract cluster information 
clusters = as.factor(km$cluster)

# use cluster info to colour PCA
#make PCA for comp 1&2 and 3&4
make_pc1_pc2(clusters, em1)
make_pc3_pc4(clusters,em1)

# discover the optimal number of clusters
install.packages("factoextra")
fviz_nbclust(t(em.s), kmeans, method = "wss")

# clusters as a dataframe
cluster_data = data.frame(km$cluster)
View(cluster_data)
# rename column
names(cluster_data) = "km.cluster"

# get the genes in each cluster 
cluster_1 = row.names(subset(cluster_data, km.cluster == 1))
cluster_2 = row.names(subset(cluster_data, km.cluster == 2))
cluster_3 = row.names(subset(cluster_data, km.cluster == 3))

# find the markers for these groups using the get_de function
de1v2 = get_de(cluster_1, cluster_2, em1)
de1v3 = get_de(cluster_1, cluster_3, em1)
de2v3 = get_de(cluster_2, cluster_3, em1)
View(de1v2)

# get sig genes from each de table
de1v2_sig = subset(de1v2, p.adj < 0.05 & abs(log2fold) > 1)
de1v3_sig = subset(de1v3, p.adj < 0.05 & abs(log2fold) > 1)
de2v3_sig = subset(de2v3, p.adj < 0.05 & abs(log2fold) > 1)

# sort the sig de tables 
de1v2_sig = de1v2_sig[order(de1v2_sig$p),] #order orders argument by ascending order
de1v3_sig = de1v3_sig[order(de1v3_sig$p),]
de2v3_sig = de2v3_sig[order(de2v3_sig$p),]

# take top 25 genes from each list of sig genes
de1v2_top25 = row.names(de1v2_sig[1:25,])
de1v3_top25 = row.names(de1v3_sig[1:25,])
de2v3_top25 = row.names(de2v3_sig[1:25,])

# create list of sig genes (unique IDs)
sig_genes = unique(c(de1v2_top25, de1v3_top25, de2v3_top25))

# create vector for ordering clusters 
ordered_clusters = c(cluster_1, cluster_2, cluster_3)

# heatmap function
make_heatmap = function(em_table, gene_list) {
  # prepares the table
  hm.matrix = as.matrix(em_table[gene_list,])
  # does the clustering
  y.dist = Dist(hm.matrix, method="spearman") 
  y.cluster = hclust(y.dist, method="average") 
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average") 
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  # melts
  hm.matrix_clustered.m = melt(hm.matrix_clustered)
  # chooses colours
  palette = colorRampPalette(c("blue","pink","red"))(100)
  # plots
  ggp = ggplot(hm.matrix_clustered.m, aes(x=Var2, y=Var1, fill=value)) + geom_tile() +
    scale_fill_gradientn(colours = palette) +
    theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0)) + labs(x="",y="")
  # returns the plot
  return(ggp) 
}

# make ordered heatmap
make_heatmap(em1[,c(cluster_1, cluster_2, cluster_3)], sig_genes)

