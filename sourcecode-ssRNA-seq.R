
# Contents:
# Part 1:
# Subset em by clinical features, e.g. factors or continuous variables
# Serform statistical tests (t-test & wilcox) singularly or in series using loops, for a few genes or all genes
# Generate de tables with p values, p adj and log2fold 

# Part 2
# get_p function: calculates the p-value and log2fold of a t-test or a Wilcox test for single or multiple genes between two predefined groups
# get_de function: calculates the p-value and log2fold of a t-test or a Wilcox test for a whole expression matrix between two predefined groups

# Part 3
# Identify biomarker signatures using PCA components and quartiles

# Part 4
# Plot biomarker signatures using heatmaps

# Part 5 
# Single cells RNA-seq analysis with K-mean clustering

# Part 1
# setwd
setwd("~/data")

# load data
em = read.table("EM.csv", header = TRUE, row.names = 1, sep = "\t")
ss = read.table("SS.csv", header = TRUE, sep = "\t")


# 1) get the names of the samples in each group
# start by dividing continuous data into two factor levels
neutro_count = cut(ss$NEUTROPHIL_COUNT, 2, labels = c("Low_Neutrophils","High_Neutrophils"))  #by indicating $SAMPLE, we specify what col of the df to subset
ss$NEUTROPHIL_LEVEL = neutro_count

group_1 = subset(ss, NEUTROPHIL_COUNT == "Low_Neutrophils")$SAMPLE
group_2 = subset(ss, NEUTROPHIL_COUNT == "High_Neutrophils")$SAMPLE

# 2) slice EM table using the vectors just created 
group_1_em = em[,group_1]
group_2_em = em[,group_2]

#test difference in specified protein abundance btw 2 groups
candidate_gene = "ABRA"
group_1_candidate = as.numeric(group_1_em[candidate_gene,]) #as.numeric transforms df into a numerical vector
group_2_candidate = as.numeric(group_2_em[candidate_gene,])

#perform statistical test

# 1) t-test
t_result = t.test(group_1_candidate, group_2_candidate, alternative = "two.sided")

#note: a two-tailed test uses both the positive and negative tails of the distribution,
# so, it tests for the possibility of positive or negative differences. 
# a one-tailed test is appropriate if you only want to determine if there is 
# a difference between groups in a specific direction.

t_result$p.value

# 2) wilcox rank sum test (non-parametric alternative to t-test)
wilcox.test(group_1_candidate, group_2_candidate, alternative = "two.sided")


# tests in a LOOP: minimal code approach 
#what genes do you want to test? 
genes_to_test = c("PAPPA","KLK3","CCL21")

for(gene in genes_to_test){
  #gets expression data for each gene for each group
  group_1_candidate = as.numeric(group_1_em[gene,]) 
  group_2_candidate = as.numeric(group_2_em[gene,])
  
  #does test & calculates p-value
  t_results = t.test(group_1_candidate, group_2_candidate)
  
  #store p-values
  # 1) create df for storage
  de = data.frame(matrix(nrow = length(genes_to_test), ncol=1))
  row.names(de) = genes_to_test
  names(de) = "p-value"
  # 2) store p-value
  de[gene, "p-value"] = w_results$p.value
}

# Test gene expression between more variable groups

#compare two groups divided by the clinical parameter M vs F
#subset em by factor
# 1) get the names of the samples in each group
group_F = subset(ss, SEX == "F")$SAMPLE  #by indicating $SAMPLE, we specify what col of the df to subset
group_M = subset(ss, SEX == "M")$SAMPLE

# 2) slice EM table using the vectors just created 
group_F_em = em[,group_F]
group_M_em = em[,group_M]

#create df to store p-values
de.sex = data.frame(matrix(nrow = length(genes_to_test), ncol=1))
row.names(de.sex) = genes_to_test
names(de.sex) = "p-value"

for(gene in genes_to_test){
  #gets expression data for each gene for each group
  group_F_candidate = as.numeric(group_F_em[gene,]) 
  group_M_candidate = as.numeric(group_M_em[gene,])
  
  #does test & calculates p-value
  t_results = t.test(group_F_candidate, group_M_candidate, alternative="two.sided")
  
  #store p-values
  de.sex[gene, "p-value"] = t_results$p.value
}


# repeat test by subsetting a continuous variable
age_discrete = cut(ss$AGE, 2, labels = c("low","high"))
ss$AGE_GROUP = age_discrete

#compare HIGH vs LOW age groups
# 1) get the names of the samples in each group
age_high = subset(ss, AGE_GROUP == "high")$SAMPLE  #by indicating $SAMPLE, we specify what col of the df to subset
age_low = subset(ss, AGE_GROUP == "low")$SAMPLE

# 2) slice EM table using the vectors just created 
age_high_em = em[,age_high]
age_low_em = em[,age_low]

# create df to store p-values
de.age = data.frame(matrix(nrow = length(genes_to_test), ncol=1))
row.names(de.age) = genes_to_test
names(de.age) = "p-value"

#perform statistical test using FOR loop
for(gene in genes_to_test){
  #gets expression data for each gene for each group
  group_age_high = as.numeric(age_high_em[gene,]) 
  group_age_low = as.numeric(age_low_em[gene,])
  
  #does test & calculates p-value
  t_results = t.test(group_age_high, group_age_low, alternative="two.sided")
  
  #store p-values
  de.age[gene, "p-value"] = t_results$p.value
}


# Now perform tests btw high_eosinophils and low_eosinophils for ALL genes in the EM table using the same loop
de = data.frame(matrix(nrow=nrow(em), ncol=3))
row.names(de) = row.names(em)
names(de) = c("p", "p.adj", "log2fold")

# Loop through each ROW of the em table
# get the protein abundance, perform a test, store the p, p.adj and log2fold in de, then move to the next row.

for (row_number in 1:nrow(em)) {
  # gets the expression data for each group, for the current gene
  group_1_candidate = as.numeric(group_1_em[row_number,]) 
  group_2_candidate = as.numeric(group_2_em[row_number,])
  print(group_1_candidate)
  print(group_2_candidate)
  # calculates P
  t_result = t.test(group_1_candidate, group_2_candidate, alternative = "two.sided")
  # stores P
  de[row_number,"p"] = t_result$p.value
  #calculates p.adj
  de[row_number,"p.adj"] = p.adjust(de[row_number,"p"], method = p.adjust.methods, n=3)
  #calculates log2fold change
  log2fold = log(mean(group_1_candidate),2)-log(mean(group_2_candidate),2)
  #stores log2folds
  de[row_number,"log2fold"] = log2fold
  }

#select sig genes
de_sig = subset(de, p.adj<0.05 & abs(log2fold)>1)
sig_genes = row.names(de_sig)


## Part 2
# Generate a function that takes the 2 groups and returns a p-value (Wilcox) for candidate genes

## use on sample_group
group_1 = subset(ss, SAMPLE_GROUP == "High_Eosinophils")$SAMPLE 
group_2 = subset(ss, SAMPLE_GROUP == "Low_Eosinophils")$SAMPLE

get_p = function(group_1, group_2, cand_genes){
  # slice EM table using the vectors just created 
  group_1_em = em[,group_1]
  group_2_em = em[,group_2]
  
  #create df to store p-values
  de = data.frame(matrix(nrow = length(cand_genes), ncol=3))
  row.names(de) = cand_genes
  names(de) = c("p-value", "p.adj","log2fold")
  
  #perform statistical test using FOR loop
  for(gene in cand_genes){
    #gets expression data for each gene for each group
    group_1_candidate = as.numeric(group_1_em[gene,]) 
    group_2_candidate = as.numeric(group_2_em[gene,])
    #does test & calculates p-value
    w_results = wilcox.test(group_1_candidate, group_2_candidate, alternative="two.sided")
    #store p-values
    de[gene, "p-value"] = w_results$p.value
    #calculates p.adj
    de[gene,"p.adj"] = p.adjust(de[gene,"p-value"], method = p.adjust.methods, n=3)
    #calculates log2fold change
    log2fold = log(mean(group_1_candidate),2)-log(mean(group_2_candidate),2)
    #stores log2folds
    de[gene,"log2fold"] = log2fold
  }
  return(de)
}
get_p(group_1, group_2, cand_genes=c("CCL21","PAPPA","KLK3"))


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

get_de(group_1, group_2, em)

# Part 3

#load library 
library("ggplot2")

# PCA: steps
# scale data
em_scaled = na.omit(data.frame(t(scale(t(em)))))
# run PCA
xx = prcomp(t(em_scaled)) 
pca_coordinates = data.frame(xx$x)
# get % variation
vars = apply(xx$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100 
prop_y = round(vars["PC2"] / sum(vars),4) * 100 
x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="") 
y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
# plot
ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss$SEX)) +
  geom_point() +
  labs(title = "PCA", x= x_axis_label, y= y_axis_label) + 
  theme_bw()

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

# see how samples cluster based on characteristics/colours
make_pc1_pc2(ss$SEX, em)
make_pc1_pc2(ss$AGE, em)
make_pc1_pc2(ss$SAMPLE_GROUP, em)

# function to plot PC3 and PC4
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

make_pc3_pc4(ss$SEX, em)
make_pc3_pc4(ss$SAMPLE_GROUP, em)
make_pc3_pc4(ss$AGE, em)

# Identify genetic signatures associated with components 
summary(pca_coordinates[,1])

#limits 1st and 3rd quartile
q1 = summary(pca_coordinates[,1])[2] #select 2 value of vector output by summary for PC1
q3 = summary(pca_coordinates[,1])[5] #select 5 value

q1_data = subset(pca_coordinates, PC1 < q1) #get the sample names + coordinates of PC1 in 1st quartile
q3_data = subset(pca_coordinates, PC1 > q3) # same but for 3rd quartile

#isolate gene names from table
q1_samples = row.names(q1_data)
q3_samples = row.names(q3_data)

#get the markers, aka get the de for the samples at high end and low end of PC1
de_pc1 = get_de(q1_samples, q3_samples, em)
pc1_markers = subset(de_pc1, p < 0.05)

summary(pca_coordinates[,2])

#component gene function 
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

markers_PC1 = get_component_genes(1, em1)  
markers_PC2 = get_component_genes(2, em)  
markers_PC3 = get_component_genes(3, em)

#access individual items of the list
markers_PC1$markers
markers_PC2$markers


## Part 4
# Heatmaps
# recommended HEATMAP code
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

#slice em to only plot markers
order_em_PC1 = c(markers_PC1$q1_samples, markers_PC1$q3_samples)
order_em_PC2 = c(markers_PC2$q1_samples, markers_PC2$q3_samples)
order_em_PC3 = c(markers_PC3$q1_samples, markers_PC3$q3_samples)

# slice em based on PCA-derived samples (different sample set for each component)
em_sliced_PC1 = em[,order_em_PC1]
em_sliced_PC2 = em[,order_em_PC2]
em_sliced_PC3 = em[,order_em_PC3]

# plot heatmaps of the biomarkers for components 1-3
make_heatmap(em_sliced_PC1, markers_PC1$markers)
make_heatmap(em_sliced_PC2, markers_PC2$markers)
make_heatmap(em_sliced_PC3, markers_PC3$markers)

# tidy up code: one-liner
make_heatmap(em[c(markers_PC1$q1_samples, markers_PC1$q3_samples)], markers_PC1$markers)
make_heatmap(em[c(markers_PC2$q1_samples, markers_PC2$q3_samples)], markers_PC2$markers)
make_heatmap(em[c(markers_PC3$q1_samples, markers_PC3$q3_samples)], markers_PC3$markers)

# Part 5: Single Cell Omics 
#new em with single cell samples
em1 = read.table("EM.csv", header = TRUE, row.names = 1, sep="\t")

#make PCA for comp 1&2 and 3&4
make_pc1_pc2("", em1)
make_pc3_pc4("", em1)

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
#repost get_de function 
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

# make ordered heatmap
make_heatmap(em1[,c(cluster_1, cluster_2, cluster_3)], sig_genes)

