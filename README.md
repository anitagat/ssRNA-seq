# Biomarker signatures of bulk-RNA seq and ssRNA-seq data
Analysis of biomarker signatures of bulk RNA-seq and ssRNA-seq data.

# Contents:
### Part 1 (code in file: biomarkersignatures.R)
- Subset expression matrix (obtained through RNA-seq, metabolomics, ChIP-seq) by clinical features, e.g. factors or continuous variables.
- Perform statistical tests (t-test & wilcox) singularly or in series using loops, for a few genes or all genes
Generate de tables with p values, p adj and log2fold 

### Part 2 (code in file: biomarkersignatures.R)
- **get_p** function: calculates the p-value and log2fold of a t-test or a Wilcox test for single or multiple genes between two predefined groups
- **get_de** function: calculates the p-value and log2fold of a t-test or a Wilcox test for a whole expression matrix between two predefined groups

### Part 3 (code in file: biomarkersignatures.R)
Identify biomarker signatures using PCA components and quartiles. 
- **make_pc1_pc2** function
- **make_pc3_pc4** function
- **get_component_genes** function

### Part 4 (code in file: biomarkersignatures.R)
Plot biomarker signatures using heatmaps.
- **make_heatmap function**

### Part 5 (code in file: ssRNA-se.R)
Single cell RNA-seq analysis with K-mean clustering (uses all functions created in biomarkersignatures.R Parts 1-4)

The same functions can be used both for bulk RNA-seq and scRNA-seq analyses. 

### Example study:
![](plots.png)


