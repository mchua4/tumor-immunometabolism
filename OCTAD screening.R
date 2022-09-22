################
# Dependencies #
################ 

packages=c("dplyr", "magrittr", "ggplot2", "doParallel", "foreach", "lme4", 
           "Rfast", "httr", "data.table")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

bioconductor_packages = c("RUVSeq", "edgeR", "DESeq2", "limma", "rhdf5", "artMS")

if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager"))
    install.pacakages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

devtools::install_github("Bin-Chen-Lab/octad")

library(octad.db)

############################
# Case and control samples #
############################

# Select case samples
head(phenoDF)
HCC_primary = subset(phenoDF, cancer=="liver hepatocellular carcinoma" & sample.type=="primary")
case_id = HCC_primary$sample.id
HCC_with_TP53_primary = subset(phenoDF, cancer=="liver hepatocellular carcinoma" & sample.type=="primary" & grepl("TP53", mutation_list))

# Select control samples
control_id = computeRefTissue(case_id, output=T, adjacent=T, source="octad", control_size=50)
HCC_adjacent = subset(phenoDF, cancer=="liver hepatocellular carcinoma" & sample.type=="adjacent" & data.source=="TCGA")
control_id = HCC_adjacent$sample.id

# t-SNE matrix
tsne$type <- "others"
tsne$type[tsne$sample.id %in% case_id] <- "case"
tsne$type[tsne$sample.id %in% control_id] <- "control"

# Plot
(p2 <- ggplot(tsne, aes(X, Y, color=type)) + geom_point(alpha=0.4)+labs(title=paste("TSNE PLOT"), x="TSNE Dim1", y="TSNE Dim2", caption="OCTAD") + theme_bw())
