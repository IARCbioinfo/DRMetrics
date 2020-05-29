rm(list=ls())
library(umap) 
library(ade4)
library(FNN)
library(reshape2)
library(ape)
library(ggpubr)
library(foreach)
library(doParallel)
library(gridExtra)
library(eulerr)
library(viridis)
library(DESeq2)
library(RColorBrewer)

# Path to the DRMetrics folder
project_path=""
# Path to output figures
output_path=""

# Load attributes and read counts --------------------------------------------------------
Attributes = read.table(unzip(paste0(project_path,'DRMetrics/data/Attributes.txt.zip'),exdir = paste0(project_path,"zip_extract"))[1], sep = '\t', header = T)
Attributes[1:5,1:10]
dim(Attributes)
summary(Attributes$Histopathology)
summary(Attributes$Molecular_clusters)

read_counts = as.matrix(read.table(unzip(paste0(project_path,'DRMetrics/data/read_counts_all.txt.zip'),exdir = paste0(project_path,"zip_extract"))[1], header = T, row.names=1))
read_counts[1:5,1:10]

# Read gene annotations ---------------------------------------------------
genespans   = read.table(unzip(paste0(project_path,'DRMetrics/data/ref_annot.gtf.gene_spans.zip'),exdir = paste0(project_path,"zip_extract"))[1],stringsAsFactors = F)
dim(genespans)
summary(duplicated(genespans$V1))
rownames(genespans)=genespans[,1]

# Remove genes on sex and mitochondrial chromosomes -----------------------
genenames   = rownames(read_counts)
summary(genenames %in% rownames(genespans))
genes.sex   = genespans[genespans$V2 %in% c("chrM", "chrX", "chrY"), ]
genes.nosex = genespans[!(genespans$V2 %in% c("chrM", "chrX", "chrY")), ]

read_counts = read_counts[which(rownames(read_counts) %in% rownames(genes.nosex)), ]
dim(read_counts)


# Normalization of the read counts ----------------------------------------

colData      = as.matrix(colnames(read_counts))
DESeq_object = DESeqDataSetFromMatrix(countData=read_counts,colData = colData,design = ~1)
VST = varianceStabilizingTransformation(DESeq_object)
VST = assay(VST)
dim(VST)
VST[1:5,1:5]


# Select most variable genes ----------------------------------------------

rvdm     = apply(VST,1,var)
ntopdm   = rev(which( cumsum( sort(rvdm,decreasing = T) )/sum(rvdm) <= 0.5 ))[1]
selectdm = order(rvdm, decreasing = TRUE)[seq_len(ntopdm)]
vst50    = t(VST[selectdm,])
vst50[1:5,1:5]
dim(vst50)

# keep only non-duplicated samples ----------------------------------------

vst50 = vst50[as.character(Attributes$Sample_ID),]
dim(vst50)


# Run PCA on the vst data  ------------------------------------------------
pca5D = dudi.pca(vst50, center = T , scale = T,  scannf = F , nf=5 )$li
head(pca5D)

# Generate the figure representing the PCA axes ---------------------------
nb.cols <- 13
distinctive_cols <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

Axis_list = list(c(2,1), c(3,1), c(4,1), c(5,1), c(3,2), c(4,2), c(5,2), c(4,3), c(5,3), c(5,4))

plist <- list()
for(i in 1:length(Axis_list)){
  x_index = Axis_list[[i]][1]
  y_index = Axis_list[[i]][2]
  xlab   = paste("PC", x_index)
  ylab   = paste("PC", y_index)
  c_data = data.frame("axisX"=  pca5D[,x_index], "axisY" =pca5D[,y_index] )
  p      = ggplot(c_data, aes(x = axisX , y = axisY, color = Attributes$Molecular_clusters)) + 
    geom_point(size=2, alpha=0.8) + scale_colour_manual(values=distinctive_cols) 
  
  if (x_index == 2 & y_index == 1){
    plist[[i]] = p + theme( legend.position="right", legend.title= element_blank()) + labs(x = xlab, y = ylab) 
  }else{
    plist[[i]] = p + theme( legend.position="none", legend.title= element_blank()) + labs(x = xlab, y = ylab) 
  }
}

pdf(paste0(output_path,"PCA_5axes.pdf"))
ggarrange(plist[[1]])
grid.arrange(grobs = lapply(2:10, function(i){plist[[i]]}),  
             layout_matrix = rbind(c(2:4),c(5:7),c(8:10)), nrow =2)
dev.off()

# Run UMAP with nn=238 ----------------------------------------------------
# run umap from scratch
UMAP_model=umap(vst50, n_neighbors = 238)
UMAP_TM = data.frame(UMAP_model$layout)
colnames(UMAP_TM)=c("UMAP1","UMAP2")
UMAP_TM <- cbind(UMAP_TM, 'Molecular_clusters' = Attributes$Molecular_clusters)
rownames(UMAP_TM)=Attributes$Sample_ID

# Note that the next two lines replace the UMAP coordinates calculated in the previous cell with our pre-computed coordinates. You should comment them if you are running a customized analysis.
UMAP_TM = read.table(unzip(paste0(project_path,'DRMetrics/data/Coords_umap_nn238.txt.zip'),exdir = paste0(project_path,"zip_extract"))[1], sep = '\t', header = T, row.names=1)
UMAP_TM <- cbind(UMAP_TM, 'Molecular_clusters' = Attributes$Molecular_clusters)

# Generate the figure representing the UMAP run with nn=238 ---------------

pdf(paste0(output_path,"UMAP_nn238.pdf"))
ggplot(UMAP_TM, aes(x = UMAP1 , y = UMAP2, color = Molecular_clusters))+
  geom_point(size=4, alpha=0.8)+ scale_colour_manual(values=distinctive_cols)+ #scale_color_brewer(palette="Spectral")+
  labs(x = "UMAP dimension 1" , y = "UMAP dimension 2", color = "")+ theme(legend.position = "bottom")+
  guides(col = guide_legend( nrow = 4)) 
dev.off()

# Validation of biological hypotheses -------------------------------------

clusters_distances_F <- function(Coords, evals, ref1, ref2, alt_hyp){
  # Coords: a data frame with the following columns: i) the samples coordinates on the dimensionality reduction axes, and ii) the molecular cluster associated to each sample
  # evals: the molecular cluster to analyse 
  # ref1: the supposed molecular cluster(s)
  # ref2: the supposed histopathological group(s)
  # alt_hyp: the alternative hypothesis for the wilcoxon test 
  Coords_evals = Coords[Coords$Molecular_clusters==evals,]
  T1           = Coords[which(Coords$Molecular_clusters %in% ref1),] 
  T1_centroid  = colMeans(T1[,1:2])
  dist_1       = sqrt( (Coords_evals[,1] - T1_centroid[1])^2 + ( Coords_evals[,2] - T1_centroid [2])^2 )
  T2           = Coords[which(Coords$Molecular_clusters %in% ref2),]
  T2_centroid  = colMeans(T2[,1:2])
  dist_2       = sqrt( (Coords_evals[,1] - T2_centroid[1])^2 + ( Coords_evals[,2] - T2_centroid [2])^2 )
  return(wilcox.test(dist_1 , dist_2, alt_hyp)$p.value) 
}
paste('Wilcoxon p-value = ', format(clusters_distances_F(UMAP_TM, "SCLC/LCNEC-like",c("LCNEC/TypeII","LCNEC/TypeI"),"SCLC/SCLC-like", 'less'), digits = 2)) 

paste('Wilcoxon p-value = ', format(clusters_distances_F(UMAP_TM, "LCNEC/SCLC-like","SCLC/SCLC-like", c("LCNEC/TypeI", "LCNEC/TypeII"), 'less'), digits = 2))

paste('Wilcoxon p-value = ', format(clusters_distances_F(UMAP_TM, "LCNEC/TypeI","LCNEC/TypeII","SCLC/SCLC-like", 'less'), digits = 2)) 

paste('Wilcoxon p-value = ', format(clusters_distances_F(UMAP_TM, "SCLC/LCNEC-like","LCNEC/TypeII","LCNEC/TypeI", 'less'), digits = 2)) 

paste('Wilcoxon p-value = ', format(clusters_distances_F(UMAP_TM, "Supra_carcinoid", c("LCNEC/TypeII","LCNEC/TypeI", 'SCLC/LCNEC-like'), c("SCLC/SCLC-like", "LCNEC/SCLC-like"), 'less'), digits = 2)) 


# Evaluation of the individual preservation the samples neighborhood ----------------------------------

source(paste0(project_path,"DRMetrics/scripts/DR_quality_metrics.r") )

distpca5D = as.matrix(dist(pca5D, method = "euclidean", diag = TRUE, upper = TRUE))
distpca5D[1:5,1:5]

k_values = seq(2,238,5)
SD_pca5D = compute_SD_allSamples(distRef = distpca5D, List_projection = list("UMAP_TM" = UMAP_TM[,1:2]),
                                 k_values = k_values,colnames_res_df = "UMAP-nn238")
head(SD_pca5D)

sd_map_lnen = SD_map_f(SD_pca5D , UMAP_TM, "bottom")
sd_map_res  = sd_map_lnen[[1]]
head(sd_map_res)
sd_map_res[order(sd_map_res$SD,decreasing = T),]

pdf(paste0(output_path,"individual_SD_umap.pdf"))
sd_map_lnen[[2]]
dev.off()

plot_dist <- function(coords,projection_coords,focal_sample){
  dist1 = as.matrix(dist(coords[, 1:(ncol(coords)-1)], method = "euclidian", diag = TRUE, upper = TRUE))
  Sample_names = rownames(coords)
  dist1        = data.frame("Sample_ID" = rownames(dist1)[order(dist1[,focal_sample])], "rank" = 1:length(rownames(dist1)) )
  merged_df    = merge(cbind(Sample_ID = rownames(projection_coords), projection_coords), dist1,  by= "Sample_ID")
  colnames(merged_df)[2:3] = c("V1","V2")
  
  proj_dist    = ggplot(merged_df, aes(x=V1, y=V2,  color=rank, label=merged_df$Sample_ID)) +
    geom_point(size=2, alpha =.8) + scale_color_distiller(palette = "Spectral") +
    geom_point(aes(merged_df$V1[which(merged_df$Sample_ID == focal_sample)],
                   merged_df$V2[which(merged_df$Sample_ID == focal_sample)] ),
               size=2, col = "black") +
    geom_text(aes(merged_df$V1[ merged_df$Sample_ID== focal_sample],
                  merged_df$V2[ merged_df$Sample_ID== focal_sample],
                  label=merged_df$Sample_ID[merged_df$Sample_ID==focal_sample]),
              col = 'black' , hjust=0, vjust=0 ) +
    labs(title="", y="dim2", x="dim1") + theme( legend.position = "bottom")
  return(proj_dist)
}

pdf(paste0(output_path,"S02340_SD.pdf"))
p1 <- plot_dist(UMAP_TM,UMAP_TM,"S02340")
p2 <- plot_dist(pca5D,UMAP_TM,"S02340")
grid.arrange(grobs = list(p1,p2), layout_matrix = rbind(c(1,2)), nrow =1)
dev.off()

# Evaluation of the global preservation the samples neighborhood ----------------------------------

n_neighborsL <- c(15,30,90,180,238)
List_projection <- list("pca2D" = pca5D[1:2], "pca5D" = pca5D) 
for (j in 1:length(n_neighborsL)){
  umap_c = umap(vst50, n_neighbors = n_neighborsL[j])$layout
  List_projection[[j+2]] <- umap_c
  names(List_projection)[j+2] <- paste("umap_nn", n_neighborsL[j])
}

distRef = as.matrix(dist(vst50, method = "euclidean", diag = TRUE, upper = TRUE))
distRef[1:5,1:5]

k_values        = seq(2,238,5)
n_neighborsL    = c(15,238)
colnames_res_df = c('PCA-2D', 'PCA-5D',  paste0('UMAP-nn',n_neighborsL) )
n_iter          = 1000
no_cores        = 20


cl = makeCluster(no_cores)
registerDoParallel(cl)
SD_means <-  foreach(i=1:n_iter, .packages = c("umap","doParallel"))%dopar%{
  List_projection = list("PCA-2D" = pca5D[,1:2], "PCA-5D" = pca5D) 
  for (j in 1:length(n_neighborsL)){
    List_projection[[j+2]] = umap(vst50, n_neighbors = n_neighborsL[j])$layout
    names(List_projection)[j+2] = paste("umap_nn", n_neighborsL[j])
  } 
  SD_comp = compute_SD_allSamples(distRef, List_projection, k_values, colnames_res_df)
  
  # For each k compute the mean value over all samples
  SD_k_mean = sapply(1:length(colnames_res_df), function(j) aggregate(SD_comp[,colnames_res_df[j]],list(SD_comp$K),mean)$x) 
  colnames(SD_k_mean) = colnames_res_df
  rownames(SD_k_mean) = k_values
  
  return(SD_k_mean)
}  
stopCluster(cl)


Mean_SD = melt( apply(simplify2array(SD_means), 1:2, mean), varnames = c("k","DRM"), value.name = "mean" )
sd_SD   = melt( apply(simplify2array(SD_means), 1:2, sd), varnames = c("k","DRM"), value.name = "sd" )
Main_SD_DF = cbind(Mean_SD, sd = sd_SD$sd)
head(Main_SD_DF)

pdf(paste0(output_path,"SD_noNorm.pdf"))
theme_set(theme_bw())
ggplot(Main_SD_DF, aes(x =  k, y = mean, color = DRM)) +
  geom_line() +  geom_point()+
  scale_color_viridis(discrete=TRUE) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd))+
  labs(title = "", y = "mean(SD)", x = "k") 
dev.off()


# Normalization of SD values 
SD_means_norm = lapply( SD_means, function(x) (x[,2] - x)/(x[,2] - x[,1] ) )

Mean_SD_norm = melt( apply(simplify2array(SD_means_norm), 1:2, mean), varnames = c("k","DRM"), value.name = "mean" )
sd_SD_norm   = melt( apply(simplify2array(SD_means_norm), 1:2, sd), varnames = c("k","DRM"), value.name = "sd" )
Main_SD_norm_DF = cbind(Mean_SD_norm, sd = sd_SD_norm$sd)

pdf(paste0(output_path,"SD_Norm.pdf"),w=9,h=7)
theme_set(theme_bw())
ggplot(Main_SD_norm_DF, aes(x = k, y = mean, color = DRM)) +
  geom_line() + geom_point()+
  scale_color_viridis(discrete=TRUE) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd))+
  labs(title = "", y = "mean(SD')", x = "k") +
  theme( axis.title.x=element_text(size=20),  # X axis title
         axis.title.y=element_text(size=20),  # Y axis title
         axis.text.x=element_text(size=16),  # X axis text
         axis.text.y=element_text(size=16),
         legend.text = element_text(size = 20) ,
         legend.title = element_blank())
dev.off()

# Pvalue for the SD metric ---------------------------------------------------------------
pwt_df=Mean_SD
colnames(pwt_df)=c("K","method","mean_seq")

if (length(k_values)< 30){
  paired_test_m  <- pairwise.wilcox.test(pwt_df$mean_seq, pwt_df$method, paired = TRUE)$p.value #p.adj = "holm",
}else{
  paired_test_m  <-  pairwise.t.test(pwt_df$mean_seq, pwt_df$method, paired = TRUE)$p.value #p.adj = "holm",
}
paired_test_m

# Compute MI values in the different spaces -------------------------------

n_genes = ncol(vst50)
ksep=25
k_step  = seq(2,238,ksep)

MI_GeneExprMat = moran_I_knn(vst50[,1:n_genes], vst50,k_step)
MI_UMAP        = moran_I_knn(vst50[,1:n_genes], UMAP_TM[,1:2],k_step)
MI_PCA5D       = moran_I_knn(vst50[,1:n_genes], pca5D,k_step)

head(MI_UMAP)


# Get overlap between the top ranked genes --------------------------------
# Compute mean MI values over k levels for each gene and order the genes based on their averaged MI values

MI_PCA5D_order       = order(colMeans(MI_PCA5D[,1,]), decreasing=T) 
MI_UMAP_order        = order(colMeans(MI_UMAP[,1,]), decreasing=T) 
MI_GeneExprMat_order = order(colMeans(MI_GeneExprMat[,1,]), decreasing=T) 

PCA5D_GeneExprMat_overlap = sapply(1:n_genes, function(i) mean( MI_PCA5D_order[1:i] %in% MI_GeneExprMat_order[1:i] ) )
PCA5D_UMAP_overlap        = sapply(1:n_genes, function(i) mean( MI_PCA5D_order[1:i] %in% MI_UMAP_order[1:i] ) )
UMAP_GeneExprMat_overlap  = sapply(1:n_genes, function(i) mean( MI_UMAP_order[1:i] %in% MI_GeneExprMat_order[1:i] ) )
PCA5D_UMAP_GeneExprMat_overlap = sapply(1:n_genes, function(i) mean( MI_PCA5D_order[1:i] %in% MI_GeneExprMat_order[1:i] & MI_PCA5D_order[1:i] %in% MI_UMAP_order[1:i] ) )
Expected = sapply(1:n_genes, function(i) (i/n_genes) )

t.test(PCA5D_GeneExprMat_overlap, UMAP_GeneExprMat_overlap, paired = T, alternative = 'greater')

cols <- c("Expected" =  "#FDE725FF", "PCA 5D - GeneExprMat" = "#440154FF","PCA 5D - UMAP 238 - GeneExprMat" = "#5DC863FF",  "PCA 5D - UMAP 238" = "#21908CFF", "UMAP 238 - GeneExprMat" = "#3B528BFF")

pdf(paste0(output_path,"MI_fig_",ksep,".pdf"),w=9,h=7)
theme_set(theme_bw())
ggplot()+
  geom_point(aes(x=seq(1:n_genes), y=PCA5D_GeneExprMat_overlap , colour="PCA 5D - GeneExprMat"), size = 0.1)+
  geom_point(aes(seq(1:n_genes),  y=UMAP_GeneExprMat_overlap, colour="UMAP 238 - GeneExprMat"), size = 0.1)+
  geom_point(aes(seq(1:n_genes),  y=PCA5D_UMAP_overlap , colour="PCA 5D - UMAP 238"), size = 0.1)+
  geom_point(aes(seq(1:n_genes),  y=PCA5D_UMAP_GeneExprMat_overlap, colour="PCA 5D - UMAP 238 - GeneExprMat"), size = 0.1)+
  labs(title="", y="Overlapping %", x="N", caption = "") +
  scale_colour_manual(values=cols, labels =c("Expected", "PCA-5D vs Original","PCA-5D vs UMAP-nn238 vs Original",  "PCA-5D vs UMAP-nn238", "UMAP-nn238 vs Original")) +
  theme( legend.position = c(.95, .05),
         legend.justification = c("right", "bottom"))+
  geom_line(aes(seq(1:n_genes),  y=Expected, colour="Expected"), size = 1)+
  geom_vline(xintercept = 1000, linetype="dashed")+
  theme( axis.title.x=element_text(size=20),  # X axis title
         axis.title.y=element_text(size=20),  # Y axis title
         axis.text.x=element_text(size=15),  # X axis text
         axis.text.y=element_text(size=15),
         legend.text = element_text(size = 15),
         legend.title = element_blank())
dev.off()

# Plot Euler diagram  -----------------------------------------------------

n_genes_sel = 1000
abc = PCA5D_UMAP_GeneExprMat_overlap[n_genes_sel]*n_genes_sel
ab  = PCA5D_UMAP_overlap[n_genes_sel]*n_genes_sel
ac  = PCA5D_GeneExprMat_overlap[n_genes_sel]*n_genes_sel
bc  = UMAP_GeneExprMat_overlap[n_genes_sel]*n_genes_sel

Venn <- c(c("A" = n_genes_sel - ((ab -abc) + abc + (ac -abc)) ,
            "B"  = n_genes_sel - ((ab -abc) + abc + (bc -abc)),
            "C"  = n_genes_sel - ((ac -abc) + abc + (bc -abc)),
            "A&B" = ab -abc,
            "A&C" = ac -abc,
            "B&C" = bc -abc ,
            "A&B&C" = abc))

pdf(paste0(output_path,"venn_ksep_",ksep,".pdf"),w=5,h=5)
PlotEuler <- plot(euler(Venn), fills = c("#5DC863FF", "#FDE725FF", "#21908CFF"),  quantities = list(cex = 1), alpha= .7, cex = 2, legend = list(labels = c("PCA 5D", "UMAP nn = 238", "GenesExprMat"), alpha=.7 ,cex =1.2))
PlotEuler 
dev.off()