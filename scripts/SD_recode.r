

# Functions ---------------------------------------------------------------

compute_SD_new=function(dist_space1,dist_space2,k){
  
  ranks_df=data.frame(Sample_ID=names(dist_space1),Space1=rank(dist_space1),Space2=rank(dist_space2),stringsAsFactors = F)
  k_subset=ranks_df[order(ranks_df$Space1),][2:k,]
  s1=sum((k-k_subset$Space1)*abs(k_subset$Space1 - k_subset$Space2))
  k_subset=ranks_df[order(ranks_df$Space2),][2:k,]
  s2=sum((k-k_subset$Space2)*abs(k_subset$Space1 - k_subset$Space2))
  SD_k=s1/2 + s2/2
  return(SD_k)
}
compute_SD_emilie=function(dist_space1,dist_space2,k){
  
  N1_dist_l <- sort(dist_space1)
  N2_dist_l <- sort(dist_space2)
  
  N1_df <- data.frame("Sample_ID" = names(N1_dist_l)[1:k] , "Rank1" = 1:k, stringsAsFactors = F)
  N2_df <- data.frame("Sample_ID" = names(N2_dist_l)[1:k] , "Rank2" = 1:k, stringsAsFactors = F)
  
  # ranks_df=merge(N1_df,N2_df,all.x=T,all.y=T)
  # ranks_df[is.na(ranks_df)]=0
  # # ranks_df=data.frame(Sample_ID=names(dist_space1),Space1=rank(dist_space1),Space2=rank(dist_space2),stringsAsFactors = F)
  # k_subset=ranks_df[order(ranks_df$Space1),][2:k,]
  # s1=sum((k-k_subset$Space1)*abs(k_subset$Space1 - k_subset$Space2))
  # k_subset=ranks_df[order(ranks_df$Space2),][2:k,]
  # s2=sum((k-k_subset$Space2)*abs(k_subset$Space1 - k_subset$Space2))
  # SD_k=s1/2 + s2/2
  
  All_neighbors <- unique(c(N1_df$Sample_ID,N2_df$Sample_ID))
  s1 = 0
  s2 = 0
  for (j in 1:length( All_neighbors)){
    if (All_neighbors[j] %in%  N1_df$Sample_ID & All_neighbors[j] %in%  N2_df$Sample_ID ){
      N1_index_j = which(N1_df$Sample_ID  == All_neighbors[j]  )
      N2_index_j = which(N2_df$Sample_ID  == All_neighbors[j]  )
      s1 = s1 + ((k - N1_df$Rank1[N1_index_j]) * abs(N1_df$Rank1[N1_index_j] - N2_df$Rank2[N2_index_j]))
      s2 = s2 + ((k - N2_df$Rank2[N2_index_j]) * abs(N1_df$Rank1[N1_index_j] - N2_df$Rank2[N2_index_j]))
    }else if (All_neighbors[j] %in%  N1_df$Sample_ID){ 
      N1_index_j = which(N1_df$Sample_ID  == All_neighbors[j]  )
      s1 = s1 + ((k - N1_df$Rank1[N1_index_j]) * abs(N1_df$Rank1[ N1_index_j]))
      s2 = s2 
    }else{ 
      N2_index_j = which(N2_df$Sample_ID  == All_neighbors[j]  )
      s1 = s1 
      s2 = s2 + ((k - N2_df$Rank2[N2_index_j]) * abs( N2_df$Rank2[N2_index_j]))
    }
  }
  S = 0.5 * s1 + 0.5 * s2
  return(S)
}
compute_SD_allSamples=function(distRef,List_projection,k_values,colnames_res_df){
  SD_comp=data.frame(Sample_ID=rep(rownames(distRef),each=length(k_values)),K=rep(k_values,nrow(distRef)))
  for(df in 1:length(List_projection)){
    
    space1 = List_projection[[df]]
    
    dist1 = as.matrix(dist(space1[, 2:dim(space1)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
    rownames(dist1) = as.character(space1[ ,1])
    colnames(dist1) = as.character(space1[ ,1])
    # check that the two matrices of distances have the same ordering 
    if(!all(rownames(dist1)==rownames(distRef))){
      print("Error: The two matrices of distances do not have the same ordering!")
    }
    
    # For each sample compute SD metric for each k value
    cl = makeCluster(2)
    registerDoParallel(cl)
    samples_names = rownames(dist1)
    SD_k_samples = foreach(i=1:nrow(dist1),.combine=rbind,.export=c("compute_SD")) %dopar% {
      SD_k_values = sapply(k_values,function(k) compute_SD(dist1[i, ],distRef[i, ],k) ) 
      data.frame(Sample_ID = rep(samples_names[i], length(k_values)),K=k_values,SD = SD_k_values)
    }
    stopCluster(cl)
    
    # merge the 
    SD_comp=cbind(SD_comp,SD_k_samples[,3])
    # SD_comp=merge(SD_comp,SD_k_samples,by="Sample_ID")
    colnames(SD_comp)[ncol(SD_comp)]=colnames_res_df[df]
  }
  return(SD_comp)
}

compute_SD_v1=F
if(compute_SD_v1){
  compute_SD=compute_SD_emilie
  out_folder="SD_compute_v1"
}else{
  compute_SD=compute_SD_new
  out_folder="SD_compute_v2"
}
# Evaluation of individual samples neighborhood ---------------------------

dataRef = pca5D
distRef = as.matrix(dist(dataRef[, 2:ncol(dataRef)], method = "euclidian", diag = TRUE, upper = TRUE))
rownames(distRef) = as.character(dataRef[ ,1])  # rownames(dataRef)
colnames(distRef) = as.character(dataRef[ ,1])
distRef[1:5,1:5]

k_values = seq(2,208,5)
Main_SQ_res_pca = compute_SD_allSamples(distRef = distRef, List_projection = list("UMAP_TM" = UMAP_TM[,1:3]), k_values = k_values,colnames_res_df = "UMAP-nn208")

sd_map_lnen = SD_map_f(Main_SQ_res_pca , UMAP_TM, "bottom")
sd_map_res = sd_map_lnen[[1]]
head(sd_map_res)
lowest_SD = sd_map_res$Sample_ID[order(sd_map_res$SD, decreasing=T)]
lowest_SD[1:5]

pdf(paste0(project_path, "results/UMAP_figures/", out_folder, "/indiv_SD_PCA5DvsUMAP208.pdf"))
sd_map_lnen[[2]]
dev.off()


# Evaluation of global samples neighborhood -------------------------------

dataRef = vst50
distRef = as.matrix(dist(dataRef[, 2:ncol(dataRef)], method = "euclidian", diag = TRUE, upper = TRUE))
rownames(distRef) = as.character(dataRef[ ,1])  # rownames(dataRef)
colnames(distRef) = as.character(dataRef[ ,1])
distRef[1:5,1:5]

k_values = seq(2,208,5)
n_neighborsL = c(15,208)
colnames_res_df = c('PCA-2D', 'PCA-5D',  sapply(n_neighborsL, function(i) paste0('UMAP-nn',i)) )
n_iter = 100

start_time = Sys.time()
no_cores = 30
cl_par = makeCluster(no_cores)
registerDoParallel(cl_par)

Main_SQ_REP <-  foreach(i=1:n_iter,.packages = c("umap","doParallel"))%dopar%{
  List_projection <- list("pca2D" = pca5D[1:3], "pca5D" = pca5D)
  for (j in 1:length(n_neighborsL)){
    umap_c <- umap(vst50[,2:dim(vst50)[2]], n_neighbors = n_neighborsL[j])
    umap_c <- as.data.frame(umap_c$layout)
    umap_c <- cbind("Sample_ID" = vst50[,1], umap_c)
    umap_c <- umap_c[order(umap_c$Sample_ID),]
    names_umap <- paste("umap_nn", n_neighborsL[j])
    List_projection[[j+2]] <- umap_c
    names(List_projection)[j+2] <- names_umap 
  }
  
  SD_comp=compute_SD_allSamples(distRef,List_projection,k_values,colnames_res_df)
  
  # For each k compute the mean value over all samples
  SD_k_mean=sapply(1:length(colnames_res_df), function(j) aggregate(SD_comp[,colnames_res_df[j]],list(SD_comp$K),mean)$x) 
  SD_k_mean=as.data.frame(SD_k_mean)
  colnames(SD_k_mean)=colnames_res_df
  SD_k_mean$K=k_values
  
  list(all_SD_comp=SD_comp, SD_k_mean=SD_k_mean)
  
}  

stopCluster(cl_par)
end_time = Sys.time()
start_time - end_time

# merge iterations outputs
SD_comp = foreach(fold.result=Main_SQ_REP, fold.num=icount(), .combine=c) %do%{
  fold.result$all_SD_comp
}
SD_k_mean = foreach(fold.result=Main_SQ_REP, fold.num=icount(), .combine=rbind) %do%{
  fold.result$SD_k_mean
}

# Nextjournal mean sd no norm -----------------------------------------------------

Mean_SD = melt(aggregate(SD_k_mean[,colnames_res_df], list(SD_k_mean$K), mean), id.vars = "Group.1")
colnames(Mean_SD)[3] = "mean"
colnames(Mean_SD)[1] = "K"
sd_SD = melt(aggregate(SD_k_mean[,colnames_res_df], list(SD_k_mean$K), sd), id.vars = "Group.1")
colnames(sd_SD)[3] = "sd"
Main_SD_DF <- cbind(Mean_SD, "sd" = sd_SD$sd)

pdf(paste0(project_path,"results/UMAP_figures/", out_folder, "/SD_noNorm.pdf"))
ggplot(Main_SD_DF, aes(x =  K, y = mean, color = variable)) +
  geom_line() +
  geom_point()+
  scale_color_viridis(discrete=TRUE) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd))+
  labs(title = "", y = "mean(SD)", x = "k") 
dev.off()

# Nextjournal mean sd norm -----------------------------------------------------

MEAN_SD_DF_norm = data.frame("K" = SD_k_mean$K)
for (i in which(colnames(SD_k_mean) != "K")){
  MEAN_SD_DF_norm[[colnames(SD_k_mean)[i]]] = (SD_k_mean$`PCA-5D` - SD_k_mean[,i])/(SD_k_mean$`PCA-5D` - SD_k_mean$`PCA-2D`)
}

Mean_SDn = melt(aggregate(MEAN_SD_DF_norm[, colnames_res_df], list(MEAN_SD_DF_norm$K), mean), id.vars = "Group.1")
colnames(Mean_SDn)[3] = "mean"
colnames(Mean_SDn)[1] = "K"
sd_SDn = melt(aggregate(MEAN_SD_DF_norm[, colnames_res_df], list(MEAN_SD_DF_norm$K), sd), id.vars = "Group.1")
colnames(sd_SDn)[3] = "sd"
Main_SD_DFn <- cbind(Mean_SDn, "sd" = sd_SDn$sd)


pdf(paste0(project_path,"results/UMAP_figures/", out_folder, "/SD_Norm.pdf"),w=9,h=7)
theme_set(theme_bw())
ggplot(Main_SD_DFn, aes(x =  K, y = mean, color = variable)) +
  geom_line() +
  geom_point()+
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

# Pvalue sd ---------------------------------------------------------------
SD_k_mean_iter = aggregate(SD_k_mean[,colnames_res_df], list(SD_k_mean$K), mean)
pwt_df=melt(SD_k_mean_iter, id.vars="Group.1")
colnames(pwt_df)=c("K","method","mean_seq")

if (nrow(SD_k_mean_iter)< 30){
  paired_test_m  <- pairwise.wilcox.test(pwt_df$mean_seq, pwt_df$method,    paired = TRUE)$p.value #p.adj = "holm",
}else{
  paired_test_m  <-  pairwise.t.test(pwt_df$mean_seq, pwt_df$method,    paired = TRUE)$p.value #p.adj = "holm",
}
paired_test_m

# to_remove ---------------------------------------------------------------


SD_map_f <- function(SD_df, Coords_df, legend_pos = "right" ){
  
  # SD_df data frame such as :
  #  col1 = Sample_ID, col2 = k, col3 = SD_values
  # Coords_df data frame such as :
  # col1 = Sample_ID, col2 = AxisX, col3 = AxisY
  
  colnames(SD_df) <- c("Sample_ID", "k", "SD")
  SD_df <- SD_df[order(SD_df$Sample_ID),]
  SD_df$Sample_ID <- as.character( SD_df$Sample_ID)
  unique_sample_id <- as.character(unique(SD_df$Sample_ID)) 
  SDmeans_ID <- unlist(lapply(1:length(unique_sample_id), function(i){
    mean(SD_df$SD[which(SD_df$Sample_ID == unique_sample_id [i])])
  }))
  colnames(Coords_df) <- c("Sample_ID", "V1", "V2")
  Coords_df <- Coords_df[order(Coords_df$Sample_ID),]
  SD_Coords_df <- cbind(  Coords_df, "SD" = SDmeans_ID)
  SD_Coords_df <- SD_Coords_df[order(SD_Coords_df$SD, decreasing = T),]
  pSD_Main <- ggplot(  SD_Coords_df, aes(x=V1, y=V2,  color=SD/max(SD))) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
  pSD_Main <- pSD_Main +  labs(title="", 
                               y="dim2", x="dim1") +
    theme( legend.position = legend_pos,
           plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
           plot.subtitle =element_text(size=14, hjust=0.5),
           plot.caption =element_text(size=12,  hjust=0.5),
           axis.title.x=element_text(size=16),  # X axis title
           axis.title.y=element_text(size=16),  # Y axis title
           axis.text.x=element_text(size=14),  # X axis text
           axis.text.y=element_text(size=14),
           legend.text = element_text(size = 10) ,
           legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text
  
  return(list(SD_Coords_df,pSD_Main))
  
  
}
