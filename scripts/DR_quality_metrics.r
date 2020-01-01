
compute_SD <- function(dist_space1,dist_space2,k){
  #computes SD for a single given point
  #dist_space1: distances in space 1
  #dist_space2: distances in space 2
  #k: number of neighbors considered
  ranks_df = data.frame(Sample_ID=names(dist_space1), Space1=rank(dist_space1)-2, Space2=rank(dist_space2)-2, stringsAsFactors = F)
  k_subset = ranks_df[order(ranks_df$Space1),][2:k,]
  s1       = sum((k-k_subset$Space1)*abs(k_subset$Space1 - k_subset$Space2))
  k_subset = ranks_df[order(ranks_df$Space2),][2:k,]
  s2       = sum((k-k_subset$Space2)*abs(k_subset$Space1 - k_subset$Space2))
  return( s1/2 + s2/2)
}

compute_SD_allSamples <- function(distRef,List_projection,k_values,colnames_res_df, threads=2){
  #computes SD values for all samples and for multiple comparisons (each reduced space listed in List_projection compared to the reference space) 
  #distRef: distances in the reference space
  #List_projection: list of reduced spaces to compare to the reference space
  #k_values: vector listing the k values corresponding to the number of neighbors considered
  SD_comp=data.frame(Sample_ID=rep(rownames(distRef),each=length(k_values)),K=rep(k_values,nrow(distRef)))
  for(df in 1:length(List_projection)){
    space1 = List_projection[[df]]
    dist1 = as.matrix(dist(space1, method = "euclidean", diag = TRUE, upper = TRUE))
    
    # check that the two matrices of distances have the same ordering 
    if(!all(rownames(dist1)==rownames(distRef))){
      print("Error: The two matrices of distances do not have the same ordering!")
    }
    # For each sample compute SD metric for each k value
    cl = makeCluster(threads)
    registerDoParallel(cl)
    samples_names = rownames(dist1)
    SD_k_samples  = foreach(i=1:nrow(dist1), .combine=c, .export=c("compute_SD")) %dopar% {
      SD_k_values = sapply(k_values,function(k) compute_SD(dist1[i, ],distRef[i, ],k) )
      SD_k_values
    }
    stopCluster(cl)
    
    # merge 
    SD_comp=cbind(SD_comp,SD_k_samples)
    colnames(SD_comp)[ncol(SD_comp)]=colnames_res_df[df]
  }
  return(SD_comp)
}

SD_map_f <- function(SD_df, Coords_df, legend_pos = "right" ){
  #SD_df: output of the compute_SD_allSamples function
  #Coords_df: data frame containing the coordinates of each sample in the projection to use for the representation of the samples
  #legend_pos: optional argument to define the position of the legend
  SD_tmp = SD_df
  SD_tmp$Sample_ID = as.character( SD_tmp$Sample_ID)
  SD_tmp           = SD_tmp[order(SD_tmp$Sample_ID),]
  unique_sample_id = unique(SD_tmp$Sample_ID)
  SDmeans_ID       = sapply(unique_sample_id, function(x){
    mean(SD_tmp[which(SD_tmp$Sample_ID == x),3])
  })
  Coords_df <- Coords_df[names(SDmeans_ID),]
  SD_Coords_df <- cbind(  Coords_df, "SD" = SDmeans_ID)
  SD_Coords_df <- SD_Coords_df[order(SD_Coords_df$SD, decreasing = T),]
  colnames(SD_Coords_df)[1:2] = c("V1","V2")
  pSD_Main <- ggplot(  SD_Coords_df, aes(x=V1,y=V2,  color=SD/max(SD))) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
  pSD_Main <- pSD_Main +  labs(title="", 
                               x=colnames(Coords_df)[1], y=colnames(Coords_df)[2]) +
    theme( legend.position = legend_pos,
           plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
           plot.subtitle =element_text(size=14, hjust=0.5),
           plot.caption =element_text(size=12,  hjust=0.5),
           axis.title.x=element_text(size=16),  # X axis title
           axis.title.y=element_text(size=16),  # Y axis title
           axis.text.x=element_text(size=14),  # X axis text
           axis.text.y=element_text(size=14),
           legend.text = element_text(size = 10) ,
           legend.title = element_blank())
  return(list(SD_Coords_df,pSD_Main))
}

moran_I_knn <-function(expr_data , spatial_data, listK){
  #expr_data: matrix containing, for each sample (in rows), the values of the features (in columns) for which the MI values will be calculated
  #spatial_data: matrix containing the coordinates of each sample in the projection used to define the samples neighborhood and thus the weight matrix needed for the MI values computation
  #listK: vector listing the k values corresponding to the number of samples considered to define samples neighborhood
  MI_array <- array(NA, dim=c(length(listK), 2, ncol(expr_data)), 
                    dimnames=list(listK, c("obs","p.value"), colnames(expr_data)) )
  for(i in 1:length(listK)){
    KNN_R  = get.knn(spatial_data, k=listK[i], algorithm=c( "brute"))$nn.index
    m_neigh <- matrix(0, ncol = nrow(KNN_R), nrow =nrow(KNN_R))
    for (j in 1:nrow(KNN_R)) m_neigh[j,KNN_R[j,]] = 1
    MI_array[i,,] = apply(expr_data,2,function(co){res=Moran.I(co, m_neigh);     return(c(obs=res$observed,p.value=res$p.value))} )
  }
  return(MI_array)
}




