
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


############################################################################################
Seq_graph_by_k  <-function (data_Seq, Names=NULL, data_diff_mean_K = NULL, log = FALSE){
  if (is.null(data_diff_mean_K) == TRUE) {
    data_diff_mean_k <- data.frame("k" =  unique(data_Seq$K))
    for (j in seq(from = 3, to = dim(data_Seq)[2], by = 1)) {
      mean_by_k <- tapply(data_Seq[, j], data_Seq$K , mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, mean_by_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_Seq)[3:dim(data_Seq)[2]]
    if (is.null(Names) == FALSE){
      if (length(Names) != (dim(data_Seq)[2] - 3)){
        warning("The list of names gave as input doesn't match with the number of curve.")
      }
      else{
        colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
      }
    }
  }
  else{
    data_diff_mean_k <- data_diff_mean_K
  }
  
  if (log == FALSE){
    data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' =( data_diff_mean_k[, 2]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
    if (dim(data_diff_mean_k)[2]>=3){
      for (i in 3:(dim(data_diff_mean_k)[2])){
        c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = (data_diff_mean_k[, i]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
        data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
      }
    }
  
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="Sequence difference metric", y=TeX("mean(SD)_k"), x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                   plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                   plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                   axis.title.x=element_text(size=12, face="italic"),  # X axis title
                                                   axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                   axis.text.x=element_text(size=12),  # X axis text
                                                   axis.text.y=element_text(size=12),
                                                   legend.title = element_blank())  # Y axis text
  print(p)
  return(p)
  }
  else {
    data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = log( data_diff_mean_k[, 2]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
    if (dim(data_diff_mean_k)[2]>=3){
      for (i in 3:(dim(data_diff_mean_k)[2])){
        c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = log(data_diff_mean_k[, i]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
        data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
      }
    }
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="Sequence difference metric",   y=TeX("log(mean(SD)_k)"), x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                     plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                     plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                     axis.title.x=element_text(size=12, face="italic"),  # X axis title
                                                     axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                     axis.text.x=element_text(size=12),  # X axis text
                                                     axis.text.y=element_text(size=12),
                                                     legend.title = element_blank())  # Y axis text
    print(p)
    return(p)
  } 
}
seq_permutation_test <- function(data, data_ref, list_K, n=30, graph = TRUE){
  
  if (n > 30){
    warning("the calcul could be long !")
  }
  colnames(data)[1] <- 'Sample_ID' ; colnames(data_ref)[1] <- 'Sample_ID'
  if (dim(data)[1] != dim(data_ref)[1]){
    warning("Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
  }
  else if( dim(data)[1] == dim(data_ref)[1] & sum(as.character(data[, 1]) == as.character(data_ref[, 1])) != length(data_ref[, 1])){
    warning("Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
  }
  data_m <- merge(data, data_ref, by=c('Sample_ID'))
  data <- data_m[, 1:dim(data)[2]]
  data_ref <- data_m[, (dim(data)[2]+1):dim(data_m)[2]]
  data_ref <- cbind(data_m[, 1], data_ref)
  colnames(data_ref)[1] <- "Sample_ID"
  global_seq_df <- Seq_calcul(list(data), dataRef = data_ref , listK = list_K)[[1]]
  mean_k <- tapply(global_seq_df$Seq, global_seq_df$K, mean)
  main_df <- data.frame('k' = unique(global_seq_df$K) , "means_ref" = mean_k)
  for (i in 1:n){
    data_shuffle <- data[,2:dim(data)[2]]
    data_shuffle <- data_shuffle[,sample(ncol(data_shuffle))]
    data_shuffle <- data_shuffle[sample(nrow(data_shuffle)),]
    data_shuffle <- cbind(data[,1], data_shuffle, row.names = NULL)
    colnames(data_shuffle)[1] <- "Sample_ID"
    Seq_data_A <- Seq_calcul(list(data_shuffle), dataRef = data_ref , listK = list_K)[[1]]
    mean_k <- tapply(Seq_data_A$Seq, Seq_data_A$K, mean)
    main_df <- cbind(main_df , mean_k)
  }
  by_k_alea <- main_df[,3:dim(main_df)[2]]
  Means_alea <- rowMeans(by_k_alea)
  WT  = wilcox.test(main_df[ ,1], Means_alea, alternative = "less")
  #print(WT)
  
  theme_set(theme_bw())
  p <- ggplot()
  for (i in 3:(dim(main_df)[2])){
    c_df <- data.frame('k' = main_df[ ,1] , 'main_df' = main_df[ ,i])
    p <- p + geom_line(data = c_df, aes(x=k, y=main_df), colour = '#848484')+geom_point(data = c_df, aes(x=k, y= main_df), colour = '#848484')
  }
  c_df <- data.frame('k' = main_df[ ,1] , 'main_df' = main_df[ ,2])
  p <- p + geom_line(data = c_df, aes(x=k, y = main_df), colour = '#B40404')+geom_point(data = c_df, aes(x=k, y=main_df), colour = '#B40404')
  
  c_MA_df <- data.frame('k' = main_df[ ,1] , 'main_df' = Means_alea)
  p <- p + geom_line(data = c_MA_df, aes(x=k, y = main_df), colour = '#388E3C')+geom_point(data = c_MA_df, aes(x=k, y=main_df), colour ='#388E3C')
  
  p <- p +  labs(title="Significance test of the Sequence difference metric",
                 y=TeX("mean(SD)_k"), x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                   plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                   plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                   axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                   axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                   axis.text.x=element_text(size=12),  # X axis text
                                                   axis.text.y=element_text(size=12))  # Y axis text
  print(p)
  

  return(list(WT,p))
}





