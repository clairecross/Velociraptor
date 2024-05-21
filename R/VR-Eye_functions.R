library(matrixStats)
library(FNN)
source("R/VR-Eye_support_functions.R")

parse_search_label <- function(ref.label) {
  split.ref.label <- str_split(ref.label, ' ')
  ref.marker.lst <- c() 
  ref.MEM.score <- as.data.frame(matrix(0, nrow=1, ncol=length(split.ref.label[[1]])))
  for (s in 1:length(split.ref.label[[1]])){ 
    split.score <- str_split(split.ref.label[[1]][s], '\\+')
    marker <- split.score[[1]][1] 
    score <- split.score[[1]][2] 
    ref.marker.lst <- c(ref.marker.lst, marker) 
    ref.MEM.score[1, s] <- as.numeric(score) 
  }
  colnames(ref.MEM.score) <- ref.marker.lst
  return(ref.MEM.score)
}

knn_MEM <- function(tsne.data,
                    transformed.data,
                    kvalue = 60,
                    verbose = FALSE) {
  
  ref.IQR <- median(colIQRs(as.matrix(transformed.data)))
  neighbor.index <- knnx.index(tsne.data, tsne.data, k=kvalue)
  all.labels <- as.data.frame(matrix(ncol=ncol(transformed.data), nrow=nrow(neighbor.index)))
  
  #add % complete readout to console 
  if (verbose) {
    for (r in 1:nrow(neighbor.index)) {
      cell.indices <- as.vector(neighbor.index[r, ])
      MEM.data <- transformed.data[cell.indices, ]
      MEM.data$cluster <- r 
      MEM.values <- MEM_edited(MEM.data, zero.ref=TRUE, input.IQR.ref=ref.IQR)
      all.labels[r, ] <- MEM.values[["prescaled_MEM_matrix"]][[1]]
      
      if (r %% 10 == 0) { 
        cat(paste0(round(r/nrow(neighbor.index), digits = 3), "% finished \n")) 
      } 
    }
  } else {
    for (r in 1:nrow(neighbor.index)) {
      cell.indices <- as.vector(neighbor.index[r, ])
      MEM.data <- transformed.data[cell.indices, ]
      MEM.data$cluster <- r 
      MEM.values <- MEM_edited(MEM.data, zero.ref=TRUE, input.IQR.ref=ref.IQR)
      all.labels[r, ] <- MEM.values[["prescaled_MEM_matrix"]][[1]]
    }
  }
  colnames(all.labels) <- colnames(MEM.values[["prescaled_MEM_matrix"]][[1]])
  
  #normalize to 0-10 scale
  scale_max <- max(abs(all.labels[,c(seq_len(ncol(all.labels)))]))
  MEM.scores <- (all.labels[,c(1:ncol(all.labels))]/scale_max)*10
  
  return(MEM.scores)
}

sim_scores <- function(common.MEM.scores,
                       common.ref.MEM.score) {
  #RMSD scores
  rmsd_vals <- as.data.frame(matrix(nrow=nrow(common.MEM.scores), ncol=1))
  colnames(rmsd_vals) <- 'RMSD score'
  
  #calculate RMSD for each cell compared to reference cell
  for (r in 1:nrow(common.MEM.scores)){ 
    sum_squares <- 0
    for (c in 1:ncol(common.MEM.scores)){ 
      sum_squares <- sum_squares + ((common.MEM.scores[r,c]-common.ref.MEM.score[1,c])^2)
    }
    rmsd_vals[r,] <- sqrt(sum_squares/ncol(common.MEM.scores))
  }
  
  #calculate similarity scores to the reference MEM label for each cell
  sim <- 100 - (rmsd_vals/10*100) 
  colnames(sim) <- "similarity_score"
  return(sim)
}

check_panel <- function(ref.label, 
                        panel.info) {
  ref.MEM.score <- parse_search_label(ref.label)
  common.ref.MEM.score <- as.data.frame(matrix(nrow = 1))
  common.markers <- c()
  common.ref.label <- c()
  
  #loop through each marker in searching label 
  for (marker in colnames(ref.MEM.score)) { 
    if (marker %in% panel.info$trimmed_names) {
      common.ref.MEM.score <- cbind(common.ref.MEM.score, ref.MEM.score[, marker])
      common.markers <- c(common.markers, marker)
      common.ref.label <- c(common.ref.label, paste0(marker, '+', ref.MEM.score[, marker]))
    } else {
      # print(paste0(marker,' is not measured in this dataset and will be excluded from similarity calculations'))
    }
  }
  #exclude first column 
  common.ref.MEM.score <- as.data.frame(common.ref.MEM.score[, 2:ncol(common.ref.MEM.score)])
  colnames(common.ref.MEM.score) <- common.markers
  return(list(common.ref.MEM.score, paste(common.ref.label, collapse=' ')))
}

marker.selection <- function(marker.list, 
                             token, 
                             position='before'){ #specify if marker name comes 'before' or 'after' the token
  #select markers to transform - include all marker channels with the token; exclude live/dead, fileID, tSNE/UMAP columns
  selected.indices <- c()
  selected.names <- c()
  for (i in 1:length(marker.list)){
    if (grepl(token, marker.list[i]) & !grepl('Event_length|Rh|Ir_|_ID|Ig_lambda|IgK|tSNE|tSNE|Dead|NA', marker.list[i])){ 
      selected.indices <- c(selected.indices, i)
      selected.names <- c(selected.names, marker.list[i])
    }
  }
  
  #exclude the channel from the marker name
  trimmed.names <- c()
  for (i in selected.names){
    marker <- str_split(i, token) #split the string by the token
    if (position =='before'){trimmed.marker <- str_split(marker[[1]][1], ' ')} #remove anything after the marker name
    if (position =='after'){trimmed.marker <- str_split(marker[[1]][2], ' ')} #remove anything before the marker name
    trimmed.names <- c(trimmed.names, trimmed.marker[[1]][1])
  }
  
  return(data.frame(indices=selected.indices, full_names=selected.names, trimmed_names=trimmed.names))
}