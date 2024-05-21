theme_VRPTR <- function() {
  theme_bw() %+replace%
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(size = 16)
    )
}

VR_Eye_sim_plot <- function(all.data, ref.label, min.sim=-1) {
  plot.data <- all.data[, c("tSNE1","tSNE2","similarity_score")]
  plot.data <- plot.data[order(plot.data$similarity_score), ]
  
  #calculate minimum similarity score and color palette ramp 
  if (min.sim==-1){min.sim <- min(plot.data$similarity_score)}
  sim.step <- (100-min.sim)/4
  
  #calculate a ratio to make the t-SNE axes a square
  range <- apply(apply(plot.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  sim.plot <- ggplot(plot.data) +
    geom_point(aes(x = tSNE2, y = tSNE1, col = similarity_score)) +
    coord_fixed(ratio = graphical.ratio) + 
    scale_color_gradientn(colors = c("lightgrey", "lightgrey", "#5E4FA2", "#88CFA4", "#FFFFBF", "#F88D52", "#9E0142"), 
                          values = rescale(c(0, min.sim-0.00001, min.sim, min.sim+sim.step, min.sim+2*sim.step, min.sim+3*sim.step, 100)), 
                          name = "Similarity", limits=c(0,100)) +
    labs(
      title = paste0("Similarity: ", ref.label), 
      x = "t-SNE2", y = "t-SNE1", 
      caption="Data from Leelatian & Sinnaeve et al., eLife. 2020"
    ) +
    theme_VRPTR()
  
  return(sim.plot)
}

VR_Eye_bin_plot <- function(all.data, ref.label) {
  plot.data <- all.data[order(all.data$sim_bin), c("tSNE1","tSNE2","sim_bin")]
  plot.data$sim_bin <- as.factor(plot.data$sim_bin)
  bin.colors <- colorRampPalette(c('lightgray','lightgreen','darkgreen'))(length(my.bins) + 1)
  names(bin.colors) <- levels(all.data$sim_bin)
  
  #calculate a ratio to make the t-SNE axes a square
  range <- apply(apply(plot.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  bin.plot <- ggplot(plot.data) +
    geom_point(aes(x = tSNE2, y = tSNE1, col = sim_bin)) +
    coord_fixed(ratio = graphical.ratio) + 
    scale_color_manual(values = bin.colors, name = ("Similarity bin")) +
    labs(
      title = paste0("Similarity binned: ", ref.label), 
      x = "t-SNE 2", y = "t-SNE 1"
    ) +
    theme_bw() +
    theme_VRPTR()
  
  return(bin.plot)
}

VR_Eye_cluster_plot <- function(all.data,
                           cluster.data, 
                           sim.plot = NULL,
                           color.scheme = "tatarize", 
                           ref.label = "") {
  
  plot.data <- all.data[, c("tSNE1","tSNE2")]
  
  num_clusters = length(unique(cluster.data$cluster))
  
  if (color.scheme == "qual") {
    set.seed(5)
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cluster.colors <- sample(col_vector)
  } else if (color.scheme == "sim") {
    cluster.colors <- sim_colors(num_clusters, sim.plot)
  } else if (color.scheme == "tatarize") { 
    cluster.colors <- tatarize_optimized(num_clusters)
  }

  #calculate a ratio to make the t-SNE axes a square
  range <- apply(apply(plot.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  cluster.plot <- ggplot(plot.data) +
    geom_point(aes(x = tSNE2, y = tSNE1), col = "lightgray") +
    geom_point(data = cluster.data, aes(x = tSNE2, y = tSNE1, col = as.factor(cluster))) +
    coord_fixed(ratio = graphical.ratio) + 
    scale_color_manual(values = cluster.colors, name = "Cluster", labels = cluster.avg$label) +
    labs(
      title = paste0("Ordered clusters: ", ref.label), 
      x = "t-SNE 2", y = "t-SNE 1"
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme_VRPTR()
  
  return(cluster.plot)
}

sim_colors <- function(num.clusters, sim.plot) {
  sim.built = ggplot_build(sim.plot)
  sim.colors = unique(sim.built$data[[1]]["colour"])[, 1]
  
  if (num.clusters > length(sim.colors)) {
    cluster.colors <- rev(colorRampPalette(sim.colors)(num.clusters))
  } else {
    set.seed(1)
    color.indxs = sample(length(sim.colors), num.clusters)
    color.indxs <- color.indxs[order(color.indxs, decreasing = TRUE)]
    cluster.colors <- sim.colors[color.indxs]
  }
  return(cluster.colors)
}