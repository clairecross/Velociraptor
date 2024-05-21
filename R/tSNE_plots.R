#plot cell density on t-SNE axes
tsne_density <- function(tsne.axes, export=TRUE, output_path){
  #plot
  tsne.by.density <- ggplot(tsne.axes, aes(x=tSNE2, y=tSNE1)) + coord_fixed() + 
    geom_point(size = 0.5) + geom_density_2d_filled(bins = 39) + 
    scale_fill_manual(values = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
                                 viridis::viridis(28,option = "A"))) +
    labs(x="t-SNE2", y="t-SNE1", title="Cell Density on t-SNE Axes", 
         caption="Data from Leelatian & Sinnaeve et al., eLife. 2020") + theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "none") 
  #export
  if(export==TRUE){
    png(paste(output_path, strftime(Sys.time(),"%Y-%m-%d"), " t-SNE density.png"), height = 1000, width = 1000, res=200)
    print(tsne.by.density)
    dev.off()
  }
}

#plot heat on markers
tsne_heat_on_markers <- function(tsne.axes, transformed.marker.data, n_col, export=TRUE, output_path){
  #calculate number of rows
  n_row <- ceiling(ncol(transformed.marker.data)/n_col)
  #plot
  tsne.by.marker <- as_tibble(tsne.axes) %>%
    bind_cols(transformed.marker.data)  %>%
    gather(channel, intensity, -tSNE1, -tSNE2) %>%
    mutate(across(channel,factor))%>%
    group_split(channel) %>%
    map(
      ~ggplot(.,aes(x=tSNE2, y=tSNE1, col=intensity)) +
        geom_point(size = 3) +
        scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(5)) +
        facet_grid(~ channel, labeller = function(x) label_value(x, multi_line = FALSE)) +
        coord_fixed() +
        theme_bw()+
        theme(strip.text.x=element_text(size=20), legend.title=element_blank(), panel.grid=element_blank()))%>%
    plot_grid(plotlist = ., align = 'hv', ncol = n_col)
  #export
  if(export==TRUE){
    png(paste(output_path, strftime(Sys.time(),"%Y-%m-%d"), " t-SNE_heat on markers.png"),height = n_row*500, width = n_col*500)
    print(tsne.by.marker)
    dev.off()
  }
}
