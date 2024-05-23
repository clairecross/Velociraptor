#VR-Claw plot 
VR_Claw_plot_fxn <- function(df=data.for.KNN, #data used for KNN
                             x_idx=1, 
                             y_idx=2,
                             survival_stats=cox.summary, 
                             plot.title='VR-Claw', 
                             xlab='dim1', 
                             ylab='dim2'){
  #set significance status for each cell
  plotting.data <- data.frame(dim1=df[,1], 
                              dim2=df[,2], 
                              cluster=c(1:nrow(df)), 
                              num.neighborhood.patients=num.patients, 
                              status="p>0.1")

  for (z in 1:length(survival_stats)){
    cox.coef <- survival_stats[[z]]$coefficients
    #each neighborhood must have the minimum number of patients per cluster to be significant
    if (plotting.data$num.neighborhood.patients[z] >= min.patients.per.cluster){
      if (cox.coef[,c(5)] <= 0.01 & cox.coef[,c(2)] <= 1){ #p<0.01, HR<1
        plotting.data$status[z] <- "p<0.01 HR<1"
      }else if (cox.coef[,c(5)] <= 0.05 & cox.coef[,c(2)] <= 1){ #p<0.05, HR<1
        plotting.data$status[z] <- "p<0.05 HR<1"
      }else if (cox.coef[,c(5)] <= 0.1 & cox.coef[,c(2)] <= 1){ #p<0.1 HR<1
        plotting.data$status[z] <- "p<0.1 HR<1"
      }else if (cox.coef[,c(5)] <= 0.01 & cox.coef[,c(2)] >= 1){ #p<0.01 HR>1
        plotting.data$status[z] <- "p<0.01 HR>1"
      }else if (cox.coef[,c(5)] <= 0.05 & cox.coef[,c(2)] >= 1){ #p<0.05 HR>1
        plotting.data$status[z] <- "p<0.05 HR>1"
      }else if (cox.coef[,c(5)] <= 0.1 & cox.coef[,c(2)] >= 1){ #p<0.1 HR>1
        plotting.data$status[z] <- "p<0.1 HR>1"
      }
    }
  }

  #get maximum x and y values to set the scales
  max_x = max(c(max(plotting.data[,x_idx]),abs(min(plotting.data[,x_idx]))))
  max_y = max(c(max(plotting.data[,y_idx]),abs(min(plotting.data[,y_idx]))))
      
  #reorder data so the significant cells are displayed on top
  light.grey.cells <- plotting.data[which(plotting.data$status=="p>0.1"),]
  dark.grey.cells <- plotting.data[which(plotting.data$status=="p<0.1 HR>1" | plotting.data$status=="p<0.1 HR<1"),]
  light.red.cells <- plotting.data[which(plotting.data$status=="p<0.05 HR>1"),]
  dark.red.cells <- plotting.data[which(plotting.data$status=="p<0.01 HR>1"),]
  light.blue.cells <- plotting.data[which(plotting.data$status=="p<0.05 HR<1"),]
  dark.blue.cells <- plotting.data[which(plotting.data$status=="p<0.01 HR<1"),]
  to.plot <- rbind(light.grey.cells, dark.grey.cells, 
                   light.red.cells, light.blue.cells,
                   dark.red.cells, dark.blue.cells)
    
  #plot
  Claw_plot <- ggplot(to.plot) + theme_bw() + coord_fixed(ratio = 1) +
                geom_point(aes(x=to.plot[,x_idx], y=to.plot[,y_idx], 
                               col = as.factor(to.plot$status))) + 
                scale_color_manual(values = col.palette) + 
                guides(colour = guide_legend(override.aes = list(size=5))) + 
                theme(panel.grid=element_blank()) + 
                labs(color = "Prognostic Subsets", title=plot.title, x=xlab, y=ylab) +
                ylim(-max_y,max_y) + xlim(-max_x,max_x)
  
  return(Claw_plot)
}

#population survival plots
