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
  if(is.data.frame(survival_stats)){plotting.data$status <- survival_stats$status}else{
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
pop_survival_analysis <- function(population, pop.name, col.pal){
  group.names <- c(paste0(pop.name, "-Low"), paste0(pop.name, "-High"))
  
  #calculate pop abundance per
  population$sum <- rowSums(population)
  pop.abundance.IQR <- IQR(population$sum)
  population$group <- ifelse(population$sum>=pop.abundance.IQR, 1, 0)
  Group <- factor(population$group, levels = c(0, 1), labels = group.names)
  factored.data <- cbind(OS.data, Group) #select OS Time, OS Status and Group
  
  #calculate cox proportional hazards model 
  model.to.plot <- survfit(Surv(OS.data[,clinical_col], OS.data[,status_col]) ~ Group, data=factored.data) #graphing data
  model.for.stats <- coxph(Surv(OS.data[,clinical_col], OS.data[,status_col]) ~ Group, data=factored.data) #Cox PH Model
  model.stats <- summary(model.for.stats) #summarize cox ph
  
  #access interesting stats from Cox PH Model
  pval <- round(model.stats[["coefficients"]][,'Pr(>|z|)'], 3)
  HR <- round(model.stats[["coefficients"]][,'exp(coef)'], 3)
  conf.ints <- c(round(model.stats[["conf.int"]][,3], 2), 
                 round(model.stats[["conf.int"]][,4], 2))
  CI <- paste0('[', min(conf.ints), ',', max(conf.ints), ']')
  
  #plot
  pop.plot <- ggsurvplot(model.to.plot, data=factored.data, risk.table=T,
                         legend.labs = group.names, legend.title = "Group",
                         censor.shape=124, palette=col.pal)
  pop.plot$plot <- pop.plot$plot + annotate("text", x=Inf, y=Inf, vjust=1, hjust=1, size=5,
                                            label=paste0('p-value = ', pval, '
                                                         HR = ', HR, '
                                                         95% CI ', CI)) #annotate with stats
  return(list(population, pop.plot))
}
