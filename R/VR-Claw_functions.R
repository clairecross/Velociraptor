#calculate the percentage of each patient's cells that are in a defined neighborhood
calc_nbh_percent <- function(pl, cp_df_r){ #pass in patient_list and a row of cluster_and_patient_df
  cp_df_r <- as.data.frame(t(cp_df_r)) #reset data type to row in a dataframe
  #for each patient, what percentage of their cells are in this neighborhood?
  nbh_percent <- length(cp_df_r[1:kval][,which(cp_df_r[1:kval] %in% pl)])*100/sample_num
  return(nbh_percent)
}

#determine how many patients are present in a defined neighborhood
calc_num_patients <- function(cp_df){
  num.patients <- length(unique(wafflecut(cp_df, intervals=cell.index.intervals)))
  return(num.patients)
}

#determine if a defined neighborhood is associated with the selected clinical columns
nbh_analysis <- function(cp_df){
  #what is each patient's cell abundance for this particular neighborhood
  nbh_percent <- sapply(patient_list, calc_nbh_percent, cp_df)
  
  #what is that neighborhood's IQR (75th percentile - 25th percentile)
  abundance.IQR <- IQR(nbh_percent) 
  
  #divide patients into high and low groups based on subset abundance and test for enrichment
  abundance.groups <- ifelse(nbh_percent > abundance.IQR, 1, 0) # 1 = high, 0 = low
  
  #Cox proportional hazards model
  Group <- factor(abundance.groups, levels = c(0,1), labels = c("Low", "High"))
  group.survival.data <- cbind(OS.data, Group)
  # high.low.groups <- split(group.survival.data, group.survival.data$Group)
  # low.median <- median(high.low.groups[["Low"]][[clinical_col]])
  # high.median <- median(high.low.groups[["High"]][[clinical_col]])
  coxph.model <- coxph(Surv(OS.data[,clinical_col], OS.data[,status_col]) ~ Group,
                       data=group.survival.data)
  return(summary(coxph.model))
}

#dbscan clustering
dbscan_fxn <- function(df, embedding_idxs, sig_level, epsilon, minPoints){
  #select embedding columns to cluster on
  embedding <- df[,embedding_idxs]
  
  #cluster
  dbscan.results <- dbscan::dbscan(embedding, eps = epsilon, minPts = minPoints)
  df$cluster <- dbscan.results$cluster
  
  #remove noise cluster 0
  df <- df[which(df$cluster!=0),]
  
  #rename clusters to include significance stats
  df$cluster <- paste0(df$cluster, sig_level)
  df$cluster <-as.numeric(df$cluster)
  
  return(df)
}

#cluster_survival_analysis
patient_cluster_abundance <- function(patient.list){
  #calculate the percentage of each patient's cells in each cluster 
  subset.abundance <- list()
  for (p in 1:length(patient.list)){ #looping over p = patients
    subset.abundance[[p]] = (summary(patient.list[[p]][,'cluster']))*100/nrow(patient.list[[p]]) #cell abundance in each cluster
  }
  all.subset.abundances = ldply(subset.abundance, rbind)  #dataframe of patient cell abundance in each cluster
  all.subset.abundances[is.na(all.subset.abundances)]<-0 
  
  return(all.subset.abundances)
}

#Cox proportional hazards on clusters
cox_ph_cluster_fxn <- function(all.subset.abundances){
  #split patients into HIGH/LOW for each cluster by comparing to the cluster's abundance IQR value
  abundance.groups <- list()
  for(c in 1:ncol(all.subset.abundances)){ # looping over c = cluster
    abundance.IQR <- IQR(all.subset.abundances[,c]) #cluster's abundance IQR value
    abundance.groups[[c]] <- (all.subset.abundances[,c]>abundance.IQR) #TRUE=abundance greater than the IQR
    abundance.groups[[c]][abundance.groups[[c]]==TRUE] = 1 # 1 = high, 0 = low
  }   
  #Cox proportional hazards model
  cox.cluster.summary <- list()
  for (s in 1:ncol(all.subset.abundances)){ #looping over s = cluster
    Group <- factor(abundance.groups[[s]], levels = c(0,1), labels = c("Low", "High"))
    group.survival.data <- cbind(OS.data, Group)
    coxph.cluster.model <- coxph(Surv(OS.data[,clinical_col], OS.data[,status_col]) ~ Group, data=group.survival.data)
    cox.cluster.summary[[colnames(all.subset.abundances)[s]]] <- summary(coxph.cluster.model)
  }
  
  #plotting info
  survival.plots <- list()
  final.prog.clusters <- c()
  pval.list <- c()
  HR.list <- c()
  for (s in 2:ncol(all.subset.abundances)){ #start at 2 to exclude non VR clusters
    if (colnames(all.subset.abundances)[s]!="(Other)"){
      cox.cluster.coefficients <- cox.cluster.summary[[s]]$coefficients
      CI <- cox.cluster.summary[[s]]$conf.int
      Group <- factor(abundance.groups[[s]], levels = c(0,1), labels = c("Low", "High"))
      group.survival.data <- cbind(OS.data, Group)
      coxph.model <- survfit(Surv(OS.data[,clinical_col], OS.data[,status_col]) ~ Group, data=group.survival.data)
      #display cluster stats
      cluster.title <- paste0("Subset #", colnames(all.subset.abundances)[s], " (p = ", 
                              round(cox.cluster.coefficients[,5],3), ", HR = ",
                              round(cox.cluster.coefficients[,2],3), ", CI[",
                              round(CI[,3],3),",",round(CI[,4],3),"])")
      pval.list <- c(pval.list, cox.cluster.summary[[s]][["coefficients"]][,5])
      HR.list <- c(HR.list, cox.cluster.summary[[s]][["coefficients"]][,2])
      #plot significant clusters
      if (cox.cluster.coefficients[,5] <= pval.cutoff){
        final.prog.clusters <- c(final.prog.clusters, colnames(all.subset.abundances)[s])
      }
      survival.plots[[s-1]] <- ggsurvplot(coxph.model, 
                                        data=group.survival.data, 
                                        conf.int=F, 
                                        pval=F, 
                                        risk.table=T, 
                                        tables.y.text = FALSE, 
                                        legend.labs = c("Low", "High"), 
                                        legend.title = "Group",
                                        censor.shape=124,
                                        title = cluster.title)
    }
  }
  return(list(survival.plots, pval.list, HR.list))
}
