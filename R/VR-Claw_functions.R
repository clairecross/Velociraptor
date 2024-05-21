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