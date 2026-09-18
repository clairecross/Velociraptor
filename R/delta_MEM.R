parse_label <- function(input_str) {
  
  # if there is a tab
  if (any(grepl("\t",  input_str))) { 
    input_str <- unlist(strsplit(input_str, "\t"))[2]
  }
  input_split = unlist(strsplit(input_str, " "))
  label_names = c()
  label_values = c()
  for (term in input_split) {
    
    # handle cluster number 
    if ((term == input_split[1] | term == input_split[2]) & suppressWarnings(!is.na(as.numeric(term)))) {
      label_names <- c("cluster")
      label_values <- as.numeric(term)
    } else if (term != ":" & term != "▲" & term != "None" & term != "▼") {
      term_split = unlist(strsplit(term, "\\+"))
      label_names[length(label_names) + 1] <- term_split[1]
      label_values[length(label_values) + 1] <- as.numeric(term_split[2])
    }
  }
  names(label_values) <- label_names
  return(label_values)
}

parse_file <- function(filepath) {
  label_set = read.csv(filepath)
  parsed = list()
  for (i in 1:nrow(label_set)) {
    parsed[[i]] <- parse_label(label_set[i, ])
    names(parsed)[i] <- parsed[[i]][1]
  }
  return(parsed)
}

compare_labels <- function(label1, label2) {
  delta.label = NULL
  delta.table = NULL
  
  # markers in both labels
  for (marker in intersect(names(label1), names(label2))) {
    if (marker == "cluster") {
      if (label1[marker] == label2[marker]) { 
        delta.label <- paste0(label1[marker]," \U1D6AB : ")
      } else { 
        delta.label <- paste0(label1[marker]," \U1D6AB ",label2[marker]," : ")
      }
    } else {
      num = label2[marker] - label1[marker]
      if (num != 0) { 
        delta.table[length(delta.table) + 1] <- num
        names(delta.table)[length(delta.table)] <- marker
      } 
    }
  }  
  
  # markers in label2, but not label1 
  for (marker in setdiff(names(label2), names(label1))) {
    if (marker == "cluster") { 
      warning("A cluster number was only found in one comparison label. Printing delta MEM label without cluster number. \n")
    } else {
      delta.table[length(delta.table) + 1] <- label2[marker]
      names(delta.table)[length(delta.table)] <- marker
    }
  }
  
  # markers in label1, but not label2
  for (marker in setdiff(names(label1), names(label2))) {
    if (marker == "cluster") { 
      warning("A cluster number was only found in one comparison label. Printing delta MEM label without cluster number. \n")
    } else {
      delta.table[length(delta.table) + 1] <- 0 - label1[marker]
      names(delta.table)[length(delta.table)] <- marker  
    }
  }
  
  # format as string
  if (is.null(delta.table)) {
    delta.label <- paste0(delta.label, "No change")
  } else {
    delta.table <- delta.table[order(abs(delta.table), decreasing = TRUE)]
    for (i in 1:length(delta.table)) {
      marker = names(delta.table)[i]
      num = delta.table[i]
      if (num >= 0) {
        num <- paste0("+", num)
      }
      delta.label <- paste0(delta.label, marker, num, " ")
    }
    delta.label <- substr(delta.label, 1, nchar(delta.label) - 1)
  }
  return(delta.label)
}  

delta_MEM <- function(X, Y = NULL) {
  all.parsed = list()
  delta.labels = vector()
  
  # detect input types
  if (dir.exists(X)) {
    filelist = list.files(path = X, pattern = ".txt")
    if (length(filelist) == 2) {
      Y <- paste0(X, "/", filelist[[2]])
      X <- paste0(X, "/", filelist[[1]])
    } 
  } 
  if (is.null(Y)) {
    stop("One of your files could not be found. Please check your folder or file paths.")
  }
  
  if (file.exists(X)) {
    if (!file.exists(Y)) {
      stop("One of your files could not be found. Please check your folder or file paths.")
    } else {
      all.parsed[[1]] <- parse_file(X)
      all.parsed[[2]] <- parse_file(Y)
      
      # get labels for matching clusters
      common.clusters = intersect(names(all.parsed[[1]]), names(all.parsed[[2]]))
      for (clust in common.clusters) {
        label1 = all.parsed[[1]][[clust]]
        label2 = all.parsed[[2]][[clust]]
        delta.labels[length(delta.labels) + 1] <- compare_labels(label1, label2)
      }
      return(delta.labels)
    }
  } else if (!file.exists(X)) {
    if (file.exists(Y)) {
      stop("One of your files could not be found. Please check your folder or file paths.")
    } else if (is.character(X) & is.character(Y)) {
      return(compare_labels(parse_label(X), parse_label(Y)))
    } 
  }
}