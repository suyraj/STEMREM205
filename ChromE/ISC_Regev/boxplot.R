reg_count_boxplot <- function(sc.list, order = NULL, split.by = "Developmental_potential", reg.cutoff = 0, frac = F){
  # Split list based on split.by criteria
  sc.list.split <- list()
  for(i in names(sc.list)){
    sc.list.split[[i]] <- SplitObject(sc.list[[i]], split.by = split.by)
  }
  if(frac == T){
    # Sum Regulons
    sum.list <- list()
    for(i in names(sc.list.split)){ # into each data set
      sum.list[[i]] <- list()
      for(j in names(sc.list.split[[i]])){ # into each developmental potential
        sum.list[[i]][[j]] <- c() # this array will hold the results
        for(cell in colnames(sc.list.split[[i]][[j]]@assays$RNA)){
          sum.list[[i]][[j]] <- c(sum.list[[i]][[j]],
                                  # Here is where we filter based on number
                                  sum(sc.list.split[[i]][[j]]@assays$RNA[,cell] > reg.cutoff) / length(sc.list.split[[i]][[j]]@assays$RNA[,cell]) )
        }
      }
    }
  } else{
    # Sum Regulons
    sum.list <- list()
    for(i in names(sc.list.split)){ # into each data set
      sum.list[[i]] <- list()
      for(j in names(sc.list.split[[i]])){ # into each developmental potential
        sum.list[[i]][[j]] <- c() # this array will hold the results
        for(cell in colnames(sc.list.split[[i]][[j]]@assays$RNA)){
          sum.list[[i]][[j]] <- c(sum.list[[i]][[j]], sum(sc.list.split[[i]][[j]]@assays$RNA[,cell] > reg.cutoff)) # Here is where we filter based on num
        }
      }
    }
  }
  # join tables
  table.join <- list()
  for(i in names(sum.list)){
    table.join[[i]] <- as.data.frame(t(do.call(cbind, sum.list[[i]])))
  }
  sc.table <- bind_rows(table.join) # the second bind to end with a single table
  # fix names of the table
  table.names <- c()
  for(i in names(table.join)){
    for(j in rownames(table.join[[i]])){
      table.names <- c(table.names, paste0(i,"_",j))
    }
  }
  rownames(sc.table) <- table.names # change the rownames to this now
  sc.table <- as.data.frame(t( sc.table )) # Transpose the Table
  # Order the table for plotting
  if(is.null(order) == FALSE){sc.table <- sc.table[ , order ] }# this orders sc.table in the correct order
  sc.stack <- stack(sc.table) # reshape table
  # Add colos
  sc.stack$group <- as.factor(sub("_[^_]+$", "", sc.stack$ind))  #add group to the stacked table
  # Generate Plot
  plot <- ggplot(sc.stack, aes(x = ind, y = values, fill = group)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Gene count") + ggtitle(paste0("Gene Count Cutoff = ", reg.cutoff))
  return(plot)
}