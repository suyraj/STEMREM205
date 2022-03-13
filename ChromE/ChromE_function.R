# ChromE: Chromosome Expression Maps
# S205 Group Project



#' ChromE
#'
#' @param seurat_object Seurat object
#' @param chr_table modified gene table in order they are found in the chromosome
#' @param annotation the meta data column that contains the celltype
#' @param order here you can put an order for the annotation that you want to see stuff in
#' @param assay the assay that you want averaged. Usually this is the RNA assay
#' @param k This is how many genes the sliding windows should use. default is 200
#' @param plot.range this dictates the range of the heatmap. defaults to using the data but you can set this range yourself
#'
#'
#' @return Returns a Complex Heatmap
#' @export
#'
#' 
ChromE <- function(seurat_object, chr_table, annotation, order = NULL, assay = "RNA", k = 250, plot.range = c(0,50)){
  
  ### Check to make sure required packages installed ###
  if(require("Seurat") == F) stop("Package Error: Please install & load Seurat") # Checks for Seurat install
  if(require("zoo") == F) stop("Package Error: Please install & load zoom") # Checks for zoom install
  if(require("ComplexHeatmap") == F) stop("Package Error: Please install / load ComplexHeatmap") # Checks for ComplexHeatmap install
  
  ### Check for errors ###
  if(class(seurat_object) != "Seurat") stop("Error: the data is not given as a Seurat Object") # checks to make sure you're giving it a seurat object
  if("expression"  %in% colnames(chr_table) == F) stop("Error: The reference table does not have a column called expression") # checks to make sure there is an expression column
  if(annotation  %in% colnames(seurat_object@meta.data) == F)  stop("Error: the annotation you provided does not exist in the seurat metadata") # checks to make sure annotation is present in Seurat metadata
  if(is.null(order) == F  & !all(order %in% seurat_object@meta.data[,annotation])) stop("Error: the order you've given contains items that do not exist in the annotation column. Check the levels() of your annotation") # checks to make sure order is present in the annotation column
  
  ### Generate Table for Plotting ###
  # Generate a table with the average expression
  avg_expression <- as.data.frame(AverageExpression(seurat_object, group.by = annotation, slot = "scale.data")[[assay]])
  
  
  # Add the gene expression from the average expression table in the order of the chromosome###
  chrom_expression_table <- list()
  
  for(cell_type in colnames(avg_expression)){
    chrom_expression_table[[cell_type]] <-  data.frame(gene = chr_table$gene,chr = chr_table$chr, order = chr_table$order, expression = chr_table$expression) # move genes from chr to list
    rownames(chrom_expression_table[[cell_type]]) <- chr_table$gene
  
    for(gene in rownames(avg_expression)){
      chrom_expression_table[[cell_type]][[gene, "expression"]] <-  avg_expression[[gene, cell_type]] # this moves the expression over
    }
  }
  
  
  ### Sliding Window to better see the expression ###
  for(cell_type in names(chrom_expression_table)){
    chrom_expression_table[[cell_type]]$expression <- rollmax(chrom_expression_table[[cell_type]]$expression,
                                                              k = k,
                                                              fill = 0,
                                                              partial = T)
  }
  
  ## Setting range of the heatmap
  if(is.null(plot.range) == FALSE){
    col.range <- plot.range
  }else{
    all.range <- c()
    for(cell_type in names(chrom_expression_table)){
      all.range <- c(all.range, range(chrom_expression_table[[cell_type]]$expression))
    }
    col.range <- range(all.range) # Get the range of numbers from the table
  }
  
  # Catch if input plot range is wrong
  if(length(col.range) != 2 & class(col.range) != "numeric") stop("Error: the plot range you've given might be in the wrong format. Input should look like: c(min, max)")
  
  ### Split each cell type based on chromosome ###
  for(cell_type in names(chrom_expression_table)){
    chrom_expression_table[[cell_type]] <- split(chrom_expression_table[[cell_type]], chrom_expression_table[[cell_type]]$chr) # split based on the chromosome
    
    for(c in names(chrom_expression_table[[cell_type]])){
      chrom_expression_table[[cell_type]][[c]] <- arrange(chrom_expression_table[[cell_type]][[c]], order) # re-order just in case
    }
    
  }
  
  ### Optionally order the table if there's a specific order you want
  if(is.null(order) == FALSE){
    chrom_expression_table <- chrom_expression_table[ order ]
  }
  

  ### Combine Based on chromosomes ###
  chr_expr_table <- list()
  
  for(c in names(chrom_expression_table[[1]])){
    chr_expr_table[[c]] <- data.frame(row.names = rownames(chrom_expression_table[[1]][[c]])) # this creates the dataframe
    
    for(cell_type in names(chrom_expression_table)){
      chr_expr_table[[c]][cell_type] <- chrom_expression_table[[cell_type]][[c]]$expression
    }
    
  }  
  
  ### Cleanup Table for Plotting ###
  for(c in names(chr_expr_table)){
    chr_expr_table[[c]] <- t( chr_expr_table[[c]] )
  }
  
 
  
  ### Plotting Table ###
  
  col_fun <- colorRamp2(seq(col.range[[1]], col.range[[2]], length.out = 50), viridis(50)) # make colors for the range of numbers
  
  # Plot all the individual chromosome heatmaps
  chr_plots <- NULL # create an empty heatmap list object
  for(c in names(chr_expr_table)){
    
    # only keep last legend
    keep.legend = FALSE
    if(c == end(names(chr_expr_table))[[1]]){keep.legend = TRUE}
    chr_title <- c
    
    # setting x y & m
    if(length(chr_expr_table) == 25){ #25 would be human dataset!
      if(c == 23){chr_title <- "X"}else if(c == 24){chr_title <- "Y"} else if(c == 25){chr_title <- "M"}
    } 
    
    if(length(chr_expr_table) == 23){ #23 would be mouse dataset!
      if(c == 21){chr_title <- "X"}else if(c == 22){chr_title <- "Y"} else if(c == 23){chr_title <- "M"}
    } 
    
    # generate heatmap
    chr_plots = chr_plots + Heatmap(chr_expr_table[[c]], 
                                    col = col_fun, 
                                    cluster_rows = F, cluster_columns = F,
                                    column_title = chr_title,
                                    show_column_names = F, name = " ",
                                    show_heatmap_legend = keep.legend,
                                    border_gp = gpar(col= "black", lty = 1, lex = 2, lineend = "round", linejoin = "round"))
  }
  
  
  # Now combine all the individual heatmaps together
  plt <- draw(chr_plots, ht_gap = unit(5, "points"))
  
  
  return(plt)  
}


#' ChromE_table
#' Generates just the ChromE output table
#'
#' @param seurat_object 
#' @param chr_table 
#' @param annotation 
#' @param order 
#' @param assay 
#'
#' @return
#' @export
#'
#' @examples
ChromE_table <- function(seurat_object, chr_table, annotation, order = NULL, assay = "RNA"){
  
  ### Check to make sure required packages installed ###
  if(require("Seurat") == F) stop("Package Error: Please install & load Seurat") # Checks for Seurat install
  if(require("zoo") == F) stop("Package Error: Please install & load zoom") # Checks for zoom install
  if(require("ComplexHeatmap") == F) stop("Package Error: Please install / load ComplexHeatmap") # Checks for ComplexHeatmap install
  
  ### Check for errors ###
  if(class(seurat_object) != "Seurat") stop("Error: the data is not given as a Seurat Object") # checks to make sure you're giving it a seurat object
  if("expression"  %in% colnames(chr_table) == F) stop("Error: The reference table does not have a column called expression") # checks to make sure there is an expression column
  if(annotation  %in% colnames(seurat_object@meta.data) == F)  stop("Error: the annotation you provided does not exist in the seurat metadata") # checks to make sure annotation is present in Seurat metadata
  if(is.null(order) == F  & !all(order %in% seurat_object@meta.data[,annotation])) stop("Error: the order you've given contains items that do not exist in the annotation column. Check the levels() of your annotation") # checks to make sure order is present in the annotation column
  
  ### Generate Table for Plotting ###
  # Generate a table with the average expression
  avg_expression <- as.data.frame(AverageExpression(seurat_object, group.by = annotation, slot = "data")[[assay]])
  
  
  # Add the gene expression from the average expression table in the order of the chromosome###
  chrom_expression_table <- list()
  
  for(cell_type in colnames(avg_expression)){
    chrom_expression_table[[cell_type]] <-  data.frame(gene = chr_table$gene,chr = chr_table$chr, order = chr_table$order, expression = chr_table$expression) # move genes from chr to list
    rownames(chrom_expression_table[[cell_type]]) <- chr_table$gene
    
    for(gene in rownames(avg_expression)){
      chrom_expression_table[[cell_type]][[gene, "expression"]] <-  avg_expression[[gene, cell_type]] # this moves the expression over
    }
  }
  
  
  ### Optionally order the table if there's a specific order you want
  if(is.null(order) == FALSE){
    chrom_expression_table <- chrom_expression_table[ order ]
  }
  
  return(chrom_expression_table)
  
}


#' ChromE_plot
#' Generates the ChromE heatmap
#' 
#' @param chrom_expression_table 
#' @param chr_table 
#' @param plot_annotations 
#' @param chr 
#' @param k 
#' @param order 
#' @param func 
#' @param plot.range 
#'
#' @return
#' @export
#'
#' @examples
ChromE_plot <- function(chrom_expression_table, chr_table, plot_annotations = NULL, chr, k = 10, order = NULL, func = "mean", plot.range = c(0,50)){
  
  if(func =="mean"){
    ### Sliding Window to better see the expression ###
    for(cell_type in names(chrom_expression_table)){
      chrom_expression_table[[cell_type]]$expression <- rollmean(chrom_expression_table[[cell_type]]$expression,
                                                                 k = k,
                                                                 fill = 0,
                                                                 partial = T)
    }
  } else if(func == "max"){
    for(cell_type in names(chrom_expression_table)){
      chrom_expression_table[[cell_type]]$expression <- rollmax(chrom_expression_table[[cell_type]]$expression,
                                                                k = k,
                                                                fill = 0,
                                                                partial = T)
    } 
  } else stop("func can only equal mean or max")
  
  chrom_expression_table <- chrom_expression_table[plot_annotations] # filter and only keep cell types of interest
  
  ### Split each cell type based on chromosome ###
  for(cell_type in names(chrom_expression_table)){
    chrom_expression_table[[cell_type]] <- split(chrom_expression_table[[cell_type]], chrom_expression_table[[cell_type]]$chr) # split based on the chromosome
    
    for(c in names(chrom_expression_table[[cell_type]])){
      chrom_expression_table[[cell_type]][[c]] <- arrange(chrom_expression_table[[cell_type]][[c]], order) # re-order just in case
    }
    
  }
  
  ### Optionally order the table if there's a specific order you want
  if(is.null(order) == FALSE){
    chrom_expression_table <- chrom_expression_table[ order ]
  }
  
  
  ### Combine Based on the chromosome we want to keep ###
  chr_expr_table <- data.frame(row.names = rownames(chrom_expression_table[[1]][[chr]]))
  
  for(cell_type in names(chrom_expression_table)){
    chr_expr_table[cell_type] <- chrom_expression_table[[cell_type]][[chr]]$expression
  }
  
  
  ### Cleanup Table for Plotting ###
  
  chr_expr_table <- t(chr_expr_table)
  
  
  ### Plotting Table ###
  
  # Catch if input plot range is wrong
  if(length(col.range) != 2 & class(col.range) != "numeric") stop("Error: the plot range you've given might be in the wrong format. Input should look like: c(min, max)")
  
  ## Setting range of the heatmap
  if(is.null(plot.range) == FALSE){
    col.range <- plot.range
  }else{
    all.range <- c()
    for(cell_type in names(chrom_expression_table)){
      all.range <- c(all.range, range(chrom_expression_table[[cell_type]]$expression))
    }
    col.range <- range(all.range) # Get the range of numbers from the table
  }
  
  
  col_fun <- colorRamp2(seq(col.range[[1]], col.range[[2]], length.out = 50), viridis(50)) # make colors for the range of numbers
  
  # Plot  the individual chromosome 
  chr_plots = Heatmap(chr_expr_table, 
                      col = col_fun, 
                      cluster_rows = F, cluster_columns = F,
                      column_title = chr,
                      show_column_names = F, name = " ",
                      show_heatmap_legend = T,
                      border_gp = gpar(col= "black", lty = 1, lex = 2, lineend = "round", linejoin = "round"))
  
  
  
  return(chr_plots)  
}
  
  
ChromE_histogram <- function(chrom_expression_table, chr, k = 10, plot_annotations = NULL, order = NULL, func = "mean", plot.range = c(0,50) ){
  
  ## Rolling window to amplify signal ##
  if(func =="mean"){
    ### Sliding Window to better see the expression ###
    for(cell_type in names(chrom_expression_table)){
      chrom_expression_table[[cell_type]]$expression <- rollmean(chrom_expression_table[[cell_type]]$expression,
                                                                 k = k,
                                                                 fill = 0,
                                                                 partial = T)
    }
  } else if(func == "max"){
    for(cell_type in names(chrom_expression_table)){
      chrom_expression_table[[cell_type]]$expression <- rollmax(chrom_expression_table[[cell_type]]$expression,
                                                                k = k,
                                                                fill = 0,
                                                                partial = T)
    } 
  } else stop("func can only equal mean or max")
  
  if(is.null(plot_annotations) == FALSE){
    chrom_expression_table <- chrom_expression_table[plot_annotations] # filter and only keep cell types of interest
  } else {plot_annotations <- names(chrom_expression_table)}
  
  ### Split each cell type based on chromosome ###
  for(cell_type in names(chrom_expression_table)){
    chrom_expression_table[[cell_type]] <- split(chrom_expression_table[[cell_type]], chrom_expression_table[[cell_type]]$chr) # split based on the chromosome
    
    for(c in names(chrom_expression_table[[cell_type]])){
      chrom_expression_table[[cell_type]][[c]] <- arrange(chrom_expression_table[[cell_type]][[c]], order) # re-order just in case
    }
    
  }
  
  ### Optionally order the table if there's a specific order you want
  if(is.null(order) == FALSE){
    chrom_expression_table <- chrom_expression_table[ order ]
  }
  
  
  ### Combine Based on the chromosome we want to keep ###
  chr_expr_table <- data.frame(row.names = rownames(chrom_expression_table[[1]][[chr]]))
  
  for(cell_type in names(chrom_expression_table)){
    chr_expr_table[cell_type] <- chrom_expression_table[[cell_type]][[chr]]$expression
  }
  
  chr_expr_table$order <- 1:nrow(chr_expr_table)
  
  return(chr_expr_table)
  
}


ChromE_single <- function(seurat_object, chr_table, annotation, ident.1, assay = "RNA", k = 200, plot.range = NULL){
  
  ### Check to make sure required packages installed ###
  if(require("Seurat") == F) stop("Package Error: Please install & load Seurat") # Checks for Seurat install
  if(require("zoo") == F) stop("Package Error: Please install & load zoom") # Checks for zoom install
  if(require("ComplexHeatmap") == F) stop("Package Error: Please install / load ComplexHeatmap") # Checks for ComplexHeatmap install
  
  ### Check for errors ###
  if(class(seurat_object) != "Seurat") stop("Error: the data is not given as a Seurat Object") # checks to make sure you're giving it a seurat object
  if("expression"  %in% colnames(chr_table) == F) stop("Error: The reference table does not have a column called expression") # checks to make sure there is an expression column
  if(annotation  %in% colnames(seurat_object@meta.data) == F)  stop("Error: the annotation you provided does not exist in the seurat metadata") # checks to make sure annotation is present in Seurat metadata
  
  ### Generate Table for Plotting ###
  # Generate a table with the average expression
  avg_expression <- as.data.frame(AverageExpression(seurat_object, group.by = annotation)[[assay]])
  
  
  # Add the gene expression from the average expression table in the order of the chromosome###
  chrom_expression_table <- list()
  
  for(cell_type in colnames(avg_expression)){
    chrom_expression_table[[cell_type]] <-  data.frame(gene = chr_table$gene,chr = chr_table$chr, order = chr_table$order, expression = chr_table$expression) # move genes from chr to list
    rownames(chrom_expression_table[[cell_type]]) <- chr_table$gene
    
    for(gene in rownames(avg_expression)){
      chrom_expression_table[[cell_type]][[gene, "expression"]] <-  avg_expression[[gene, cell_type]] # this moves the expression over
    }
    
  }
  
  # Generate a table to hold the delta
  delta.table <- chrom_expression_table[[ident.1]] #set to the first ident
  
  ## Setting range of the heatmap
  if(is.null(plot.range) == FALSE){
    col.range <- plot.range
  }else{
    col.range <- range(delta.table$expression) # Get the range of numbers from the table
  }
  # Catch if input plot range is wrong
  if(length(col.range) != 2 & class(col.range) != "numeric") stop("Error: the plot range you've given might be in the wrong format. Input should look like: c(min, max)")
  
  
  ### Sliding Window to better see the expression ###
  delta.table$expression <- rollmax(delta.table$expression, k = k, fill = 0)
  
  
 
  ## Split By chromosome
  delta.table <- split(delta.table, delta.table$chr)
  
  # Order in case
  chr_expr_table <- list()
  for(c in names(delta.table)){
    delta.table[[c]] <- arrange(delta.table[[c]], order)
    
    # move to the chr_expr_table
    chr_expr_table[[c]] <- data.frame(row.names = rownames(delta.table[[c]]))
    chr_expr_table[[c]][ident.1] <- delta.table[[c]]$expression
  }
  
  
  ### Cleanup Table for Plotting ###
  for(c in names(chr_expr_table)){
    chr_expr_table[[c]] <- t( chr_expr_table[[c]] )
  }
  
  ### Plotting Table ###
  col_fun <- colorRamp2(seq(col.range[[1]], col.range[[2]], length.out = 50), viridis(50)) # make colors for the range of numbers
  
  # Plot all the individual chromosome heatmaps
  chr_plots <- NULL # create an empty heatmap list object
  for(c in names(chr_expr_table)){
    
    # only keep last legend
    keep.legend = FALSE
    if(c == end(names(chr_expr_table))[[1]]){keep.legend = TRUE}
    chr_title <- c
    
    # setting x y & m
    if(length(chr_expr_table) == 25){ #25 would be human dataset!
      if(c == 23){chr_title <- "X"}else if(c == 24){chr_title <- "Y"} else if(c == 25){chr_title <- "M"}
    } 
    
    if(length(chr_expr_table) == 23){ #23 would be mouse dataset!
      if(c == 21){chr_title <- "X"}else if(c == 22){chr_title <- "Y"} else if(c == 23){chr_title <- "M"}
    }   
      
    # generate heatmap
    chr_plots = chr_plots + Heatmap(chr_expr_table[[c]], 
                                    col = col_fun, 
                                    cluster_rows = F, cluster_columns = F,
                                    column_title = chr_title,
                                    show_column_names = F, name = " ",
                                    show_heatmap_legend = keep.legend,
                                    border_gp = gpar(col= "black", lty = 1, lex = 2, lineend = "round", linejoin = "round"))
  }
  
  
  # Now combine all the individual heatmaps together
  plt <- draw(chr_plots, ht_gap = unit(5, "points"))
  
  
  return(plt)  
  
}

gene_count_boxplot <- function(seurat_object, split.by, order = NULL){
  
  # split based on the categories
  seurat_object.split <- SplitObject(seurat_object, split.by = split.by)
  
  if(is.null(order) == FALSE){
    seurat_object.split <- seurat_object.split[order] 
  }
  
  # Sum
  sum.list <- list()
  for(i in names(seurat_object.split)){
    for(cell in colnames(seurat_object.split[[i]]@assays$RNA)){
      sum.list[[i]] <- c(sum.list[[i]], sum(seurat_object.split[[i]]@assays$RNA[,cell] > 0))
    }
  }
  
  # join the lists together
  epi.table <- as.data.frame( do.call(cbind,sum.list) )
  
  epi.stack <- stack(epi.table)
  
  return(
    ggplot(epi.stack, aes(x = ind, y = values, fill = ind)) +
      geom_boxplot()
  )
}




