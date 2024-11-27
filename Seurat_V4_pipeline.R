#!/usr/bin/Rscript --vanilla

options(warn=-1)

suppressMessages(library(optparse))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(harmony))
suppressMessages(library(sctransform))
suppressMessages(library(MAST))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggprism))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(future.apply))
suppressMessages(library(parallel))

plan("multicore", workers = parallel::detectCores()-2)
# plan("multisession", workers = parallel::detectCores()-2)
# sprintf("FUTURE CURRENT WORKERS = %s", nbrOfWorkers())
options(future.globals.maxSize = 1000 * 1024^3)

# Define arguments
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Path to the parent folder containing CellRanger output sample directories. Required!"),
  make_option(c("-m","--meta_data"), type = "character", help = "Path to meta_data txt file. Required! The first column must be file name."),
  make_option(c("-o","--output_dir"), type = "character", default = "./", help = "Path to save the output RData file. Default: the current folder."),
  make_option(c("-s","--samplename"), type = "character", default = "sample", help = "Name of the sample. Default: sample."),
  make_option(c("-c","--copykat_path"), type = "character", default = NULL, help = "Path to the folder containing copykat files. Default: NULL."),
  make_option("--min_cells", type = "integer", default = 3, help = "Minimum cells for CreateSeuratObject. Default: 3."),
  make_option("--min_features", type = "integer", default = 200, help = "Minimum features for CreateSeuratObject. Default: 200."),
  make_option("--nFeature_RNA_min", type = "integer", default = 200, help = "Minimum nFeature_RNA for filtering. Default: 200."),
  make_option("--nFeature_RNA_max", type = "integer", default = 6000, help = "Maximum nFeature_RNA for filtering. Default: 6000."),
  make_option("--nCount_RNA_min", type = "integer", default = 300, help = "Minimum nCount_RNA for filtering. Default: 300."),
  make_option("--nCount_RNA_max", type = "integer", default = 30000, help = "Maximum nCount_RNA for filtering. Default: 30000."),
  make_option(c("-p","--percent_mt"), type = "integer", default = 10, help = "Maximum percent.mt for filtering. Default: 10."),
  make_option(c("-d","--dim_use"), type = "integer", default = 30, help = "Number of dimensions to use. Default: 30."),
  make_option(c("-v","--var_features"), type = "integer", default = 2000, help = "Number of variable features for FindVariableFeatures. Default: 2000."),
  make_option(c("-n","--npcs"), type = "integer", default = 30, help = "Number of PCs for RunPCA. Default: 30."),
  make_option(c("-r","--resolution"), type = "numeric", default = 0.5, help = "Resolution for FindClusters. Default: 0.5."),
  make_option(c("-g","--harmony_group_by"), type = "character", default = "orig.ident", help = "Group by variable for RunHarmony. Default: orig.ident.")
)

opt <- parse_args(OptionParser(usage="Rscript %prog -i input_dir -m /your/path/to/metadata.txt -o output_dir -p 25 -d 30 -v 2000 -n 30 -r 0.5",
                               description="\nSeurat_V4 to QC and Cluster scRNA-seq data.",
                               option_list = option_list))

# Read sample directories
sample_inputdir <- list.dirs(opt$input_dir, full.names = TRUE, recursive = FALSE)

print("Sample directories detected: ")
print(sample_inputdir)

# Initialize a list to store Seurat objects
seurat_list <- list()

# Go through each sample directory
# Metadata if available, Process each sample
if(is.null(opt$meta_data)){
  
  print("Metadata is not provided! Required!")
  
  } else if (!is.null(opt$meta_data) && file.exists(opt$meta_data)) {
  metadata <- read.table(opt$meta_data, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(metadata) <- c("file","study","sample","tissue","disease")
  for (i in 1:nrow(metadata)) {
    sample_name <- metadata[i, 1]
    
    # Create Seurat object
    sample_dir <- file.path(opt$input_dir, sample_name)
    
    if (dir.exists(sample_dir)) {
      message(paste("Processing sample:", sample_name))
      
      seurat_data <- Read10X(data.dir = sample_dir)
      seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                       project = sample_name,
                                       min.cells = opt$min_cells,
                                       min.features = opt$min_features)
    
      # Add metadata information
      for (j in 2:ncol(metadata)) {
        seurat_obj@meta.data[[colnames(metadata)[j]]] <- metadata[i, j]
      }
      
      # Add copykat data if available
      if(is.null(opt$copykat_path)){
        print("copykat path not provided!")
      } else if (!is.null(opt$copykat_path)){
        copykat_file <- file.path(opt$copykat_path, paste0(sample_name, "_copykat_prediction.txt"))
        copykat_data <- read.table(copykat_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        seurat_obj@meta.data$copykat.pred <- "unknown"
        copykat_cells <- intersect(rownames(seurat_obj@meta.data), copykat_data$cell.names)
        seurat_obj@meta.data[copykat_cells, "copykat.pred"] <- copykat_data[match(copykat_cells, copykat_data$cell.names), "copykat.pred"]
      }
      
      # Replace underscores in feature names to prevent warnings
      rownames(seurat_obj) <- gsub("_", "-", rownames(seurat_obj))
      
      # Add percent.mt
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      
      # Filter cells
      seurat_obj <- subset(
        seurat_obj,
        subset = nFeature_RNA > opt$nFeature_RNA_min & 
          nFeature_RNA < opt$nFeature_RNA_max & 
          nCount_RNA > opt$nCount_RNA_min & 
          nCount_RNA < opt$nCount_RNA_max & 
          percent.mt < opt$percent_mt
      )
      
      seurat_list[[sample_name]] <- seurat_obj
      
    } else {
      # print(paste0("Warning: Sample directory for", sample_name, "not found!"))
      warning(paste("Sample directory not found for:", sample_name))
    }
  }
}

print(paste0(opt$samplename," has ",length(seurat_list)," samples!"))

# Combine all Seurat objects into one

# combined_seurat <- Reduce(function(x, y) merge(x, y, add.cell.ids = names(seurat_obj)), seurat_obj)
# combined_seurat <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

if (length(seurat_list) == 1) {
  # If seurat_list has only one sample, assign it directly to combined_seurat and give the cell number
  seurat_list[[1]] <- RenameCells(seurat_list[[1]], add.cell.id = names(seurat_list))
  combined_seurat <- seurat_list[[1]]
} else {
  # If there are multiple samples, use the merge function to merge them and add numbers at the same time
  combined_seurat <- merge(
    x = seurat_list[[1]], 
    y = seurat_list[-1], 
    add.cell.ids = names(seurat_list)
  )
}

print("Merged sample before standard process: ")
print(combined_seurat)

# Process combined Seurat object
# check the number of kind of orig.ident
if (length(unique(combined_seurat$orig.ident)) != 1) {
  # if type of orig.ident >1，run all steps
  combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = opt$var_features, verbose = F)
  combined_seurat <- ScaleData(combined_seurat, features = rownames(combined_seurat)) # features = VariableFeatures(.) ,vars.to.regress = c("nCount_RNA","percent.mt")
  combined_seurat <- RunPCA(combined_seurat, npcs = opt$npcs, features = VariableFeatures(combined_seurat), verbose = F)
  combined_seurat <- RunHarmony(combined_seurat, group.by.vars = opt$harmony_group_by, verbose = F)
  combined_seurat <- FindNeighbors(combined_seurat, dims = 1:opt$dim_use, verbose = F)
  combined_seurat <- FindClusters(combined_seurat, resolution = opt$resolution, verbose = F)
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:opt$dim_use, reduction = "harmony", seed.use = 1, verbose = F) # umap.method = "umap-learn",metric = "correlation"
  combined_seurat <- RunTSNE(combined_seurat, dims = 1:opt$dim_use, reduction = "harmony", seed.use = 1, verbose = F)
} else {
  # if type of orig.ident =1，run to RunPCA step
  combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = opt$var_features, verbose = F)
  combined_seurat <- ScaleData(combined_seurat, features = rownames(combined_seurat))
  combined_seurat <- RunPCA(combined_seurat, npcs = opt$npcs, features = VariableFeatures(combined_seurat), verbose = F)
}

print("Merged sample after standard process: ")
print(combined_seurat)

# Save Seurat objects and combined object
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

assign(paste0(opt$samplename, "_seurat_list"), seurat_list)
assign(paste0(opt$samplename, "_combined_seurat"), combined_seurat)

save(
  list = c(paste0(opt$samplename, "_seurat_list"), 
           paste0(opt$samplename, "_combined_seurat")),
  file = file.path(opt$output_dir, paste0(opt$samplename, ".rda"))
)

plan("sequential")

print(paste0("Congratulations! ",opt$samplename," has finished!"))


