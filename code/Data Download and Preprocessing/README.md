library(Seurat)
library(dplyr)
library(Azimuth)
library(SeuratData)


base_dir <- "/datos/lymphocytes/monocytes/ad/fastqfiles"
sample_dirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)


seurat_list <- lapply(sample_dirs, function(dir) {
  sample_name <- basename(dir)
  counts <- Read10X(data.dir = file.path(dir, "filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  seurat_obj$Sample <- sample_name
  return(seurat_obj)
})

seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1],
                    add.cell.ids = sapply(seurat_list, function(x) x$Sample))

saveRDS(seurat_obj, file = "seurat_obj.rds")
#------------------QC-------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

#---------------------------------------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

#--------------------------------------------------------------

metadata_df <- read.csv("metadata.csv", stringsAsFactors = FALSE)
metadata_df <- read.csv("metadata.csv", stringsAsFactors = FALSE)

metadata <- data.frame(
  row.names = colnames(seurat_obj),
  APOE = metadata_df$APOE[match(seurat_obj$Sample, metadata_df$Sample)],
  Condition = metadata_df$Condition[match(seurat_obj$Sample, metadata_df$Sample)]
)

seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)



# ----  Azimuth ----
SeuratData::InstallData("pbmcref")
reference <- SeuratData::LoadData("pbmcref")

seurat_obj <- RunAzimuth(seurat_obj, reference = "pbmcref")

saveRDS(seurat_obj, file = "seurat_obj_final.rds")
cat("✅ Final Seurat object saved as seurat_obj_final.rds\n")
                                          





                                          



