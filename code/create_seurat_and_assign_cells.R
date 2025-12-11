library(Seurat)
library(dplyr)

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



