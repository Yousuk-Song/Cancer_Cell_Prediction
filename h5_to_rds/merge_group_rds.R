#!/usr/bin/env Rscript

# 필요한 패키지 로드
library(Seurat)
library(dplyr)

# 명령줄 인자로 그룹명 및 RDS 파일 리스트 받기
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript merge_rds_v5.R <group_name> <rds_file1> <rds_file2> ...")
}

group_name <- args[1]  # 그룹명 (ex. GSE141383)
rds_files <- args[-1]  # 나머지는 RDS 파일 리스트

# 🔹 첫 번째 RDS 파일을 읽어 초기 Seurat 객체 설정
print(paste("📂 Loading first RDS:", rds_files[1]))
seurat_list <- list()
seurat_list[[1]] <- readRDS(rds_files[1])

# 🔹 나머지 RDS 파일을 읽고 데이터 병합 준비
for (i in 2:length(rds_files)) {
  print(paste("📂 Loading and merging RDS:", rds_files[i]))
  seurat_list[[i]] <- readRDS(rds_files[i])
}

# 🔹 Expression Matrix 병합 (Seurat v5 호환)
#merged_counts <- do.call(cbind, lapply(seurat_list, function(obj) GetAssayData(obj, slot = "counts")))
merged_counts <- do.call(cbind, lapply(seurat_list, function(obj) GetAssayData(obj, layer = "counts")))


# 🔹 Metadata 병합 (단순 row-wise merge)
merged_metadata <- do.call(rbind, lapply(seurat_list, function(obj) obj@meta.data))

# 🔹 병합된 데이터로 새로운 Seurat 객체 생성
merged_seurat <- CreateSeuratObject(counts = merged_counts, meta.data = merged_metadata)

# 🔹 저장
output_rds <- paste0(group_name, "_merged.rds")
saveRDS(merged_seurat, file = output_rds)

print(paste("✅ Merged Seurat object saved as:", output_rds))
