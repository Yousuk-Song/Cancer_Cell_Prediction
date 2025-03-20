#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)

# ✅ Command-line arguments 처리
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("❌ 사용법: Rscript change_gene_order.R <rds_path> <gene_list_path>")
}

rds_path <- args[1]  # 첫 번째 인자: RDS 파일 경로
gene_list_path <- args[2]  # 두 번째 인자: 정렬할 유전자 리스트 경로

cat("📌 입력 RDS 파일:", rds_path, "\n")
cat("📌 유전자 리스트 파일:", gene_list_path, "\n")

# ✅ 1. RDS 파일 로드
seurat_obj <- readRDS(rds_path)
seurat_obj <- UpdateSeuratObject(seurat_obj)  # Seurat v5 구조 업데이트

# ✅ 2. Count Matrix 가져오기 (Seurat v5 방식, 유전자가 colnames에 있음)
count_matrix <- GetAssayData(seurat_obj, layer = "counts")

# ✅ 3. 새로운 유전자 순서 로드
gene_order <- readLines(gene_list_path)

# ✅ 4. 현재 매트릭스와 유전자 리스트 비교
existing_genes <- colnames(count_matrix)  # 기존 count matrix의 유전자 목록
common_genes <- intersect(gene_order, existing_genes)  # 매트릭스에 존재하는 유전자
missing_genes <- setdiff(gene_order, existing_genes)  # 매트릭스에 없는 유전자 (추가 대상)
extra_genes <- setdiff(existing_genes, gene_order)  # gene list에 없는 유전자 (제거 대상)

# ✅ 5. 매트릭스에서 제거해야 할 유전자 제거
count_matrix <- count_matrix[, common_genes, drop = FALSE]

# ✅ 6. 없는 유전자 추가 (발현량 0으로 설정)
if (length(missing_genes) > 0) {
  zero_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = length(missing_genes))
  colnames(zero_matrix) <- missing_genes
  rownames(zero_matrix) <- rownames(count_matrix)

  # 기존 count matrix와 결합
  count_matrix <- cbind(count_matrix, zero_matrix)
}

# ✅ 7. 최종적으로 matrix 기준으로 유전자 순서 정렬
sorted_matrix <- count_matrix[, gene_order, drop = FALSE]  # 최종 정렬

# ✅ 8. 임시 Assay 생성 후 기본 Assay 변경 (기존 RNA 삭제 문제 해결)
seurat_obj[["temp"]] <- CreateAssayObject(counts = sorted_matrix)
DefaultAssay(seurat_obj) <- "temp"  # 기본 Assay를 임시 Assay로 변경
seurat_obj[["RNA"]] <- NULL  # 기존 RNA Assay 삭제

# ✅ 9. 새로운 RNA Assay 생성 및 적용
new_assay <- CreateAssayObject(counts = sorted_matrix)
seurat_obj[["RNA"]] <- new_assay  # 새로운 RNA Assay 적용
DefaultAssay(seurat_obj) <- "RNA"  # 기본 Assay를 다시 RNA로 변경
seurat_obj[["temp"]] <- NULL  # 임시 Assay 제거


# ✅ 11. Seurat 내부 메모리 최적화 (불필요한 데이터 삭제)
seurat_obj <- DietSeurat(seurat_obj, assays = "RNA")

# ✅ 12. 새로운 RDS 파일 저장 (원본 파일 이름 유지)
output_rds_path <- gsub(".rds", "_sorted.rds", rds_path)
saveRDS(seurat_obj, output_rds_path)

# ✅ 13. 최종 matrix에서 추출한 gene list 저장 (원본 파일 이름 유지)
final_gene_list <- colnames(seurat_obj@assays$RNA$counts)
output_gene_list_path <- gsub(".txt", "_sorted.txt", gene_list_path)
writeLines(final_gene_list, output_gene_list_path)

cat("✅ 유전자 순서 강제 적용 완료!\n")
cat(paste0("📌 정렬된 RDS 파일 저장: ", output_rds_path, "\n"))
cat(paste0("📌 정렬된 유전자 리스트 저장: ", output_gene_list_path, "\n"))
cat(paste0("📌 추가된 유전자 수: ", length(missing_genes), "\n"))
cat(paste0("📌 제거된 유전자 수: ", length(extra_genes), "\n"))
