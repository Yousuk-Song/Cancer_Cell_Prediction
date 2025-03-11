#!/usr/bin/env Rscript

# 🔹 필요한 패키지 로드
library(hdf5r)
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(data.table)

# 🔹 명령줄 인자로 HDF5 파일 및 메타데이터 TSV 파일 받기
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript h5_to_seurat_full.R <h5_file> <metadata_tsv> <sample_name>")
}

h5_file <- args[1]        # HDF5 파일 경로
metadata_tsv <- args[2]   # 메타데이터 TSV 파일
sample_name <- args[3]    # 샘플명 (Python 코드에서 name에 해당)

#h5_file <- "/data/processed_data/scRSEQ_AML/DISCO/BATCH/batch_1/GSM4476485.h5"        # HDF5 파일 경로
#metadata_tsv <- "/data/processed_data/scRSEQ_AML/DISCO/BATCH/MetaData/BRCA_GSE148673_CellMetainfo_table.tsv"   # 메타데이터 TSV 파일
#sample_name <- "GSM4476485"    # 샘플명 (Python 코드에서 name에 해당)


# 🔹 HDF5 파일 읽기
h5 <- H5File$new(h5_file, mode = "r")

# 🔹 유전자 및 세포 바코드 불러오기
gene_names <- h5[["matrix/features/name"]][]
cell_barcodes <- h5[["matrix/barcodes"]][]

# 🔹 Sparse Matrix (CSC 형식) 불러오기
data <- h5[["matrix/data"]][]
indices <- h5[["matrix/indices"]][]
indptr <- h5[["matrix/indptr"]][]
num_genes <- h5[["matrix/shape"]][1]
num_cells <- h5[["matrix/shape"]][2]

# 🔹 희소 행렬 변환 (dgCMatrix 형식, 행: 유전자, 열: 세포)
expr_matrix <- sparseMatrix(
  i = indices + 1,
  p = indptr,
  x = data,
  dims = c(num_genes, num_cells),
  dimnames = list(gene_names, cell_barcodes)
)

# 🔹 메타데이터 불러오기
metadata <- fread(metadata_tsv, data.table = FALSE)

# 🔹 메타데이터 컬럼 이름 정리 (공백 및 특수문자 대체)
colnames(metadata) <- colnames(metadata) %>%
  str_replace_all("[()]", "") %>%  # 괄호 제거
  str_replace_all("\\s+", "_")      # 공백을 언더바(_)로 변환

# 🔹 Cell ID 정리 (`-1` 제거 및 `@` 처리)
metadata <- metadata %>%
  mutate(
    Cell.x = Cell,  # 원본 Cell ID 저장
    Sample = str_extract(Cell, "^[^@]+"),  # @ 앞부분 추출 (샘플 정보)
    Cell_suffix = str_extract(Cell, "(?<=@).*"),  # @ 뒤의 값만 유지
    Cell = gsub("-1$", "", Cell_suffix)  # 최종 Cell ID 정리
  )

# 🔹 Expression Matrix의 Cell 바코드도 동일하게 정리
cell_barcodes_clean <- gsub("-1$", "", sub(".*@", "", cell_barcodes))
colnames(expr_matrix) <- cell_barcodes_clean

# 🔹 모든 Cell ID에 `sample_name_` 추가
cell_barcodes_renamed <- paste0(sample_name, "_", cell_barcodes_clean)
metadata$Cell <- paste0(sample_name, "_", metadata$Cell)

# 🔹 Expression Matrix 컬럼 이름 변경
colnames(expr_matrix) <- cell_barcodes_renamed

# 🔹 Expression Matrix와 메타데이터에서 공통된 Cell만 필터링
common_cells <- intersect(colnames(expr_matrix), metadata$Cell)
common_cells

expr_matrix_filtered <- expr_matrix[, colnames(expr_matrix) %in% common_cells]
metadata_filtered <- metadata %>% filter(Cell %in% common_cells)

# 🔹 불필요한 열 제거
metadata_filtered <- metadata_filtered %>%
  select(-Cell.x, -Cell_suffix)


# 🔹 Seurat 객체 생성
seurat_obj <- CreateSeuratObject(counts = expr_matrix_filtered, meta.data = metadata_filtered)

seurat_obj@meta.data$orig.ident <- sample_name

# 🔹 RDS 파일 저장 (H5 파일명 기반)
output_rds <- gsub("\\.h5$", ".rds", h5_file)
saveRDS(seurat_obj, file = output_rds)

# 결과 메시지 출력
print(paste("✅ Seurat object has been successfully created and saved as", output_rds))

# HDF5 파일 닫기
h5$close()

