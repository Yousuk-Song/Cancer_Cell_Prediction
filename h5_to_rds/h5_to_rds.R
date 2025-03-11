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
if (length(args) != 4) {
  stop("Usage: Rscript h5_to_seurat_full.R <h5_file> <metadata_tsv> <sample_name> <cancer_type>")
}

h5_file <- args[1]        # HDF5 파일 경로
metadata_tsv <- args[2]   # 메타데이터 TSV 파일
sample_name <- args[3]    # 샘플명 (Python 코드에서 name에 해당)
cancer_type <- args[4]    # 암종 ex) BRCA, AML, AEL

#h5_file <- "/data/processed_data/scRSEQ_AML/DISCO/BATCH/batch_1/GSM2758472.h5"        # HDF5 파일 경로
#metadata_tsv <- "/data/processed_data/scRSEQ_AML/DISCO/BATCH/MetaData/Glioma_GSE103224_CellMetainfo_table.tsv"   # 메타데이터 TSV 파일
#sample_name <- "GSM2758471"    # 샘플명 (Python 코드에서 name에 해당)
#cancer_type <- "Glioma"    # 암종 ex) BRCA, AML, AEL

# 파일명만 추출

meta_filename <- basename(metadata_tsv)

# "_" 기준으로 split
split_parts <- strsplit(meta_filename, "_")[[1]]

# "_" 기준으로 split하여 두 번째 요소 추출
project_id <- strsplit(meta_filename, "_")[[1]][2]


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

# 🔹 `Cell`에서 샘플 정보(`Sample`)와 `Cell ID` 분리
metadata <- metadata %>%
  mutate(
    Sample = ifelse(str_detect(Cell, "@"), str_extract(Cell, "^[^@]+"), sample_name),  # `@` 앞부분을 `Sample` 컬럼으로 저장 (없으면 sample_name 사용)
    Cell = ifelse(str_detect(Cell, "@"), str_extract(Cell, "(?<=@).*"), Cell)  # `@` 이후 부분만 `Cell`에 저장 (없으면 그대로 유지)
  )

# 🔹 `Sample`이 `sample_name`과 일치하는 경우만 필터링
metadata <- metadata %>%
  filter(Sample == sample_name)

# 🔹 Expression Matrix의 Cell 바코드 정리 (`sample_name_` 추가)
cell_barcodes_clean <- gsub("-1$", "", sub(".*@", "", cell_barcodes))
cell_barcodes_renamed <- paste0(sample_name, "_", cell_barcodes_clean)
colnames(expr_matrix) <- cell_barcodes_renamed

# 🔹 메타데이터의 `Cell`에도 `sample_name_` 추가
metadata$Cell <- paste0(sample_name, "_", metadata$Cell)

# 🔹 Expression Matrix와 메타데이터에서 공통된 Cell만 필터링
common_cells <- intersect(colnames(expr_matrix), metadata$Cell)
expr_matrix_filtered <- expr_matrix[, common_cells]
metadata_filtered <- metadata %>% filter(Cell %in% common_cells)

# 🔹 불필요한 열 제거
metadata_filtered <- metadata_filtered %>%
  select(-Sample)  # `Sample`은 이제 필요 없으면 제거 가능

# 🔹 Cancer Type 컬럼 추가 (모든 셀에 동일한 값)
metadata_filtered$CancerType <- cancer_type

# 🔹 Seurat 객체 생성
rownames(metadata_filtered) <- metadata_filtered$Cell
metadata_filtered <- metadata_filtered %>% select(-Cell)  # `Cell` 컬럼 삭제 (이미 rownames로 이동)

seurat_obj <- CreateSeuratObject(counts = expr_matrix_filtered, meta.data = metadata_filtered)

# 🔹 orig.ident 값 설정
seurat_obj@meta.data$orig.ident <- sample_name

# 🔹 RDS 파일 저장 (H5 파일명 기반)
# .h5 확장자를 .rds로 변경
output_rds <- gsub("\\.h5$", ".rds", h5_file)

# 경로와 파일명 분리
output_dir <- dirname(output_rds)  # 경로 추출
output_filename <- basename(output_rds)  # 파일명 추출

# 최종 경로 조합
output_rds <- file.path(output_dir, paste0(project_id, "_", output_filename))

# 결과 출력
print(output_rds)

saveRDS(seurat_obj, file = output_rds)

# 결과 메시지 출력
print(paste("✅ Seurat object has been successfully created and saved as", output_rds))

# HDF5 파일 닫기
h5$close()
