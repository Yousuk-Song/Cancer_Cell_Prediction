#!/usr/bin/python

import pandas as pd
import sys
import re

def replace_cell_ids(input_file, output_file, replacements):
    # 파일 읽기
    df = pd.read_csv(input_file, sep='\t')

    # Cell 컬럼의 @ 앞부분만 변경
    def modify_cell(cell):
        parts = cell.split('@')
        if len(parts) == 2 and parts[0] in replacements:
            parts[0] = replacements[parts[0]]
        return '@'.join(parts)

    df.iloc[:, 0] = df.iloc[:, 0].apply(modify_cell)

    # 파일 저장
    df.to_csv(output_file, sep='\t', index=False)
    print(f"변경 완료: {output_file}")

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = input_file.replace("_CellMetainfo_table.tsv", "_CellMetainfo_table.modified.tsv")

    # 치환할 문자열 딕셔너리 생성
    replacements = {}
    for arg in sys.argv[2:]:
        if arg.startswith("--replace="):
            old, new = arg.replace("--replace=", "").split("=")
            replacements[old] = new

    replace_cell_ids(input_file, output_file, replacements)


# example: python script.py AEL_GSE142213_CellMetainfo_table.tsv --replace=OX1164=GSM4222796
