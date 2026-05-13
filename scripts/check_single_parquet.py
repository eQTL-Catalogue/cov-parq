import pyarrow.parquet as pq
import pandas as pd
import argparse
import sys
import os

def check_parquet_nas(file_path):
    """
    Checks a parquet file for NA/null values by iterating over row groups
    to avoid high memory usage.
    """
    print(f"Checking file: {file_path}", flush=True)
    
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}", flush=True)
        return

    try:
        parquet_file = pq.ParquetFile(file_path)
        num_row_groups = parquet_file.num_row_groups
        print(f"File has {num_row_groups} row groups.", flush=True)
        
        total_rows = 0
        na_found = False
        na_columns = {}

        for i in range(num_row_groups):
            # Read one row group at a time
            table = parquet_file.read_row_group(i)
            df = table.to_pandas()
            
            total_rows += len(df)
            
            # Check for NAs in this chunk
            if df.isna().any().any():
                na_found = True
                chunk_na_counts = df.isna().sum()
                chunk_na_cols = chunk_na_counts[chunk_na_counts > 0]
                
                for col, count in chunk_na_cols.items():
                    na_columns[col] = na_columns.get(col, 0) + count
                
                # Print first few rows with NAs from this chunk if it's the first time finding them
                if len(na_columns) == len(chunk_na_cols): # First time finding NAs
                     print(f"\n[Row Group {i}] First 5 rows with NAs:", flush=True)
                     print(df[df.isna().any(axis=1)].head(), flush=True)

            # Print progress every 10 row groups
            if (i + 1) % 10 == 0:
                print(f"Processed {i + 1}/{num_row_groups} row groups...", end='\r', flush=True)

        print(f"\nProcessed {num_row_groups} row groups.", flush=True)
        
        if na_found:
            print("\nWARNING: NA values found!", flush=True)
            print("-" * 30, flush=True)
            print(f"Total rows checked: {total_rows}", flush=True)
            print(f"Columns with NAs and total count:\n", flush=True)
            for col, count in na_columns.items():
                print(f"{col}: {count}", flush=True)
        else:
            print("\nSUCCESS: No NA values found in the file.", flush=True)
            print(f"Checked {total_rows} rows.", flush=True)

    except Exception as e:
        print(f"\nError reading parquet file: {e}", flush=True)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check a parquet file for NA values (memory efficient).")
    parser.add_argument("file_path", nargs="?", 
                        default="/gpfs/helios/projects/eQTLCatalogue/coverage_parquet/MacroMap/CIL_6.parquet",
                        help="Path to the parquet file to check")
    
    args = parser.parse_args()
    check_parquet_nas(args.file_path)
