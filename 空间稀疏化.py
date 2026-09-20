import pandas as pd
import numpy as np
import os

def scientific_spatial_thinning(input_csv, output_csv, decimalLatitude, decimalLongitude, resolution_min=2.5):
    """
    基於環境柵格分辨率的空間稀疏化腳本
    
    參數:
    input_csv: 原始 CSV 文件路徑
    output_csv: 稀疏化後保存的路徑
    lat_col: 緯度列名
    lon_col: 經度列名
    resolution_min: 分辨率（分鐘），預設為 2.5
    """
    
    # 1. 讀取數據
    df = pd.read_csv(input_csv)
    original_count = len(df)
    
    # 2. 將分鐘分辨率轉換為十進制度數
    # 2.5 minutes = 0.0416666667 degrees
    res_degree = resolution_min / 60.0
    
    # 3. 創建網格索引 (Grid Binning)
    # 這是最科學的方法：將坐標除以分辨率後取整，相同網格的點會得到相同的索引值
    df['grid_lat'] = np.floor(df[decimalLatitude] / res_degree)
    df['grid_lon'] = np.floor(df[decimalLongitude] / res_degree)
    
    # 4. 執行稀疏化：在每個網格組中只保留第一個點
    # drop_duplicates 會保留重複項中的第一行，移除其餘行
    thinned_df = df.drop_duplicates(subset=['grid_lat', 'grid_lon'])
    
    # 5. 移除輔助列並保存
    thinned_df = thinned_df.drop(columns=['grid_lat', 'grid_lon'])
    thinned_df.to_csv(output_csv, index=False)
    
    # 6. 輸出報告
    thinned_count = len(thinned_df)
    reduction = original_count - thinned_count
    print(f"--- 稀疏化報告 ({input_csv}) ---")
    print(f"原始點位數: {original_count}")
    print(f"保留點位數: {thinned_count}")
    print(f"移除重複點: {reduction} (約 {reduction/original_count:.1%})")
    print(f"使用的分辨率網格: {resolution_min} arc-minutes ({res_degree:.6f} degrees)")
    print(f"結果已保存至: {output_csv}\n")

# ================= 使用示例 =================
# 請確保你的列名與 CSV 中完全一致（注意大小寫）
# 如果你的 CSV 列名是 'decimalLatitude'，請相應修改

file_to_process = r"E:\distribution\Occurence_test\Solenopsis_invicta.csv"
output_file = r"E:\distribution\Occurence_test\Solenopsis_invicta_thinned_2.5m.csv"

scientific_spatial_thinning(
    input_csv=file_to_process, 
    output_csv=output_file, 
    decimalLatitude='decimalLatitude', 
    decimalLongitude='decimalLongitude', 
    resolution_min=2.5
)