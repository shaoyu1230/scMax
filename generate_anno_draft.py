#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
利用 scMax_experience 数据库进行单细胞智能辅助注释草图生成器
"""

import os
import sys
import yaml
import sqlite3
import argparse
import pandas as pd

def safe_db_connect(db_path):
    if not os.path.exists(db_path):
        print(f"警告: 经验数据库文件未找到 -> {db_path}")
        print("将生成一份全部为空的模板文件。")
        return None
    try:
        conn = sqlite3.connect(db_path)
        return conn
    except sqlite3.Error as e:
        print(f"数据库连接失败: {e}")
        return None

def fetch_experience(conn, species=None, tissue_type=None, annotation_level=None):
    """
    基于给定的物种/组织从经验库拉取所有的已知 Marker 矩阵。
    """
    if not conn:
        return []
    
    query = "SELECT cell_type, cell_marker, reference, description FROM marker_experience WHERE 1=1"
    params = []
    
    if species and species != "Unknown":
        query += " AND species = ?"
        params.append(species)
        
    if tissue_type and tissue_type != "Unknown":
        query += " AND tissue_type = ?"
        params.append(tissue_type)
        
    if annotation_level:
        query += " AND annotation_level = ?"
        params.append(annotation_level)
    
    cursor = conn.cursor()
    cursor.execute(query, params)
    rows = cursor.fetchall()
    
    exp_list = []
    for r in rows:
        c_type, markers_str, ref, desc = r
        # 解析出这个细胞类型下面挂载的所有特征基因集
        markers = [m.strip().upper() for m in str(markers_str).split(',') if m.strip()]
        exp_list.append({
            'cell_type': c_type,
            'markers': set(markers),
            'reference': ref if ref else "",
            'description': desc if desc else ""
        })
    return exp_list

def infer_celltype(cluster_markers, exp_list):
    """
    打分算法：用当前簇的高变基因去碰撞经验库。
    返回最高分数的那个匹配记录。
    """
    if not exp_list or not cluster_markers:
        return None
    
    best_match = None
    max_score = 0
    
    c_markers_set = set([m.strip().upper() for m in cluster_markers])
    
    for exp in exp_list:
        # P-value 最高的 Top N 基因与经验库里的基因交集数
        overlap = c_markers_set.intersection(exp['markers'])
        score = len(overlap)
        
        if score > max_score:
            max_score = score
            best_match = exp
            
    # 如果完全没有交集，就不瞎猜了
    if max_score > 0:
        return best_match
    return None

def generate_draft(markers_file, config_file, output_file, top_n=10):
    # 1. 解析核心配置
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
        
    db_conf = config.get("Database_Info", {})
    db_path = db_conf.get("db_path", "")
    # 如果没配的话尝试找相对目录默认备胎
    if not db_path:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        db_path = os.path.join(base_dir, "WebUI", "scMax_projects.db")
        
    species = db_conf.get("species", "Unknown")
    tissue = db_conf.get("tissue_type", "Unknown")
    anno_level = db_conf.get("annotation_level", "")
    
    print(f"--- 智能辅助注释初始化 ---")
    print(f"经验库引擎: {db_path}")
    print(f"过滤条件 -> 物种: {species} | 组织: {tissue} | 层次缩影: {anno_level}")
    
    # 2. 读取差异基因表 (一般由 Seurat FindAllMarkers 生成, 必有 cluster 和 gene列)
    try:
        df_markers = pd.read_csv(markers_file, sep=None, engine='python')
    except Exception as e:
        print(f"读取 Markers 文件失败 {markers_file}: {e}")
        sys.exit(1)
        
    # 检查核心列名是否存在
    if 'cluster' not in df_markers.columns or 'gene' not in df_markers.columns:
        # Seurat 转 CSV 偶尔会把 gene 放 index，查一查是不是没有叫 gene 的列
        # 退回使用首列为基因的方式尝试挽救
        if 'cluster' in df_markers.columns:
             df_markers['gene'] = df_markers.iloc[:, 0]
        else:
             print("Error: 输入文件中必须含有 'cluster' 分型列与 'gene' 基因列！")
             sys.exit(1)
            
    # 3. 按簇切分并获取每个簇打分最高的前 Top N 基因
    conn = safe_db_connect(db_path)
    experience_pool = fetch_experience(conn, species, tissue, anno_level) if conn else []
    if experience_pool:
        print(f"成功加载 {len(experience_pool)} 条相关经验数据...")
    else:
        print("当前条件下无任何经验记录可供参考比对，将退化为纯白板模板。")
        
    draft_rows = []
    
    # 按 p_val_adj 排序以防它没有排序
    if 'p_val_adj' in df_markers.columns:
        df_markers = df_markers.sort_values(by=['cluster', 'p_val_adj'])
        
    clusters = df_markers['cluster'].unique()
    
    for c in clusters:
        c_df = df_markers[df_markers['cluster'] == c].head(top_n)
        c_top_genes = c_df['gene'].tolist()
        
        # 将待画图的 Markers 以逗号拼成一列
        display_markers = ",".join(c_top_genes[:5]) # 默认提取前5个最显著画点图即可，放太多图难看
        
        best_exp = infer_celltype(c_top_genes, experience_pool)
        
        if best_exp:
            pred_type = best_exp['cell_type']
            ref = best_exp['reference']
            desc = best_exp['description']
        else:
            pred_type = ""
            ref = ""
            desc = ""
            
        draft_rows.append({
            "Cluster": c,
            "CellType": pred_type,
            "Markers": display_markers,
            "reference": ref,
            "description": desc
        })
        
    # 4. 生成组装稿
    df_draft = pd.DataFrame(draft_rows)
    
    try:
        if output_file.endswith('.csv'):
            df_draft.to_csv(output_file, index=False)
        else:
            df_draft.to_excel(output_file, index=False)
        print(f"✅ 智能辅助注释草稿已生成成功: {output_file}")
        print("请在查阅修改后，将此文件作为 05_celltype 的 `annofile` 输入。")
    except Exception as e:
        print(f"生写输出文件错误: {e}")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据经验库智能预测出具一版注释草图 (.xls 格式)。")
    parser.add_argument("-i", "--input_markers", required=True, help="集群特异基因文件 (一般为 03或04的 markers.csv)")
    parser.add_argument("-c", "--config", required=True, help="YAML 配置文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出草稿名，如: my_draft_anno.xlsx")
    parser.add_argument("-n", "--topn", type=int, default=15, help="取前多少个差异基因用来跟库里的比对，默认 15 个")
    
    args = parser.parse_args()
    
    generate_draft(args.input_markers, args.config, args.output, top_n=args.topn)
