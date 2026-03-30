import os
import sqlite3
import argparse
import yaml
import glob
from datetime import datetime
import csv
import shutil

def init_db(db_path):
    """初始化数据库表结构"""
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # 核心结果库：projects
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS projects (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            std_project_id TEXT,
            custom_project_id TEXT,
            crm_id TEXT,
            annotation_level TEXT,
            species TEXT,
            species_zh TEXT,
            tissue_type TEXT,
            tissue_type_zh TEXT,
            disease_type TEXT,
            seq_tech TEXT,
            analyst TEXT,
            notes TEXT,
            outdir TEXT,
            umap_path TEXT,
            fraction_path TEXT,
            marker_table_path TEXT,
            anno_png_path TEXT,
            report_path TEXT,
            bk_umap TEXT,
            bk_fraction TEXT,
            bk_anno_png TEXT,
            bk_report TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # 动态热更新旧数据库
    for col in ['disease_type', 'seq_tech', 'anno_png_path', 'bk_anno_png']:
        try:
            cursor.execute(f"ALTER TABLE projects ADD COLUMN {col} TEXT")
        except sqlite3.OperationalError:
            pass
    
    # 经验沉淀库：marker_experience
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS marker_experience (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            custom_project_id TEXT,
            crm_id TEXT,
            annotation_level TEXT,
            species TEXT,
            tissue_type TEXT,
            cluster TEXT,
            cell_type TEXT,
            markers TEXT,
            reference TEXT,
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    conn.commit()
    return conn

def main():
    parser = argparse.ArgumentParser(description="scMax Database Manager for Cell Annotation Results")
    parser.add_argument('-c', '--config', required=True, help="Path to config_template.yaml")
    parser.add_argument('-o', '--outdir', required=True, help="Cell type output directory (e.g. 05_celltype_out)")
    args = parser.parse_args()

    # 读取配置文件
    with open(args.config, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)
        
    db_info = cfg.get('Database_Info', {})
    if not db_info.get('run', False):
        print("Database_Info run flag is false. Skipping DB write.")
        return

    # 提取信息
    db_path = db_info.get('db_path', 'scMax_projects.db') # 用户要求如果没配，可以在本地
    std_project_id = db_info.get('std_project_id', '')
    custom_project_id = db_info.get('custom_project_id', '')
    crm_id = db_info.get('crm_id', '')
    annotation_level = db_info.get('annotation_level', 'MajorType')
    species = db_info.get('species', 'Human')
    species_zh = db_info.get('species_zh', '')
    tissue_type = db_info.get('tissue_type', '')
    tissue_type_zh = db_info.get('tissue_type_zh', '')
    disease_type = db_info.get('disease_type', '')
    seq_tech = db_info.get('seq_tech', '10x scRNA-seq')
    analyst = db_info.get('analyst', '')
    notes = db_info.get('notes', '')
    
    # 解析输出目录中的图片和表格
    outdir = os.path.abspath(args.outdir)
    # 计算对应的对外挂载 delivery 目录
    upload_dir = outdir.replace("/output", "/upload") if outdir.endswith("/output") else outdir
    
    def find_file(pattern):
        # 兼容前缀，允许更灵活匹配结构
        matches = glob.glob(os.path.join(outdir, pattern), recursive=True)
        if not matches:
            matches = glob.glob(os.path.join(outdir, '**', pattern), recursive=True)
            
        if matches:
            real_path = matches[0]
            # 返回逻辑 upload 路径给数据库（确保WebUI上的链接可用）以及物理路径供拷贝使用
            upload_path = real_path.replace(outdir, upload_dir) if outdir.endswith("/output") else real_path
            return real_path, upload_path
        return "", ""
        
    # 优先找 Keep，如果没有再找 All (放宽匹配前缀以兼容带有 1_ 2_ 等前缀的老版本目录)
    umap_real, umap_upload = find_file('**/*annotation/*CellType_Keep/*.UMAP-blank.png')
    if not umap_real:
        umap_real, umap_upload = find_file('*annotation/*CellType_All/*.UMAP-blank.png')

    # 细胞占比图
    fraction_real, fraction_upload = find_file('**/*celltype_fraction/*Keep.in.Sample-bar.png')
    if not fraction_real:
        fraction_real, fraction_upload = find_file('*celltype_fraction/*All.in.Sample-bar.png')

    marker_table_real, marker_table_upload = find_file('**/*annotation/*CellType_All/*Cluster-CellType.anno.xls')
    anno_png_real, anno_png_upload = find_file('*annotation/*CellType_Keep/*CellType.dotplot.png')

    # HTML 报告 (在 pipeline 中最终会被重命名跑在 outdir 根目录)
    report_real, report_upload = find_file('*CellType.Annotation_report.html')
    
    # === 开始独立图床静态文件备份复制过程 ===
    webui_dir = os.path.dirname(os.path.abspath(db_path))
    assets_dir = os.path.join(webui_dir, 'assets', f"{crm_id}_{annotation_level}")
    os.makedirs(assets_dir, exist_ok=True)
    
    bk_umap = ""
    bk_fraction = ""
    bk_anno_png = ""
    bk_report = ""

    if umap_real and os.path.exists(umap_real):
        bk_umap = os.path.join(assets_dir, 'UMAP_preview.png')
        shutil.copy2(umap_real, bk_umap)
    else:
        raise FileNotFoundError(f"Error: 无法找到或拷贝 UMAP 图片到 assets，缺少源文件! 检查目录: {outdir}")
        
    if fraction_real and os.path.exists(fraction_real):
        bk_fraction = os.path.join(assets_dir, 'Fraction_bar.png')
        shutil.copy2(fraction_real, bk_fraction)
    else:
        raise FileNotFoundError(f"Error: 无法找到或拷贝细胞占比图到 assets，缺少源文件! 检查目录: {outdir}")
        
    if anno_png_real and os.path.exists(anno_png_real):
        bk_anno_png = os.path.join(assets_dir, 'DotPlot.png')
        shutil.copy2(anno_png_real, bk_anno_png)
    else:
        raise FileNotFoundError(f"Error: 无法找到或拷贝 DotPlot 气泡图到 assets，缺少源文件! 检查目录: {outdir}")
        
    if report_real and os.path.exists(report_real):
        bk_report = os.path.join(assets_dir, 'report.html')
        shutil.copy2(report_real, bk_report)
    else:
        raise FileNotFoundError(f"Error: 无法找到或拷贝 HTML 报告到 assets，缺少源文件! 检查目录: {outdir}")
    print(f"[*] Assests successfully backed up to: {assets_dir}")

    print(f"[*] Initializing database at {db_path}...")
    conn = init_db(db_path)
    cursor = conn.cursor()

    # 清除旧有数据避免同一次项目同一分群层次的重复运行堆叠
    if crm_id and annotation_level:
        print(f"[*] Cleaning up previous records for crm_id: {crm_id} and annotation_level: {annotation_level}...")
        cursor.execute("DELETE FROM projects WHERE crm_id = ? AND annotation_level = ?", (crm_id, annotation_level))
        cursor.execute("DELETE FROM marker_experience WHERE crm_id = ? AND annotation_level = ?", (crm_id, annotation_level))

    # 1. 写入主记录
    print("[*] Inserting project result record...")
    cursor.execute('''
        INSERT INTO projects (
            std_project_id, custom_project_id, crm_id, annotation_level, species, species_zh, tissue_type, tissue_type_zh, disease_type, seq_tech, analyst, notes,
            outdir, umap_path, fraction_path, marker_table_path, anno_png_path, report_path, bk_umap, bk_fraction, bk_anno_png, bk_report
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        std_project_id, custom_project_id, crm_id, annotation_level, species, species_zh, tissue_type, tissue_type_zh, disease_type, seq_tech, analyst, notes,
        outdir, umap_upload, fraction_upload, marker_table_upload, anno_png_upload, report_upload, bk_umap, bk_fraction, bk_anno_png, bk_report
    ))
    
    # 2. 写入经验库 (解析 marker table)
    if marker_table_real and os.path.exists(marker_table_real):
        print(f"[*] Parsing marker table {marker_table_real} for experience database...")
        try:
            with open(marker_table_real, 'r', encoding='utf-8') as mf:
                reader = csv.DictReader(mf, delimiter='\t')
                
                # Check for required fields
                if reader.fieldnames and all(col in reader.fieldnames for col in ('Cluster', 'CellType', 'Markers')):
                    records = []
                    for row in reader:
                        cluster = str(row['Cluster'])
                        cell_type = str(row['CellType'])
                        markers = str(row['Markers'])
                        reference = str(row.get('reference', ''))
                        description = str(row.get('description', ''))
                        records.append((custom_project_id, crm_id, annotation_level, species, tissue_type, cluster, cell_type, markers, reference, description))
                    
                    cursor.executemany('''
                        INSERT INTO marker_experience (custom_project_id, crm_id, annotation_level, species, tissue_type, cluster, cell_type, markers, reference, description)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', records)
                    print(f"[+] Successfully registered {len(records)} marker experiences.")
                else:
                    print("[-] Marker table is missing required columns (Cluster, CellType, Markers).")
        except Exception as e:
            print(f"[-] Failed to parse marker table: {e}")

    conn.commit()
    conn.close()
    print("[*] Database update finished successfully.")

if __name__ == '__main__':
    main()
