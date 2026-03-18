#!/bin/bash

set -e

echo 'Starting scMax Pipeline with config: /Users/angela/Documents/00_annoroad/06_Development/scMax/config_template.yaml'
OUTDIR=/Users/angela/Documents/00_annoroad/06_Development/scMax/test_output_v2
mkdir -p $OUTDIR

#========================================
# Step 1: 01_scQC (数据整合与质控)
#========================================
if [ -d $OUTDIR/01_scQC_out ]; then rm -rf $OUTDIR/01_scQC_out; fi
mkdir -p $OUTDIR/01_scQC_out
Rscript /Users/angela/Documents/00_annoroad/06_Development/scMax/01_scQC/script/scQC.R \
  --samples 'samples_info.csv' \
  --method 'merge' \
  --meta 'matadata.csv' \
  --col_sample 'Sample' \
  --col_group 'Group' \
  --species 'mouse' \
  --doub_rate 0.076 \
  --outdir $OUTDIR/01_scQC_out

#========================================
# Step 2: 02_filter (根据极值和质控结果对死细胞修剪脱落)
#========================================
if [ -d $OUTDIR/02_filter_out ]; then rm -rf $OUTDIR/02_filter_out; fi
mkdir -p $OUTDIR/02_filter_out
Rscript /Users/angela/Documents/00_annoroad/06_Development/scMax/02_filter/script/scFilter.R \
  --rds '$OUTDIR/01_scQC_out/final_obj.rds' \
  --col_sample 'Sample' \
  --min_feat 200 --max_feat 100000 \
  --min_count 200 --max_count 1000000 \
  --max_mt 20.0 --max_hb 5.0 \
  --outdir $OUTDIR/02_filter_out

#========================================
# Step 3: 03_cluster (去批次整合与多策略+多分辨率聚类扫参)
#========================================
if [ -d $OUTDIR/03_cluster_out ]; then rm -rf $OUTDIR/03_cluster_out; fi
mkdir -p $OUTDIR/03_cluster_out
#========================================
# Step 4: 04_subcluster (高级亚群细分探针与特化整合策略)
#========================================
if [ -d $OUTDIR/04_subcluster_out ]; then rm -rf $OUTDIR/04_subcluster_out; fi
mkdir -p $OUTDIR/04_subcluster_out
Rscript /Users/angela/Documents/00_annoroad/06_Development/scMax/04_subcluster/script/scSubCluster.R \
  --rds '$OUTDIR/02_filter_out/final_obj.rds' \
  --outdir '$OUTDIR/04_subcluster_out' \
  --subset_col 'res0.4' --subset_val '1,3,5' \
  --methods 're_harmony,original_integrated' \
  --resolutions '0.2,0.4' \
  --col_sample 'Sample' \
  --col_group 'Group' \
  --integrate_by 'Sample' \
  --nfeatures 1500 \
  --npcs 20

#========================================
# Step 5: 05_celltype (专家组级人工标注、格式多输出与大屏构建)
#========================================
if [ -d $OUTDIR/05_celltype_out ]; then rm -rf $OUTDIR/05_celltype_out; fi
mkdir -p $OUTDIR/05_celltype_out
Rscript /Users/angela/Documents/00_annoroad/06_Development/scMax/05_celltype/script/scCellType.R --config '/Users/angela/Documents/00_annoroad/06_Development/scMax/config_template.yaml' --outdir $OUTDIR/05_celltype_out --inputrds '$OUTDIR/02_filter_out/final_obj.rds'

Rscript -e "suppressMessages(library(SeuratDisk)); data <- readRDS('$OUTDIR/05_celltype_out/Rdata/Data-Annotation_CellType.rds'); SaveH5Seurat(data, filename = '$OUTDIR/05_celltype_out/Data_CellAnnotated.h5Seurat', overwrite=TRUE); Convert('$OUTDIR/05_celltype_out/Data_CellAnnotated.h5Seurat', dest = 'h5ad', overwrite=TRUE); unlink('$OUTDIR/05_celltype_out/Data_CellAnnotated.h5Seurat')"
Rscript /Users/angela/Documents/00_annoroad/06_Development/scMax/05_celltype/script/rds2cloupe.R -i $OUTDIR/05_celltype_out/Rdata/Data-Annotation_CellType.rds -o $OUTDIR/05_celltype_out -n Data_CellAnnotated_Cloupe

# >>> 自动化组装动态可视化报告 <<<
python3 /Users/angela/Documents/00_annoroad/06_Development/scMax/05_celltype/script/generate_dynamic_report.py --outdir $OUTDIR/05_celltype_out --do_cluster --do_refmarker --do_annotation --do_celltype --no_group
cd $OUTDIR/05_celltype_out && jupyter nbconvert --no-input --template pj --to html report_custom.ipynb && mv report_custom.html CellType.Annotation_report.html


# ========================================
# Execute DB Registration (Central UI)
# ========================================
python3 -c "
import sqlite3, os
db_path = r'/Users/angela/Documents/00_annoroad/06_Development/scMax/WebUI/scMax_projects.db'
os.makedirs(os.path.dirname(db_path), exist_ok=True)
conn = sqlite3.connect(db_path)
c = conn.cursor()
c.execute('''
    CREATE TABLE IF NOT EXISTS projects (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        std_project_id TEXT,
        custom_project_id TEXT,
        species TEXT,
        result_path TEXT,
        umap_path TEXT,
        fraction_path TEXT,
        anno_table_path TEXT,
        anno_html_path TEXT,
        notes TEXT,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
''')
# Add column for older sqlite versions dynamically if needed (Optional fault tolerance)
try:
    c.execute('ALTER TABLE projects ADD COLUMN anno_html_path TEXT')
except Exception:
    pass
try:
    c.execute('ALTER TABLE projects ADD COLUMN fraction_path TEXT')
except Exception:
    pass
try:
    c.execute('ALTER TABLE projects ADD COLUMN anno_table_path TEXT')
except Exception:
    pass

c.execute('INSERT INTO projects (std_project_id, custom_project_id, species, result_path, umap_path, fraction_path, anno_table_path, anno_html_path, notes) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
    ('{std_id}', '{cus_id}', '{species}', r'{outdir}', '{umap_to_record}', '{fraction_to_record}', '{anno_table_to_record}', '{html_to_record}', '{notes}'))
conn.commit()
conn.close()
print('===> 极速投递：结果集及元数据已成功归档录入 Web UI 中央查询系统!')
"

echo 'scMax Pipeline completed successfully!'
