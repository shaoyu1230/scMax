import yaml
import argparse
import os
from datetime import datetime

def get_base_params(config):
    global_outdir = config.get("global_outdir", "./scMax_Out")
    rscript_cmd = config.get("Rscript_path", "Rscript")
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return global_outdir, rscript_cmd, base_dir

def build_bash_header(config_file, global_outdir):
    script_content = "#!/bin/bash\n\n"
    script_content += "set -e\n\n"
    script_content += f"echo 'Starting scMax Pipeline with origin config: {config_file}'\n"
    script_content += f"OUTDIR={global_outdir}\n"
    script_content += f"mkdir -p $OUTDIR\n\n"
    return script_content

def generate_qc_cmd(config, rscript_cmd, base_dir, next_input_rds):
    qc_conf = config.get("01_scQC") or config.get("scQC") or {}
    if not qc_conf.get("run", False):
        return "# Skip 01_scQC\n\n", next_input_rds
        
    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 1: 01_scQC (数据整合与质控)\n"
    script_content += "#" + "="*40 + "\n"
    script_dir = os.path.join(base_dir, "01_scQC", "script")
    qc_out = os.path.join("$OUTDIR", "01_scQC_out")
    script_content += f"if [ -d {qc_out} ]; then rm -rf {qc_out}; fi\n"
    script_content += f"mkdir -p {qc_out}\n"
    
    samples_file = qc_conf.get("samples_file", "")
    method = qc_conf.get("method", "merge")
    rdsfile = qc_conf.get("rdsfile", "")
    
    meta = qc_conf.get("metadatafile", "")
    col_sample = qc_conf.get("col_sample", "Sample")
    col_group = qc_conf.get("col_group", "Group")
    species = qc_conf.get("species", "mouse")
    doub_rate = qc_conf.get("doub_rate", 0.076)
    mt_genes = qc_conf.get("mt_genes", "")
    hb_genes = qc_conf.get("hb_genes", "")
    
    cmd = f"{rscript_cmd} {script_dir}/scQC.R \\\n"
    if samples_file:
        cmd += f"  --samples \"{samples_file}\" \\\n"
        cmd += f"  --method \"{method}\" \\\n"
    elif rdsfile:
        cmd += f"  --rds \"{rdsfile}\" \\\n"
    else:
        return "echo 'Error:必须指定 samples_file 或 rdsfile！'\nexit 1\n", ""
        
    if meta:
        cmd += f"  --meta \"{meta}\" \\\n"
    cmd += f"  --col_sample \"{col_sample}\" \\\n"
    cmd += f"  --col_group \"{col_group}\" \\\n"
    cmd += f"  --species \"{species}\" \\\n"
    cmd += f"  --doub_rate {doub_rate} \\\n"
    if mt_genes:
        cmd += f"  --mt_genes \"{mt_genes}\" \\\n"
    if hb_genes:
        cmd += f"  --hb_genes \"{hb_genes}\" \\\n"
    if qc_conf.get("run_doublet", False):
        cmd += f"  --run_doublet \\\n"
    if qc_conf.get("run_decontx", False):
        cmd += f"  --run_decontx \\\n"
    cmd += f"  --outdir {qc_out}\n\n"
    script_content += cmd
    return script_content, os.path.join(qc_out, "final_obj.rds")

def generate_filter_cmd(config, rscript_cmd, base_dir, next_input_rds):
    filter_conf = config.get("02_filter") or {}
    if not filter_conf.get("run", False):
        return "# Skip 02_filter\n\n", next_input_rds

    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 2: 02_filter (根据极值和质控结果对死细胞修剪脱落)\n"
    script_content += "#" + "="*40 + "\n"
    script_dir = os.path.join(base_dir, "02_filter", "script")
    filter_out = os.path.join("$OUTDIR", "02_filter_out")
    script_content += f"if [ -d {filter_out} ]; then rm -rf {filter_out}; fi\n"
    script_content += f"mkdir -p {filter_out}\n"
    
    col_sample = filter_conf.get("col_sample", "Sample")
    min_feat = filter_conf.get("min_feat", 200)
    max_feat = filter_conf.get("max_feat", 1000000000)
    min_count = filter_conf.get("min_count", 200)
    max_count = filter_conf.get("max_count", 1000000000)
    max_mt = filter_conf.get("max_mt", 20.0)
    max_hb = filter_conf.get("max_hb", 5.0)
    
    user_rds = filter_conf.get("rdsfile", "")
    if user_rds:
        next_input_rds = user_rds
        
    cmd = f"{rscript_cmd} {script_dir}/scFilter.R \\\n"
    cmd += f"  --rds \"{next_input_rds}\" \\\n"
    cmd += f"  --col_sample \"{col_sample}\" \\\n"
    cmd += f"  --min_feat {min_feat} --max_feat {max_feat} \\\n"
    cmd += f"  --min_count {min_count} --max_count {max_count} \\\n"
    cmd += f"  --max_mt {max_mt} --max_hb {max_hb} \\\n"
    if filter_conf.get("rm_doublet", False):
        cmd += f"  --rm_doublet \\\n"
    max_decontX = filter_conf.get("max_decontX", 1.0)
    if max_decontX < 1.0:
        cmd += f"  --max_decontX {max_decontX} \\\n"
        
    cmd += f"  --outdir {filter_out}\n\n"
    script_content += cmd
    return script_content, os.path.join(filter_out, "final_obj.rds")

def generate_cluster_cmd(config, rscript_cmd, base_dir, next_input_rds):
    cluster_conf = config.get("03_cluster") or {}
    if not cluster_conf.get("run", False):
        return "# Skip 03_cluster\n\n", next_input_rds

    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 3: 03_cluster (去批次整合与多策略+多分辨率聚类扫参)\n"
    script_content += "#" + "="*40 + "\n"
    script_dir = os.path.join(base_dir, "03_cluster", "script")
    cluster_out = os.path.join("$OUTDIR", "03_cluster_out")
    script_content += f"if [ -d {cluster_out} ]; then rm -rf {cluster_out}; fi\n"
    script_content += f"mkdir -p {cluster_out}\n"
    
    user_rds = cluster_conf.get("rdsfile", "")
    if user_rds:
        next_input_rds = user_rds
        
    subset_col = cluster_conf.get("subset_col", "")
    subset_val = cluster_conf.get("subset_val", "")
    methods = cluster_conf.get("methods", "harmony,rpca")
    resolutions = str(cluster_conf.get("resolutions", "0.2,0.4,0.6"))
    col_sample = cluster_conf.get("col_sample", "Sample")
    col_group = cluster_conf.get("col_group", "Group")
    integrate_by = cluster_conf.get("integrate_by", col_sample)
    nfeatures = cluster_conf.get("nfeatures", 2000)
    npcs = cluster_conf.get("npcs", 30)
    do_cluster = cluster_conf.get("do_cluster", False)
    do_refmarker = cluster_conf.get("do_refmarker", False)
    refmarker_file = cluster_conf.get("refmarker_file", "")
    
    db_conf = config.get("Database_Info", {})
    species = db_conf.get("species", "Human")
    if species.lower() == "mouse":
        auto_orgdb = "org.Mm.eg.db"
        auto_kegg = "mmu"
    else:
        auto_orgdb = "org.Hs.eg.db"
        auto_kegg = "hsa"
    
    cmd = f"{rscript_cmd} {script_dir}/scCluster.R \\\n"
    cmd += f"  --rds \"{next_input_rds}\" \\\n"
    cmd += f"  --outdir \"{cluster_out}\" \\\n"
    if subset_col and subset_val:
        cmd += f"  --subset_col \"{subset_col}\" --subset_val \"{subset_val}\" \\\n"
    cmd += f"  --methods \"{methods}\" \\\n"
    cmd += f"  --resolutions \"{resolutions}\" \\\n"
    cmd += f"  --col_sample \"{col_sample}\" \\\n"
    cmd += f"  --col_group \"{col_group}\" \\\n"
    cmd += f"  --integrate_by \"{integrate_by}\" \\\n"
    if do_cluster:
        cmd += f"  --do_cluster \\\n"
        cmd += f"  --orgdb \"{auto_orgdb}\" \\\n"
        cmd += f"  --organism_kegg \"{auto_kegg}\" \\\n"
    if do_refmarker:
        cmd += f"  --do_refmarker \\\n"
    if refmarker_file:
        cmd += f"  --refmarker_file \"{refmarker_file}\" \\\n"
    cmd += f"  --nfeatures {nfeatures} \\\n"
    cmd += f"  --npcs {npcs}\n\n"
    script_content += cmd
    next_res = os.path.join(cluster_out, methods.split(',')[0], f"res{resolutions.split(',')[0]}", f"{methods.split(',')[0]}_integrate_clustered.rds")
    return script_content, next_res

def generate_subcluster_cmd(config, rscript_cmd, base_dir, next_input_rds):
    subcluster_conf = config.get("04_subcluster") or {}
    if not subcluster_conf.get("run", False):
        return "# Skip 04_subcluster\n\n", next_input_rds

    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 4: 04_subcluster (高级亚群细分探针与特化整合策略)\n"
    script_content += "#" + "="*40 + "\n"
    script_dir = os.path.join(base_dir, "04_subcluster", "script")
    
    user_rds = subcluster_conf.get("rdsfile", "")
    if user_rds:
        next_input_rds = user_rds
        
    subset_col = subcluster_conf.get("subset_col", "res0.4")
    subset_val = str(subcluster_conf.get("subset_val", ""))
    safe_val_name = subset_val.replace(",", "_").replace(" ", "") if subset_val else "unknown"
    subcluster_out = os.path.join("$OUTDIR", f"04_subcluster_{safe_val_name}")
    
    script_content += f"if [ -d {subcluster_out} ]; then rm -rf {subcluster_out}; fi\n"
    script_content += f"mkdir -p {subcluster_out}\n"
    
    methods = subcluster_conf.get("methods", "re_harmony,original_integrated")
    resolutions = str(subcluster_conf.get("resolutions", "0.2,0.4"))
    col_sample = subcluster_conf.get("col_sample", "Sample")
    col_group = subcluster_conf.get("col_group", "Group")
    integrate_by = subcluster_conf.get("integrate_by", col_sample)
    nfeatures = subcluster_conf.get("nfeatures", 1500)
    npcs = subcluster_conf.get("npcs", 20)
    do_cluster = subcluster_conf.get("do_cluster", False)
    do_refmarker = subcluster_conf.get("do_refmarker", False)
    refmarker_file = subcluster_conf.get("refmarker_file", "")
    
    db_conf = config.get("Database_Info", {})
    species = db_conf.get("species", "Human")
    if species.lower() == "mouse":
        auto_orgdb = "org.Mm.eg.db"
        auto_kegg = "mmu"
    else:
        auto_orgdb = "org.Hs.eg.db"
        auto_kegg = "hsa"
    
    if not subset_col or not subset_val:
        return "echo 'Error: 启用 04_subcluster 但未填写必需指定的提取列或值！'\nexit 1\n", ""
    
    cmd = f"{rscript_cmd} {script_dir}/scSubCluster.R \\\n"
    cmd += f"  --rds \"{next_input_rds}\" \\\n"
    cmd += f"  --outdir \"{subcluster_out}\" \\\n"
    cmd += f"  --subset_col \"{subset_col}\" --subset_val \"{subset_val}\" \\\n"
    cmd += f"  --methods \"{methods}\" \\\n"
    cmd += f"  --resolutions \"{resolutions}\" \\\n"
    cmd += f"  --col_sample \"{col_sample}\" \\\n"
    cmd += f"  --col_group \"{col_group}\" \\\n"
    cmd += f"  --integrate_by \"{integrate_by}\" \\\n"
    if do_cluster:
        cmd += f"  --do_cluster \\\n"
        cmd += f"  --orgdb \"{auto_orgdb}\" \\\n"
        cmd += f"  --organism_kegg \"{auto_kegg}\" \\\n"
    if do_refmarker:
        cmd += f"  --do_refmarker \\\n"
    if refmarker_file:
        cmd += f"  --refmarker_file \"{refmarker_file}\" \\\n"
    cmd += f"  --nfeatures {nfeatures} \\\n"
    cmd += f"  --npcs {npcs}\n\n"
    script_content += cmd
    next_res = os.path.join(subcluster_out, methods.split(',')[0], f"res{resolutions.split(',')[0]}", f"{methods.split(',')[0]}_integrate_clustered.rds")
    return script_content, next_res

def generate_differential_cmd(config, rscript_cmd, base_dir, next_input_rds):
    diff_conf = config.get("06_differential") or {}
    if not diff_conf.get("run", False):
        return "# Skip 06_differential\n\n", next_input_rds

    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 6: 06_differential (差异表达分析挖掘模式 A/B)\n"
    script_content += "#" + "="*40 + "\n"
    script_dir = os.path.join(base_dir, "06_differential", "script")
    
    user_rds = diff_conf.get("rdsfile", "")
    if user_rds:
        next_input_rds = user_rds
    
    method = diff_conf.get("method", "all_markers")
    groupby = diff_conf.get("groupby", "CellType")
    split_by = diff_conf.get("split_by", "CellType")
    cmp_file = diff_conf.get("cmp_file", "")
    do_enrich = diff_conf.get("do_enrich", True)
    
    diff_out = os.path.join("$OUTDIR", f"06_differential_{method}")
    script_content += f"if [ -d {diff_out} ]; then rm -rf {diff_out}; fi\n"
    script_content += f"mkdir -p {diff_out}\n"
    
    db_conf = config.get("Database_Info", {})
    species = db_conf.get("species", "Human")
    if species.lower() == "mouse":
        auto_orgdb = "org.Mm.eg.db"
        auto_kegg = "mmu"
    else:
        auto_orgdb = "org.Hs.eg.db"
        auto_kegg = "hsa"
        
    cmd = f"{rscript_cmd} {script_dir}/scDifferential.R \\\n"
    cmd += f"  --rds \"{next_input_rds}\" \\\n"
    cmd += f"  --outdir \"{diff_out}\" \\\n"
    cmd += f"  --method \"{method}\" \\\n"
    cmd += f"  --groupby \"{groupby}\" \\\n"
    cmd += f"  --split_by \"{split_by}\" \\\n"
    if cmp_file:
        cmd += f"  --cmp_file \"{cmp_file}\" \\\n"
    if do_enrich:
        cmd += f"  --do_enrich \\\n"
    cmd += f"  --orgdb \"{auto_orgdb}\" \\\n"
    cmd += f"  --organism_kegg \"{auto_kegg}\"\n\n"
    
    script_content += cmd
    return script_content, next_input_rds

def generate_celltype_cmd(config, rscript_cmd, base_dir, next_input_rds, config_file):
    celltype_conf = config.get("05_celltype") or {}
    if not celltype_conf.get("run", False):
        return "# Skip 05_celltype\n\n", next_input_rds, ""

    script_content = "#" + "="*40 + "\n"
    script_content += "# Step 5: 05_celltype (专家组级人工标注、格式多输出与大屏构建)\n"
    script_content += "#" + "="*40 + "\n"
    db_conf = config.get("Database_Info", {})
    anno_level = db_conf.get("annotation_level", "MajorType").replace(" ", "_").replace("-", "_")
    
    script_dir = os.path.join(base_dir, "05_celltype", "script")
    celltype_base = os.path.join("$OUTDIR", f"05_annotation_{anno_level}")
    celltype_out = os.path.join(celltype_base, "output")
    celltype_upload = os.path.join(celltype_base, "upload")
    
    if celltype_conf.get("force_clean", False):
        script_content += f"if [ -d {celltype_base} ]; then rm -rf {celltype_base}; fi\n"
        
    script_content += f"mkdir -p {celltype_out}\n"
    script_content += f"mkdir -p {celltype_upload}\n"
    
    user_rds = celltype_conf.get("rdsfile", "")
    if user_rds:
        next_input_rds = user_rds
        
    annofile = celltype_conf.get("annofile", "")
    louper_path = celltype_conf.get("louper_path", "")
    
    do_cluster = celltype_conf.get("do_cluster", True)
    do_refmarker = celltype_conf.get("do_refmarker", True)
    do_annotation = celltype_conf.get("do_annotation", True)
    do_celltype = celltype_conf.get("do_celltype", True)
    do_celltype_de = celltype_conf.get("do_celltype_de", False)
    do_report = celltype_conf.get("do_report", True)
    
    if do_annotation and not annofile:
        return "echo 'Error: 启用 do_annotation 注释阶段但未传入映射注释表(annofile)！'\nexit 1\n", "", ""
        
    rds_dir = os.path.dirname(next_input_rds) if next_input_rds else ""
    dyn_args = f"--outdir {celltype_out}"
    rscript_skip_flags = ""
    
    if do_cluster:
        script_content += f"""
# 智能跳过: 如果数据源的上一级（03/04）已经做了 cluster_characterization，则直接软链，不再重复执行
if [ -d "{rds_dir}/cluster_characterization" ] && [ -n "{rds_dir}" ]; then
    echo "[Smart-Link] Found pre-computed cluster_characterization in {rds_dir}, linking directly..."
    ln -sfn "$(cd "{rds_dir}/cluster_characterization" && pwd)" "{celltype_out}/cluster_characterization"
    SKIP_CLUSTER_FLAG=""
    RSCRIPT_SKIP_CLUSTER="--skip_cluster"
else
    SKIP_CLUSTER_FLAG="--do_cluster"
    RSCRIPT_SKIP_CLUSTER=""
fi
"""
        dyn_args += " ${SKIP_CLUSTER_FLAG}"
        rscript_skip_flags += " ${RSCRIPT_SKIP_CLUSTER}"
        
    if do_refmarker:
        script_content += f"""
# 智能跳过: 如果数据源的上一级做了 marker_expression，则软链
if [ -d "{rds_dir}/marker_expression" ] && [ -n "{rds_dir}" ]; then
    echo "[Smart-Link] Found pre-computed marker_expression in {rds_dir}, linking directly..."
    ln -sfn "$(cd "{rds_dir}/marker_expression" && pwd)" "{celltype_out}/marker_expression"
    SKIP_REFMARKER_FLAG=""
    RSCRIPT_SKIP_REFMARKER="--skip_refmarker"
else
    SKIP_REFMARKER_FLAG="--do_refmarker"
    RSCRIPT_SKIP_REFMARKER=""
fi
"""
        dyn_args += " ${SKIP_REFMARKER_FLAG}"
        rscript_skip_flags += " ${RSCRIPT_SKIP_REFMARKER}"

    cmd = f"{rscript_cmd} {script_dir}/scCellType.R --config '{config_file}' --outdir \"{celltype_out}\" --inputrds \"{next_input_rds}\"{rscript_skip_flags}\n\n"
    
    python_path = celltype_conf.get("python_path", "python3")
    jupyter_path = celltype_conf.get("jupyter_path", "jupyter")
    
    if do_annotation: dyn_args += " --do_annotation"
    if do_celltype: dyn_args += " --do_celltype"
    if do_celltype_de: dyn_args += " --do_celltype_de"
        
    col_group_check = celltype_conf.get("col_group", "").strip()
    if not col_group_check or col_group_check.lower() == "none" or col_group_check == "NA":
        dyn_args += " --no_group"
    
    if do_report:
        cmd += "# >>> do_report 为 true，确认最终细胞类型注释完毕，组装出结题动态可视化报告及转化其他各平台交付格式(h5ad/cloupe) <<<\n"
        cmd += f"{python_path} {script_dir}/generate_dynamic_report.py {dyn_args}\n"
        cmd += f"cd {celltype_out} && {jupyter_path} nbconvert --no-input --template pj --HTMLExporter.embed_images=True --to html report_custom.ipynb && mv report_custom.html CellType.Annotation_report.html\n\n"
        
        cmd += f"{rscript_cmd} -e \"suppressMessages(library(SeuratDisk)); data <- readRDS('{celltype_out}/Rdata/Data-Annotation_CellType.rds'); SaveH5Seurat(data, filename = '{celltype_out}/Data_CellAnnotated.h5Seurat', overwrite=TRUE); Convert('{celltype_out}/Data_CellAnnotated.h5Seurat', dest = 'h5ad', overwrite=TRUE); unlink('{celltype_out}/Data_CellAnnotated.h5Seurat')\"\n"
        louper_arg = f" --louper_path \"{louper_path}\"" if louper_path else ""
        cmd += f"{rscript_cmd} {script_dir}/rds2cloupe.R -i \"{celltype_out}/Rdata/Data-Annotation_CellType.rds\" -o \"{celltype_out}\" -n Data_CellAnnotated_Cloupe{louper_arg}\n\n"
    else:
        cmd += "# >>> do_report 为 false，代表注释处于检查确认的中间态，在此跳过耗时的动态图文 HTML 报告渲染及 h5ad/cloupe 格式转换步骤 <<<\n\n"
    
    cmd += f"(\ncd {celltype_upload}\n"
    cmd += f"ln -sf ../output/cluster_characterization .\n"
    cmd += f"ln -sf ../output/marker_expression .\n"
    cmd += f"ln -sf ../output/annotation .\n"
    cmd += f"ln -sf ../output/celltype_fraction .\n"
    cmd += f"ln -sf ../output/celltype_characterization .\n"
    cmd += f"ln -sf ../output/CellType.Annotation_report.html .\n"
    cmd += f"ln -sf ../output/Rdata/Data-Annotation_CellType.rds .\n"
    readme_content = (
        "├── Data-Annotation_CellType.rds       分析结果的Seurat对象，需要用 R 语言读取\n"
        "├── cluster_characterization         cluster特征：包括cluster分布，样本/分组分布，cluster占比等\n"
        "├── marker_expression                参考marker表达绘图\n"
        "├── annotation                       细胞类型注释结果\n"
        "│   ├── CellType_All                   所有细胞注释结果\n"
        "│   └── CellType_Keep                  去除低质量、doublet等不感兴趣的细胞后，细胞注释结果\n"
        "├── celltype_fraction                细胞占比分析\n"
        "├── celltype_characterization        细胞类型特征，包括差异分析和富集分析\n"
        "└── CellType.Annotation_report.html    细胞类型注释报告\n"
    )
    cmd += f"cat << 'EOF' > 结果文件说明.txt\n{readme_content}EOF\n"
    cmd += ")\n\n"
    script_content += cmd
    return script_content, os.path.join(celltype_out, "Rdata", "Data-Annotation_CellType.rds"), celltype_out

def generate_db_cmd(config, base_dir, celltype_out, config_file):
    db_conf = config.get("Database_Info", {})
    if not db_conf.get("run", False):
        return ""
    script_content = f'''
# ========================================
# Execute DB Registration (Central UI & Experience DB)
# ========================================
python3 {os.path.join(base_dir, "scMax_db_manager.py")} -c '{os.path.abspath(config_file)}' -o "{celltype_out}"
'''
    return script_content

def generate_bash_script_from_dict(config, script_outdir, config_file="base_config.yaml", run_steps=None):
    os.makedirs(script_outdir, exist_ok=True)
    global_outdir, rscript_cmd, base_dir = get_base_params(config)
    os.makedirs(global_outdir, exist_ok=True)
    
    script_content = build_bash_header(config_file, global_outdir)
    
    next_input_rds = ""
    
    step_content, next_input_rds = generate_qc_cmd(config, rscript_cmd, base_dir, next_input_rds)
    script_content += step_content
    
    step_content, next_input_rds = generate_filter_cmd(config, rscript_cmd, base_dir, next_input_rds)
    script_content += step_content
    
    step_content, next_input_rds = generate_cluster_cmd(config, rscript_cmd, base_dir, next_input_rds)
    script_content += step_content
    
    step_content, next_input_rds = generate_subcluster_cmd(config, rscript_cmd, base_dir, next_input_rds)
    script_content += step_content
    
    step_content, next_input_rds, celltype_out = generate_celltype_cmd(config, rscript_cmd, base_dir, next_input_rds, config_file)
    script_content += step_content
    
    step_content, next_input_rds = generate_differential_cmd(config, rscript_cmd, base_dir, next_input_rds)
    script_content += step_content
    
    if celltype_out:
        script_content += generate_db_cmd(config, base_dir, celltype_out, config_file)
        
    script_content += "\necho 'scMax Pipeline completed successfully!'\n"
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    step_suffix = "_" + "_".join(run_steps) if run_steps else ""
    sh_file = os.path.join(script_outdir, f"run_scMax_pipeline{step_suffix}_{timestamp}.sh")
    with open(sh_file, "w", encoding='utf-8') as f:
        f.write(script_content)
    os.chmod(sh_file, 0o755)
    return sh_file, global_outdir

def generate_bash_script(config_file, script_outdir):
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
        
    run_steps = []
    if ("01_scQC" in config or "scQC" in config) and config.get("01_scQC", config.get("scQC", {})).get("run", True):
        run_steps.append("01")
    if "02_filter" in config and config.get("02_filter", {}).get("run", True):
        run_steps.append("02")
    if "03_cluster" in config and config.get("03_cluster", {}).get("run", True):
        run_steps.append("03")
    if "04_subcluster" in config and config.get("04_subcluster", {}).get("run", False):
        run_steps.append("04")
    
    if "06_differential" in config and config.get("06_differential", {}).get("run", False):
        run_steps.append("06")
    
    celltype_conf = config.get("05_celltype", {})
    if "05_celltype" in config and celltype_conf.get("run", False):
        sub_opts = []
        if celltype_conf.get("do_cluster", True): sub_opts.append("clus")
        if celltype_conf.get("do_refmarker", True): sub_opts.append("refm")
        if celltype_conf.get("do_annotation", True): sub_opts.append("anno")
        if celltype_conf.get("do_celltype", True): sub_opts.append("ctype")
        if celltype_conf.get("do_celltype_de", False): sub_opts.append("de")
        if sub_opts:
            run_steps.append("05_" + "-".join(sub_opts))
        else:
            run_steps.append("05")
            
    sh_file, global_outdir = generate_bash_script_from_dict(config, script_outdir, config_file, run_steps)
    print(f"生成的 bash 分析脚本已保存在投递口: {sh_file}")
    print(f"后台数据的实际归档将送往: {global_outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="生成 scMax 分析流程的 bash 脚本")
    parser.add_argument("-c", "--config", required=True, help="YAML 配置文件")
    parser.add_argument("-o", "--outdir", default=".", help="生成的 bash 运行脚本所存放的目录位置 (默认当前目录)")
    args = parser.parse_args()
    
    generate_bash_script(args.config, args.outdir)
