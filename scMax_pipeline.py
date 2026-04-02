import yaml
import argparse
import os

def generate_bash_script(config_file, script_outdir):
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    # 建立生成 bash 投递脚本的目录
    os.makedirs(script_outdir, exist_ok=True)
    
    # 从 config 中独取实际的分析落地主目录
    global_outdir = config.get("global_outdir", "./scMax_Out")
    os.makedirs(global_outdir, exist_ok=True)
    
    script_content = "#!/bin/bash\n\n"
    script_content += "set -e\n\n"
    script_content += f"echo 'Starting scMax Pipeline with config: {config_file}'\n"
    script_content += f"OUTDIR={global_outdir}\n"
    script_content += f"mkdir -p $OUTDIR\n\n"
    
    # 获取 Rscript 环境路径，默认使用系统环境里的 Rscript
    rscript_cmd = config.get("Rscript_path", "Rscript")
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 初始化全局入库变量，防止步骤跳过导致未定义
    global_final_umap = "NA"
    global_final_html = "NA"
    global_final_fraction = "NA"
    global_final_anno_table = "NA"
    
    # 贯穿全程的中间文件输入传递变量
    next_input_rds = ""
    
    # --- Step 1: 01_scQC ---
    script_content += "#" + "="*40 + "\n"
    script_content += "# Step 1: 01_scQC (数据整合与质控)\n"
    script_content += "#" + "="*40 + "\n"
    qc_in_config = "01_scQC" in config or "scQC" in config
    qc_conf = config.get("01_scQC") or config.get("scQC") or {}
    if qc_in_config and qc_conf.get("run", True):
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
            print("Error: 必须在配置文件中指定 samples_file 或者 rdsfile 作为输入数据来源！")
            return
            
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
        # 将下一步的默认输入覆盖为刚才的产出文件
        next_input_rds = os.path.join(qc_out, "final_obj.rds")
    else:
        script_content += "# Skip 01_scQC\n\n"

    # --- Step 2: 02_filter ---
    script_content += "#" + "="*40 + "\n"
    script_content += "# Step 2: 02_filter (根据极值和质控结果对死细胞修剪脱落)\n"
    script_content += "#" + "="*40 + "\n"
    filter_in_config = "02_filter" in config
    filter_conf = config.get("02_filter") or {}
    if filter_in_config and filter_conf.get("run", True):
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
        
        # 允许通过配置显式传递的独立 rds
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
        # 兼容了decontX容差数值指定
        max_decontX = filter_conf.get("max_decontX", 1.0)
        if max_decontX < 1.0:
            cmd += f"  --max_decontX {max_decontX} \\\n"
            
        cmd += f"  --outdir {filter_out}\n\n"
        script_content += cmd
        # 将下一步的默认输入覆盖为过滤后产出的文件
        next_input_rds = os.path.join(filter_out, "final_obj.rds")
    else:
        script_content += "# Skip 02_filter\n\n"
        
    # --- Step 3: 03_cluster ---
    script_content += "#" + "="*40 + "\n"
    script_content += "# Step 3: 03_cluster (去批次整合与多策略+多分辨率聚类扫参)\n"
    script_content += "#" + "="*40 + "\n"
    cluster_in_config = "03_cluster" in config
    cluster_conf = config.get("03_cluster") or {}
    if cluster_in_config and cluster_conf.get("run", True):
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
        integrate_by = cluster_conf.get("integrate_by", col_sample) # 默认使用 Sample 级别的列
        nfeatures = cluster_conf.get("nfeatures", 2000)
        npcs = cluster_conf.get("npcs", 30)
        do_cluster = cluster_conf.get("do_cluster", False)
        do_refmarker = cluster_conf.get("do_refmarker", False)
        refmarker_file = cluster_conf.get("refmarker_file", "")
        
        # 自动推断 orgdb 用于早期聚类的直接富集
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
        
    # --- Step 4: 04_subcluster ---
    script_content += "#" + "="*40 + "\n"
    script_content += "# Step 4: 04_subcluster (高级亚群细分探针与特化整合策略)\n"
    script_content += "#" + "="*40 + "\n"
    subcluster_in_config = "04_subcluster" in config
    subcluster_conf = config.get("04_subcluster") or {}
    if subcluster_in_config and subcluster_conf.get("run", False):
        script_dir = os.path.join(base_dir, "04_subcluster", "script")
        user_rds = subcluster_conf.get("rdsfile", "")
        # 如果亚群挖掘未指明提取的母体文件，则尝试延续上一环节最后输出的一份
        if user_rds:
            next_input_rds = user_rds
            
        subset_col = subcluster_conf.get("subset_col", "res0.4")
        subset_val = str(subcluster_conf.get("subset_val", ""))
        
        # 优化目录结构：带上有标识性的子群参数作为文件夹名称防覆盖
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
        
        # 强制亚群参数必须有界定
        if not subset_col or not subset_val:
            cmd = "echo 'Error: 启用 04_subcluster 但未在配置文件中填写必需指定的提取列(subset_col)或值(subset_val)！'\n"
            script_content += cmd
        else:
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
            
    else:
        script_content += "# Skip 04_subcluster\n\n"
        
    # --- Step 5: 05_celltype & 专家注释与大图产出 ---
    script_content += "#" + "="*40 + "\n"
    script_content += "# Step 5: 05_celltype (专家组级人工标注、格式多输出与大屏构建)\n"
    script_content += "#" + "="*40 + "\n"
    celltype_in_config = "05_celltype" in config
    celltype_conf = config.get("05_celltype") or {}
    if celltype_in_config and celltype_conf.get("run", False):
        db_conf = config.get("Database_Info", {})
        anno_level = db_conf.get("annotation_level", "MajorType").replace(" ", "_").replace("-", "_")
        
        script_dir = os.path.join(base_dir, "05_celltype", "script")
        celltype_base = os.path.join("$OUTDIR", f"05_annotation_{anno_level}")
        celltype_out = os.path.join(celltype_base, "output")
        celltype_upload = os.path.join(celltype_base, "upload")
        
        # 移除如果重新运行部分分析时整个结果目录被清空的预处理行为
        # 新增 force_clean 开关允许用户在第二次全新分析时自主选择完全清空目录
        if celltype_conf.get("force_clean", False):
            script_content += f"if [ -d {celltype_base} ]; then rm -rf {celltype_base}; fi\n"
            
        script_content += f"mkdir -p {celltype_out}\n"
        script_content += f"mkdir -p {celltype_upload}\n"
        
        user_rds = celltype_conf.get("rdsfile", "")
        # 如果未提供 rdsfile，默认读取 03 的结果
        if user_rds:
            next_input_rds = user_rds
            
        annofile = celltype_conf.get("annofile", "")
        louper_path = celltype_conf.get("louper_path", "")
        
        # === 可选参数提前读取 ===
        do_cluster = celltype_conf.get("do_cluster", True)
        do_refmarker = celltype_conf.get("do_refmarker", True)
        do_annotation = celltype_conf.get("do_annotation", True)
        do_celltype = celltype_conf.get("do_celltype", True)
        do_celltype_de = celltype_conf.get("do_celltype_de", False)
        
        if do_annotation and not annofile:
            script_content += "echo 'Error: 启用 do_annotation 注释阶段但未传入映射注释表(annofile)！'\n"
            script_content += "exit 1\n\n"
            
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

        # === 直接调用分析脚本（不依赖 Makefile），读取主 YAML 配置 ===
        cmd = f"{rscript_cmd} {script_dir}/scCellType.R --config '{config_file}' --outdir \"{celltype_out}\" --inputrds \"{next_input_rds}\"{rscript_skip_flags}\n\n"
        
        # === 动态组装 HTML 报告 ===
        python_path = celltype_conf.get("python_path", "python3")
        jupyter_path = celltype_conf.get("jupyter_path", "jupyter")
        
        if do_annotation: dyn_args += " --do_annotation"
        if do_celltype: dyn_args += " --do_celltype"
        if do_celltype_de: dyn_args += " --do_celltype_de"
            
        # 检测并传递无分组标记
        col_group_check = celltype_conf.get("col_group", "").strip()
        if not col_group_check or col_group_check.lower() == "none" or col_group_check == "NA":
            dyn_args += " --no_group"
        
        do_report = celltype_conf.get("do_report", True)
        if do_report:
            cmd += "# >>> do_report 为 true，确认最终细胞类型注释完毕，组装出结题动态可视化报告及转化其他各平台交付格式(h5ad/cloupe) <<<\n"
            cmd += f"{python_path} {script_dir}/generate_dynamic_report.py {dyn_args}\n"
            cmd += f"cd {celltype_out} && {jupyter_path} nbconvert --no-input --template pj --HTMLExporter.embed_images=True --to html report_custom.ipynb && mv report_custom.html CellType.Annotation_report.html\n\n"
            
            # Seurat 对象多格式输出：rds / cloupe / h5ad
            cmd += f"{rscript_cmd} -e \"suppressMessages(library(SeuratDisk)); data <- readRDS('{celltype_out}/Rdata/Data-Annotation_CellType.rds'); SaveH5Seurat(data, filename = '{celltype_out}/Data_CellAnnotated.h5Seurat', overwrite=TRUE); Convert('{celltype_out}/Data_CellAnnotated.h5Seurat', dest = 'h5ad', overwrite=TRUE); unlink('{celltype_out}/Data_CellAnnotated.h5Seurat')\"\n"
            louper_arg = f" --louper_path \"{louper_path}\"" if louper_path else ""
            cmd += f"{rscript_cmd} {script_dir}/rds2cloupe.R -i \"{celltype_out}/Rdata/Data-Annotation_CellType.rds\" -o \"{celltype_out}\" -n Data_CellAnnotated_Cloupe{louper_arg}\n\n"
        else:
            cmd += "# >>> do_report 为 false，代表注释处于检查确认的中间态，在此跳过耗时的动态图文 HTML 报告渲染及 h5ad/cloupe 格式转换步骤 <<<\n\n"
        
        # 构建 delivery 制品
        cmd += f"(\ncd {celltype_upload}\n"
        cmd += f"ln -sf ../output/cluster_characterization .\n"
        cmd += f"ln -sf ../output/marker_expression .\n"
        cmd += f"ln -sf ../output/annotation .\n"
        cmd += f"ln -sf ../output/celltype_fraction .\n"
        cmd += f"ln -sf ../output/celltype_characterization .\n"
        cmd += f"ln -sf ../output/CellType.Annotation_report.html .\n"
        # 交付rds结果
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
        # 记录最终挂载位置供入库
        global_final_umap = str(os.path.join(celltype_out, "annotation", "CellType_All", "1_CellType.UMAP.png"))
        global_final_html = str(os.path.join(celltype_out, "CellType.Annotation_report.html"))
        global_final_fraction = str(os.path.join(celltype_out, "celltype_fraction", "2_CellType_Keep.in.Sample-bar.png"))
        global_final_anno_table = str(os.path.join(celltype_out, "annotation", "CellType_All", "0_Cluster-CellType.anno.xls"))

    else:
        script_content += "# Skip 05_celltype\n\n"

    # --- Step 6: Database Integration ---
    db_conf = config.get("Database_Info", {})
    if db_conf.get("run", False):
        script_content += f'''
# ========================================
# Execute DB Registration (Central UI)
# ========================================
# ========================================
# Execute DB Registration (Central UI & Experience DB)
# ========================================
python3 {os.path.join(base_dir, "scMax_db_manager.py")} -c '{os.path.abspath(config_file)}' -o "{celltype_out}"

'''
    script_content += "\necho 'scMax Pipeline completed successfully!'\n"
    
    run_steps = []
    qc_in_config = "01_scQC" in config or "scQC" in config
    qc_conf = config.get("01_scQC") or config.get("scQC") or {}
    if qc_in_config and qc_conf.get("run", True):
        run_steps.append("01")
        
    filter_in_config = "02_filter" in config
    filter_conf = config.get("02_filter") or {}
    if filter_in_config and filter_conf.get("run", True):
        run_steps.append("02")
        
    cluster_in_config = "03_cluster" in config
    cluster_conf = config.get("03_cluster") or {}
    if cluster_in_config and cluster_conf.get("run", True):
        run_steps.append("03")
        
    subcluster_in_config = "04_subcluster" in config
    subcluster_conf = config.get("04_subcluster") or {}
    if subcluster_in_config and subcluster_conf.get("run", False):
        run_steps.append("04")
        
    celltype_in_config = "05_celltype" in config
    celltype_conf = config.get("05_celltype") or {}
    if celltype_in_config and celltype_conf.get("run", False):
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
        
    step_suffix = "_" + "_".join(run_steps) if run_steps else ""
    sh_file = os.path.join(script_outdir, f"run_scMax_pipeline{step_suffix}.sh")
    with open(sh_file, "w", encoding='utf-8') as f:
        f.write(script_content)
    os.chmod(sh_file, 0o755)
    
    print(f"生成的 bash 分析脚本已保存在投递口: {sh_file}")
    print(f"后台数据的实际归档将送往: {global_outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="生成 scMax 分析流程的 bash 脚本")
    parser.add_argument("-c", "--config", required=True, help="YAML 配置文件")
    parser.add_argument("-o", "--outdir", required=True, help="生成的 bash 运行脚本所存放的目录位置")
    args = parser.parse_args()
    
    generate_bash_script(args.config, args.outdir)
