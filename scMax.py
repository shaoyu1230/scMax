#!/usr/bin/env python3
import argparse
import sys
import os
import yaml
import scMax_pipeline

def merge_configs(base, extra):
    """递归合并两个字典配置"""
    for key, value in extra.items():
        if isinstance(value, dict) and key in base and isinstance(base[key], dict):
            merge_configs(base[key], value)
        else:
            base[key] = value

def load_base_config(config_path):
    import copy
    # 1. Level 1: 内置极简兜底 (从 pipeline 模块读取)
    final_config = copy.deepcopy(scMax_pipeline.DEFAULT_CONFIG)
    
    # 2. Level 2: 全局环境配置 (scMax.config)
    # 优先在 scMax.py 所在目录及其当前运行目录寻找
    base_dir = os.path.dirname(os.path.abspath(__file__))
    global_conf_candidates = [
        os.path.join(base_dir, "scMax.config"),
        "scMax.config"
    ]
    for g_path in global_conf_candidates:
        if os.path.exists(g_path):
            try:
                with open(g_path, 'r', encoding='utf-8') as f:
                    g_conf = yaml.safe_load(f)
                    if g_conf:
                        merge_configs(final_config, g_conf)
                break # 找到第一个即停止
            except Exception as e:
                print(f"Warning: 加载全局配置 {g_path} 失败: {e}")

    # 3. Level 3: 项目业务配置 (通过 -c 指定)
    if not os.path.exists(config_path):
        if config_path == "base_config.yaml":
            # 如果是默认指定的 base_config 不存在，则依赖 L1/L2 继续
            return final_config
        print(f"Error: 找不到指定的项目配置文件: '{config_path}'")
        sys.exit(1)
        
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            user_config = yaml.safe_load(f)
            if user_config:
                merge_configs(final_config, user_config)
    except Exception as e:
        print(f"Error: 加载项目配置 {config_path} 失败: {e}")
        sys.exit(1)
                    
    return final_config

def run_step(step_name, args, step_conf):
    base_config = load_base_config(args.config)
    if args.outdir != ".":
        base_config['global_outdir'] = args.outdir

    # Disable all steps except the current one
    for s in ["01_scQC", "scQC", "02_filter", "03_cluster", "04_subcluster", "05_celltype", "06_merge_sub", "07_differential", "06_differential", "07_merge_sub", "merge_sub", "Database_Info"]:
        if s in base_config and isinstance(base_config[s], dict):
            if s != step_name and s != "Database_Info":
                base_config[s]['run'] = False
                
    # Update or create the specific step config
    step_conf['run'] = True
    if step_name not in base_config:
        base_config[step_name] = {}
    
    # 自动将 CLI 传入的 species 注入到 Database_Info (如果存在)
    if hasattr(args, 'species') and args.species:
        if 'Database_Info' not in base_config:
            base_config['Database_Info'] = {}
        base_config['Database_Info']['species'] = args.species
        
    base_config[step_name].update(step_conf)

    run_opts = [step_name.split("_")[0]]
    if step_name == "05_celltype":
        sub_opts = []
        if step_conf.get("do_cluster", True): sub_opts.append("clus")
        if step_conf.get("do_refmarker", True): sub_opts.append("refm")
        if step_conf.get("do_annotation", True): sub_opts.append("anno")
        if step_conf.get("do_celltype", True): sub_opts.append("ctype")
        if step_conf.get("do_celltype_de", False): sub_opts.append("de")
        
        if sub_opts:
            run_opts = ["05_" + "-".join(sub_opts)]

    sh_file, _ = scMax_pipeline.generate_bash_script_from_dict(
        config=base_config, 
        script_outdir=args.outdir, 
        config_file=args.config,
        run_steps=run_opts
    )
    print(f"-> 成功根据子命令 `{args.command}` 生成执行脚本：{sh_file}")
    print(f"您可以运行：sh {sh_file}")

def main():
    parser = argparse.ArgumentParser(
        description="scMax v2.0 - 高级单细胞流水线 (Command-Line Interface)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest="command", help="可选子阶段或多阶打包运行模式 (run_all)")

    # ---------------------------------------------------------
    # 01_scQC
    # ---------------------------------------------------------
    parser_qc = subparsers.add_parser("qc", help="数据整合与初步质控 (01_scQC)")
    parser_qc.add_argument("-c", "--config", default="base_config.yaml", help="基础环境变量配置 (Base Config)")
    parser_qc.add_argument("-o", "--outdir", default=".", help="脚本投递位置和输出前缀")
    
    group_qc = parser_qc.add_mutually_exclusive_group(required=True)
    group_qc.add_argument("--samples_file", help="样本信息与路径对照表 xls/txt")
    group_qc.add_argument("--rds", help="直接传入一个Seurat RDS对象进行质控")
    
    parser_qc.add_argument("--method", default="merge", help="数据整合模式 (merge 等)")
    parser_qc.add_argument("--meta", default="", help="额外的元信息对照表(可选)")
    parser_qc.add_argument("--col_sample", default="Sample", help="样本分群界定列名")
    parser_qc.add_argument("--col_group", default="Group", help="大组分群界定列名")
    parser_qc.add_argument("--species", default="mouse", help="内置物种 (mouse/human) 用于质控线")
    parser_qc.add_argument("--doub_rate", type=float, default=0.076, help="DoubletFinder 丢弃率")
    parser_qc.add_argument("--mt_genes", default="", help="线粒体基因开头，例：'MT-'")
    parser_qc.add_argument("--hb_genes", default="", help="红细胞基因开头，例：'HBA|HBB'")
    parser_qc.add_argument("--run_doublet", action="store_true", help="开启基于 DoubletFinder 去除双细胞")
    parser_qc.add_argument("--run_decontx", action="store_true", help="开启基于 decontX 去除环境污染背景")

    # ---------------------------------------------------------
    # 02_filter
    # ---------------------------------------------------------
    parser_filter = subparsers.add_parser("filter", help="硬过滤修剪死细胞 (02_filter)")
    parser_filter.add_argument("-c", "--config", default="base_config.yaml", help="基础配置")
    parser_filter.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_filter.add_argument("--rds", default="", help="来自上一级的 RDS 结果路径；留空则从配置文件读取")
    parser_filter.add_argument("--col_sample", default="Sample", help="样本列")
    parser_filter.add_argument("--min_feat", type=int, default=200, help="最小基因数")
    parser_filter.add_argument("--max_feat", type=int, default=1000000000, help="最大基因数")
    parser_filter.add_argument("--min_count", type=int, default=200, help="最小计数")
    parser_filter.add_argument("--max_count", type=int, default=1000000000, help="最大计数")
    parser_filter.add_argument("--max_mt", type=float, default=20.0, help="线粒体阈值百分比")
    parser_filter.add_argument("--max_hb", type=float, default=5.0, help="红细胞阈值百分比")
    parser_filter.add_argument("--rm_doublet", action="store_true", help="物理删除判定为Doublet的细胞")
    parser_filter.add_argument("--max_decontX", type=float, default=1.0, help="保留最大环境污染背景分（1.0代表不以此过滤）")

    # ---------------------------------------------------------
    # 03_cluster
    # ---------------------------------------------------------
    parser_cluster = subparsers.add_parser("cluster", help="降维聚类分析 (03_cluster)")
    parser_cluster.add_argument("-c", "--config", default="base_config.yaml", help="基础配置")
    parser_cluster.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_cluster.add_argument("--rds", default="", help="需要提取聚类的RDS对象路径；留空则从配置文件读取")
    parser_cluster.add_argument("--subset_col", default="", help="提前抽取的列 (可选)")
    parser_cluster.add_argument("--subset_val", default="", help="提前抽取的值 (可选)")
    parser_cluster.add_argument("--methods", default="harmony,rpca", help="去批次重整合的方法")
    parser_cluster.add_argument("--resolutions", default="0.2,0.4,0.6", help="探索的多个精细分别率")
    parser_cluster.add_argument("--col_sample", default="Sample", help="统计对照样本列")
    parser_cluster.add_argument("--col_group", default="Group", help="统计大组列")
    parser_cluster.add_argument("--integrate_by", default="Sample", help="执行批次整合依赖依据的列")
    parser_cluster.add_argument("--nfeatures", type=int, default=2000, help="高变基因数")
    parser_cluster.add_argument("--npcs", type=int, default=30, help="主成分数")
    parser_cluster.add_argument("--do_cluster", action="store_true", help="当场渲染高级特异性热图及执行全局GO/KEGG富集")
    parser_cluster.add_argument("--do_refmarker", action="store_true", help="根据指定 Marker 执行点图打靶")
    parser_cluster.add_argument("--refmarker_file", default="", help="带有靶点标记的 txt 文件路径")

    # ---------------------------------------------------------
    # 04_subcluster
    # ---------------------------------------------------------
    parser_subcluster = subparsers.add_parser("subcluster", help="高分辨亚群分离提取再聚类 (04_subcluster)")
    parser_subcluster.add_argument("-c", "--config", default="base_config.yaml", help="基础配置")
    parser_subcluster.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_subcluster.add_argument("--rds", default="", help="母集 RDS 对象；留空则从配置文件读取")
    parser_subcluster.add_argument("--subset_col", required=True, help="代表抽取亚群依据所在列名，例如 res0.4")
    parser_subcluster.add_argument("--subset_val", required=True, help="必须填入提取的值，如 '1' 或 '1,5'")
    parser_subcluster.add_argument("--methods", default="re_harmony,original_integrated", help="重聚类使用的数据整合方式")
    parser_subcluster.add_argument("--resolutions", default="0.2,0.4", help="为这个亚群重新跑的多个分别率探索值")
    parser_subcluster.add_argument("--col_sample", default="Sample", help="统计对照样本列")
    parser_subcluster.add_argument("--col_group", default="Group", help="统计大组列")
    parser_subcluster.add_argument("--integrate_by", default="Sample", help="整合去除批次时的依据列")
    parser_subcluster.add_argument("--nfeatures", type=int, default=1500, help="亚群变异点数(往往稍微减小)")
    parser_subcluster.add_argument("--npcs", type=int, default=20, help="重做PCA维数(往往稍微减小)")
    parser_subcluster.add_argument("--do_cluster", action="store_true", help="当场渲染这部分亚群的特异性热图及功能富集")
    parser_subcluster.add_argument("--do_refmarker", action="store_true", help="根据亚群指定特化Marker探查")
    parser_subcluster.add_argument("--refmarker_file", default="", help="针对该亚群特质准备的探针 Marker 列表位置")

    # ---------------------------------------------------------
    # 05_celltype
    # ---------------------------------------------------------
    parser_celltype = subparsers.add_parser("celltype", help="细胞学专家组注释映射分发与结题报告 (05_celltype)")
    parser_celltype.add_argument("-c", "--config", default="base_config.yaml", help="基础配置")
    parser_celltype.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_celltype.add_argument("--rds", default="", help="终态降维后未注释对象RDS；留空则从配置文件读取")
    parser_celltype.add_argument("--annofile", default="", help="关键！打签校准对照映射文件")
    parser_celltype.add_argument("--louper_path", default="", help="转化为 10X cloupe 浏览器数据的工具位置(可留空跳过)")
    
    # 因为开关逻辑在 Python 用 flags 非常顺畅，默认它们全开，可以用 flag 手动关闭
    parser_celltype.add_argument("--skip_cluster", action="store_false", dest="do_cluster", help="跳过自动渲染群组属性(节约时间/或被智能软链替代)")
    parser_celltype.add_argument("--skip_refmarker", action="store_false", dest="do_refmarker", help="跳过参考Marker点阵回图")
    parser_celltype.add_argument("--skip_annotation", action="store_false", dest="do_annotation", help="完全忽略读取映射字典，不做注释")
    parser_celltype.add_argument("--skip_celltype", action="store_false", dest="do_celltype", help="忽略细胞大类占比分化与UMAP产出")
    parser_celltype.add_argument("--do_celltype_de", action="store_true", help="对所有最终命名的细胞做一次大规模型号比对富集")
    parser_celltype.add_argument("--skip_report", action="store_false", dest="do_report", help="跳过最终结题炫酷动态长网页 HTML 渲染与多格式转化步骤")
    parser_celltype.add_argument("--force_clean", action="store_true", help="是否暴力格式化目标输出位点")

    # ---------------------------------------------------------
    # 07_differential
    # ---------------------------------------------------------
    parser_diff = subparsers.add_parser("differential", help="深度差异分析挖掘 (07_differential)")
    parser_diff.add_argument("-c", "--config", default="base_config.yaml", help="基础配置")
    parser_diff.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_diff.add_argument("--rds", default="", help="输入 Seurat RDS 对象；留空则从配置文件读取")
    parser_diff.add_argument("--method", default="", choices=["", "all_markers", "group_comparison"], help="差异模式: all_markers 或 group_comparison；留空则从配置文件读取")
    parser_diff.add_argument("--groupby", default="", help="比较列 (如 Group)；留空则从配置文件读取")
    parser_diff.add_argument("--split_by", default="", help="[模式B] 亚群拆分依据列；留空则从配置文件读取")
    parser_diff.add_argument("--cmp_file", default="", help="[模式B] 比较列表文件 (.xls/.txt)")
    parser_diff.add_argument("--species", default="", help="物种信息 (决定差异分析会自动挂载哪个富集库)；留空则从配置文件读取")
    parser_diff.add_argument("--do_enrich", action="store_true", default=None, help="是否执行 GO/KEGG 富集")

    # ---------------------------------------------------------
    # 06_merge_sub
    # ---------------------------------------------------------
    parser_merge = subparsers.add_parser("merge_sub", help="亚群注释回溯整合与大群清洗 (Step 06)")
    parser_merge.add_argument("-c", "--config", default="base_config.yaml", help="配置文件路径")
    parser_merge.add_argument("-m", "--main_rds", default="", help="主对象 RDS 路径；留空则从配置文件读取")
    parser_merge.add_argument("-s", "--sub_map", default="", help="亚群映射 (Major:Path, 逗号分隔)")
    parser_merge.add_argument("-o", "--outdir", default=".", help="输出目录")
    parser_merge.add_argument("--major_col", default="", help="大群列名；留空则从配置文件读取")
    parser_merge.add_argument("--sub_col", default="", help="亚群对象中的注释列名；留空则从配置文件读取")
    parser_merge.add_argument("--subtype_col", default="", help="回填到主对象的新列名；留空则从配置文件读取")
    parser_merge.add_argument("--sample_col", default="", help="样本列名；留空则从配置文件读取")
    parser_merge.add_argument("--group_col", default="", help="分组列名，多个用逗号分隔；留空则从配置文件读取")

    # ---------------------------------------------------------
    # init
    # ---------------------------------------------------------
    parser_init = subparsers.add_parser("init", help="初始化：拷贝配置文件模板到当前目录")
    parser_init.add_argument("-n", "--name", default="base_config.yaml", help="生成的配置文件名 (默认 base_config.yaml)")

    # ---------------------------------------------------------
    # run_all 向上兼容与统一 YAML
    # ---------------------------------------------------------
    parser_all = subparsers.add_parser("run_all", help="完全保留传统全栈配置文件驱动的串联流水执行")
    parser_all.add_argument("-c", "--config", required=True, help="含有全部参数的全量 config_template.yaml")
    parser_all.add_argument("-o", "--outdir", default=".", help="生成的 bash 运行脚本所存放的目录位置 (默认当前目录)")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)
        
    elif args.command == "run_all":
        # Call the classic generator
        scMax_pipeline.generate_bash_script(args.config, args.outdir)
        
    elif args.command == "qc":
        conf = {
            "rdsfile": args.rds,
            "samples_file": args.samples_file,
            "method": args.method,
            "metadatafile": args.meta,
            "col_sample": args.col_sample,
            "col_group": args.col_group,
            "species": args.species,
            "doub_rate": args.doub_rate,
            "mt_genes": args.mt_genes,
            "hb_genes": args.hb_genes,
            "run_doublet": args.run_doublet,
            "run_decontx": args.run_decontx
        }
        run_step("01_scQC", args, conf)

    elif args.command == "filter":
        conf = {
            "col_sample": args.col_sample,
            "min_feat": args.min_feat,
            "max_feat": args.max_feat,
            "min_count": args.min_count,
            "max_count": args.max_count,
            "max_mt": args.max_mt,
            "max_hb": args.max_hb,
            "rm_doublet": args.rm_doublet,
            "max_decontX": args.max_decontX
        }
        if args.rds != "":
            conf["rdsfile"] = args.rds
        run_step("02_filter", args, conf)
        
    elif args.command == "cluster":
        conf = {
            "subset_col": args.subset_col,
            "subset_val": args.subset_val,
            "methods": args.methods,
            "resolutions": args.resolutions,
            "col_sample": args.col_sample,
            "col_group": args.col_group,
            "integrate_by": args.integrate_by,
            "nfeatures": args.nfeatures,
            "npcs": args.npcs,
            "do_cluster": args.do_cluster,
            "do_refmarker": args.do_refmarker,
            "refmarker_file": args.refmarker_file
        }
        if args.rds != "":
            conf["rdsfile"] = args.rds
        run_step("03_cluster", args, conf)

    elif args.command == "subcluster":
        conf = {
            "subset_col": args.subset_col,
            "subset_val": args.subset_val,
            "methods": args.methods,
            "resolutions": args.resolutions,
            "col_sample": args.col_sample,
            "col_group": args.col_group,
            "integrate_by": args.integrate_by,
            "nfeatures": args.nfeatures,
            "npcs": args.npcs,
            "do_cluster": args.do_cluster,
            "do_refmarker": args.do_refmarker,
            "refmarker_file": args.refmarker_file
        }
        if args.rds != "":
            conf["rdsfile"] = args.rds
        run_step("04_subcluster", args, conf)

    elif args.command == "celltype":
        conf = {
            "annofile": args.annofile,
            "louper_path": args.louper_path,
            "do_cluster": args.do_cluster,
            "do_refmarker": args.do_refmarker,
            "do_annotation": args.do_annotation,
            "do_celltype": args.do_celltype,
            "do_celltype_de": args.do_celltype_de,
            "do_report": args.do_report,
            "force_clean": args.force_clean
        }
        if args.rds != "":
            conf["rdsfile"] = args.rds
        run_step("05_celltype", args, conf)

    elif args.command == "differential":
        conf = {}
        if args.rds != "":
            conf["rdsfile"] = args.rds
        if args.method != "":
            conf["method"] = args.method
        if args.groupby != "":
            conf["groupby"] = args.groupby
        if args.split_by != "":
            conf["split_by"] = args.split_by
        if args.cmp_file != "":
            conf["cmp_file"] = args.cmp_file
        if args.do_enrich is not None:
            conf["do_enrich"] = args.do_enrich
        run_step("07_differential", args, conf)

    elif args.command == "merge_sub":
        conf = {}
        if args.main_rds != "":
            conf["main_rds"] = args.main_rds
        if args.sub_map != "":
            conf["sub_map"] = args.sub_map
        if args.major_col != "":
            conf["major_col"] = args.major_col
        if args.sub_col != "":
            conf["sub_col"] = args.sub_col
        if args.subtype_col != "":
            conf["subtype_col"] = args.subtype_col
        if args.sample_col != "":
            conf["sample_col"] = args.sample_col
        if args.group_col != "":
            conf["group_col"] = args.group_col
        run_step("06_merge_sub", args, conf)

    elif args.command == "init":
        base_dir = os.path.dirname(os.path.abspath(__file__))
        import shutil
        
        # 定义所有可用的配置模板映射 (源文件名: 目标建议名)
        templates = {
            "config_template.yaml":             "base_config.yaml",
            "config_template_01_qc.yaml":       "01_scQC.yaml",
            "config_template_02_filter.yaml":   "02_filter.yaml",
            "config_template_03_cluster.yaml":  "03_cluster.yaml",
            "config_template_04_subcluster.yaml":"04_subcluster.yaml",
            "config_template_05_celltype.yaml": "05_celltype.yaml",
            "config_template_07_diff.yaml":     "07_differential.yaml",
            "config_template_06_merge_sub.yaml":"06_merge_sub.yaml",
            "config_template_06_diff.yaml":     "07_differential.yaml"
        }
        
        generated_files = []
        for src, dest in templates.items():
            src_path = os.path.join(base_dir, src)
            if os.path.exists(src_path):
                # 如果用户通过 -n 指定了名字，则将主模板命名为该名字
                target_name = args.name if (src == "config_template.yaml" and args.name != "base_config.yaml") else dest
                
                if not os.path.exists(target_name):
                    try:
                        shutil.copy(src_path, target_name)
                        generated_files.append(target_name)
                    except Exception as e:
                        print(f"Warning: 拷贝 {target_name} 失败: {e}")
                else:
                    # 如果已存在，则不覆盖，但后续需要提示
                    pass

        if generated_files:
            print("✨ 成功为您初始化 scMax 步进式精简配置全家桶：")
            for f in sorted(generated_files):
                # 简单映射一下中文描述
                desc_map = {
                    "base_config.yaml": "全流程标准版 (一站式控制)",
                    "01_scQC.yaml": "Step 01: 数据整合与预质控",
                    "02_filter.yaml": "Step 02: 深度过滤与质控",
                    "03_cluster.yaml": "Step 03: 降维、去批次与聚类",
                    "04_subcluster.yaml": "Step 04: 特定亚群提取与细分",
                    "05_celltype.yaml": "Step 05: 专家注释、分发与 HTML 报告",
                    "06_merge_sub.yaml": "Step 06: 亚群注释回填整合",
                    "07_differential.yaml": "Step 07: 差异分析挖掘与组间比较"
                }
                desc = desc_map.get(f, "")
                print(f"  - {f:<20} {desc}")
            print("\n💡 建议：您可以根据当前分析进度，仅编辑并执行对应的配置文件，消除冗余参数干扰。")
            print("   例如：python3 scMax.py cluster -c 03_cluster.yaml")
        else:
            print("Notice: 未生成任何新文件（目标分步配置文件均已存在）。")

if __name__ == "__main__":
    main()
