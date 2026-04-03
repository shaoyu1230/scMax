#!/usr/bin/env python3
import argparse
import sys
import os
import yaml
import scMax_pipeline

def load_base_config(config_path):
    if not os.path.exists(config_path):
        print(f"Error: Base configuration file '{config_path}' not found.")
        sys.exit(1)
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def run_step(step_name, args, step_conf):
    base_config = load_base_config(args.config)
    if args.outdir != ".":
        base_config['global_outdir'] = args.outdir

    # Disable all steps except the current one
    for s in ["01_scQC", "scQC", "02_filter", "03_cluster", "04_subcluster", "05_celltype", "Database_Info"]:
        if s in base_config and isinstance(base_config[s], dict):
            if s != step_name and s != "Database_Info":
                base_config[s]['run'] = False
                
    # Update or create the specific step config
    step_conf['run'] = True
    if step_name not in base_config:
        base_config[step_name] = {}
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
    parser_filter.add_argument("--rds", required=True, help="来自上一级的 RDS 结果路径")
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
    parser_cluster.add_argument("--rds", required=True, help="需要提取聚类的RDS对象路径")
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
    parser_subcluster.add_argument("--rds", required=True, help="母集 RDS 对象")
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
    parser_celltype.add_argument("--rds", required=True, help="终态降维后未注释对象RDS")
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
            "rdsfile": args.rds,
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
        run_step("02_filter", args, conf)
        
    elif args.command == "cluster":
        conf = {
            "rdsfile": args.rds,
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
        run_step("03_cluster", args, conf)

    elif args.command == "subcluster":
        conf = {
            "rdsfile": args.rds,
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
        run_step("04_subcluster", args, conf)

    elif args.command == "celltype":
        conf = {
            "rdsfile": args.rds,
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
        run_step("05_celltype", args, conf)

if __name__ == "__main__":
    main()
