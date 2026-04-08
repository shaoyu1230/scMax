#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
scMax 第 5 步 - 动态报告 Jupyter Notebook 组装器
作者: Antigravity
说明: 本脚本根据 5 个处理开关生成相对路径的高保真自适应 HTML 前置 ipynb。
"""

import json
import argparse
import os

def create_markdown_cell(source_lines):
    """辅助函数：创建标准的 Jupyter Markdown Cell JSON 结构"""
    # 确保除了最后一行，前面的行都带有 \n
    formatted_source = []
    for i, line in enumerate(source_lines):
        if i < len(source_lines) - 1 and not line.endswith('\n'):
            formatted_source.append(line + '\n')
        else:
            formatted_source.append(line)
            
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": formatted_source
    }

def main():
    parser = argparse.ArgumentParser(description="动态组装 scMax 细胞注释结题报告 (ipynb格式)")
    parser.add_argument("--outdir", required=True, help="分析结果储存目录(影响相对路径基点)")
    parser.add_argument("--do_cluster", action="store_true", help="是否含聚类分析")
    parser.add_argument("--do_refmarker", action="store_true", help="是否含参考Marker表达分析")
    parser.add_argument("--do_annotation", action="store_true", help="是否含细胞类型初阶注释")
    parser.add_argument("--do_celltype", action="store_true", help="是否含高级细胞类型比例统计")
    parser.add_argument("--do_celltype_de", action="store_true", help="是否含细胞类型差异基因富集")
    parser.add_argument("--no_group", action="store_true", help="是否无分组信息")
    
    args = parser.parse_args()
    
    # === 构建基础骨架 ===
    notebook = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "R",
                "language": "R",
                "name": "ir"
            },
            "language_info": {
                "codemirror_mode": "r",
                "file_extension": ".r",
                "mimetype": "text/x-r-source",
                "name": "R",
                "pygments_lexer": "r",
                "version": "4.3.1"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 5
    }
    
    from typing import List, Dict, Any, cast
    cells = cast(List[Dict[str, Any]], notebook["cells"])
    
    # --- 引言部分 ---
    cells.append({
        "cell_type": "raw",
        "metadata": {},
        "source": ["title: 安诺优达单细胞转录组细胞类型注释报告"]
    })
    
    cells.append(create_markdown_cell([
        "高通量单细胞转录组的出现，让转录组研究正式从组织维度下探到单个细胞（群）精度，在各类研究领域中得到了广泛的应用。而细胞注释是理解单细胞数据的第一步，它通过区分细胞类型和状态来奠定后续分析的基础。  ",
        "",
        "因此，围绕细胞注释结果解读，本报告主要分为五部分：‌  ",
        "1、降维聚类可视化‌：直观体现各类群细胞特征。  ",
        "‌2、特异性标志物定位‌：验证客户提供Marker基因的亚群表达特征。 （如未提供参考 marker 基因则无此部分结果）  ",
        "3、人工&数据库校准‌：人工校正和安诺数据库辅助验证，得出准确注释结果。  ",
        "‌4、组间差异比较：组间细胞类型比例比较和偏好性分析，找出关键差异亚群。  ",
        "‌5、细胞功能分析‌：富集分析不同细胞类型的功能特征，找出关键分子通路。  ",
        "其中第三、四、五部分生成的图表可直接用于论文发表，显著提升分析的深度与实用性。  "
    ]))
    
    # --- 模块 1: 聚类分析 (do_cluster) ---
    has_cluster = args.do_cluster or os.path.exists(os.path.join(args.outdir, "cluster_characterization"))
    if has_cluster:
        cells.append(create_markdown_cell([
            "# 一、聚类分析",
            "",
            "通过对单细胞转录组数据进行降维聚类，量化亚群规模（细胞数量/占比）及分布特征，同时采用‌Top10差异基因热图‌可视化亚群异质性，初步锚定类群特异性基因表达谱，为后续细胞类型注释奠定基础。"
        ]))
        cells.append(create_markdown_cell([
            "## 1.1 细胞聚类分析",
            "",
            "（1）对数据进行降维聚类得到多个细胞Cluster，进行UMAP分析，将聚类结果在UMAP映射中展示：",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"cluster_characterization/1_seurat_clusters.UMAP.png\" alt=\"图片1\" width=\"400\"></td>",
            "    <td><img src=\"cluster_characterization/1_seurat_clusters.UMAP-split.png\" alt=\"图片2\" width=\"400\"></td>",
            "  </tr>",
            "</table>",
            "",
            "<p style=\"font-size: 12px;\">注：左图为降维聚类后的UMAP结果图，右图为各细胞群的UMAP分布图，不同颜色代表不同的Cluster。</p>"
        ]))
        if args.no_group:
            cells.append(create_markdown_cell([
                "（2）统计Cluster中每个细胞中的nCount_RNA（number of UMI），nFeature_RNA（number of Gene）和percent.mt（线粒体基因占nFeature_RNA的比例），将统计结果用小提琴图展示：",
                "",
                "<img src=\"cluster_characterization/0_QC.png\" alt=\"质控\" width=\"500\">",
                "<p style=\"font-size: 12px;\">注：小提琴图的横坐标代表样本，从上至下分别为unique UMI总数、基因表达数目和线粒体基因比例的结果。</p>",
                "",
                "（3）将样本分布映射到UMAP上，展示整体分布以及各样本分布情况：",
                "",
                "<img src=\"cluster_characterization/1_Sample.UMAP.png\" alt=\"样本\" width=\"600\">"
            ]))
        else:
            cells.append(create_markdown_cell([
                "（2）统计Cluster中每个细胞中的nCount_RNA（number of UMI），nFeature_RNA（number of Gene）和percent.mt（线粒体基因占nFeature_RNA的比例），将统计结果用小提琴图展示：",
                "",
                "<img src=\"cluster_characterization/0_QC.png\" alt=\"质控\" width=\"500\">",
                "<p style=\"font-size: 12px;\">注：小提琴图的横坐标代表样本，从上至下分别为unique UMI总数、基因表达数目和线粒体基因比例的结果。</p>",
                "",
                "（3）将样本分布映射到UMAP上，展示整体分布以及各样本分布情况：",
                "",
                "<table>",
                "  <tr>",
                "    <td><img src=\"cluster_characterization/1_Sample.UMAP.png\" alt=\"样本\" width=\"400\"></td>",
                "    <td><img src=\"cluster_characterization/2_Sample.UMAP-split.png\" alt=\"散开\" width=\"400\"></td>",
                "  </tr>",
                "</table>",
                "<p style=\"font-size: 12px;\">注：左图为分样本UMAP结果图，不同颜色代表不同的样本，右图分别展示不同样本的UMAP结果图，不同颜色代表不同的Cluster。</p>",
                "",
                "（4）将分组细胞映射到UMAP上，展示整体分布及各组分布情况：",
                "",
                "<table>",
                "  <tr>",
                "    <td><img src=\"cluster_characterization/1_Group.UMAP.png\" alt=\"分组\" width=\"400\"></td>",
                "    <td><img src=\"cluster_characterization/2_Group.UMAP-split.png\" alt=\"散开\" width=\"400\"></td>",
                "  </tr>",
                "</table>",
                "<p style=\"font-size: 12px;\">注：左图为分组UMAP结果图，不同颜色代表不同的分组，右图分别展示不同分组的UMAP结果图,不同颜色代表不同的Cluster。</p>"
            ]))
        
        cells.append(create_markdown_cell([
            "## 1.2 细胞群数量和比例统计",
            "",
            "（1）统计样本中各Cluster细胞数量和比例，绘制数量和比例分布热图：",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"cluster_characterization/3_Count-seurat_clusters.in.Sample-heatmap.png\" alt=\"图片1\" width=\"400\"></td>",
            "    <td><img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Sample-heatmap.png\" alt=\"图片2\" width=\"400\"></td>",
            "  </tr>",
            "</table>",
            "<p style=\"font-size: 12px;\">注：横轴为细胞Cluster号，纵轴为样本名；左图为样本中各细胞Cluster的细胞数目，右图为样本中细胞Cluster所占的百分比。</p>"
        ]))
        
        if args.no_group:
            cells.append(create_markdown_cell([
                "（2）统计样本中各Cluster细胞占比绘制比例柱状图和趋势图：",
                "",
                "<img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Sample-bar.png\" alt=\"图1 构成\" width=\"600\">",
                "<p style=\"font-size: 12px;\">注：每个样本中各Cluster的细胞百分比堆积柱形图和趋势图，横轴为样本名，纵轴为同一样本中不同Cluster的细胞百分比。</p>"
            ]))
        else:
            cells.append(create_markdown_cell([
                "（2）统计样本中各Cluster细胞占比绘制比例柱状图和趋势图：",
                "",
                "<table>",
                "  <tr>",
                "    <td><img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Sample-bar.png\" alt=\"图片1\" width=\"400\"></td>",
                "    <td><img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Sample-trend.png\" alt=\"图片2\" width=\"400\"></td>",
                "  </tr>",
                "</table>",
                "<p style=\"font-size: 12px;\">注：每个样本中各Cluster的细胞百分比堆积柱形图和趋势图，横轴为样本名，纵轴为同一样本中不同Cluster的细胞百分比。</p>",
                "",
                "（3）统计分组中各Cluster细胞占比绘制比例柱状图，如下：",
                "",
                "<table>",
                "  <tr>",
                "    <td><img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Group-bar.png\" alt=\"图片1\" width=\"400\"></td>",
                "    <td><img src=\"cluster_characterization/3_Frac-seurat_clusters.in.Group-trend.png\" alt=\"图片2\" width=\"400\"></td>",
                "  </tr>",
                "</table>",
                "<p style=\"font-size: 12px;\">注：每个分组中不同Cluster的细胞百分比堆积柱形图和趋势图，横轴为不同的分组，纵轴为同一分组中不同Cluster的细胞百分比。不同颜色代表不同的Cluster。</p>"
            ]))
            
        cells.append(create_markdown_cell([
            "## 1.3 细胞群差异基因分析",
            "",
            "聚类分析得到了不同的细胞群，为了进一步研究每个细胞群的特征，使用Seurat包的FindAllMarkers函数进行了细胞群差异基因分析。差异基因分析过程如下（以Cluster1为例）：计算Cluster1中某个基因的平均表达量；计算所有Cluster中除了Cluster1的对应基因的平均表达量；采用wilcox统计方法进行显著性检验。显著差异基因筛选的标准为：min.pct>=0.10，avg_log2FC>=0.25，p_val_adj<0.05。",
            "",
            "（1）统计每个细胞群中的Marker基因；（2）绘制各细胞群Top10差异基因气泡图：",
            "",
            "<img src=\"cluster_characterization/4_Cluster.DE/2_ClusterTop10DEGs.png\" alt=\"Markdown Logo\" width=\"1000\">",
            "<p style=\"font-size: 12px;\">注：选择所有聚类（Cluster）中avg_log2FC数值top10的marker基因（去重后），绘制在不同Cluster中的表达量dotplot图。  ",
            "横轴为不同基因，纵轴为聚类（Cluster）。点的颜色表示表达量高低，从蓝到红表示表达量从低到高，即越红表示表达量越高。点的大小表示某Cluster中有该基因表达的细胞占比，点越大，说明细胞占比越高。</p>",
            "",
            "## 1.4 细胞群差异基因富集分析",
            "取各细胞群中显著差异（p_val_adj<0.05，avg_log2FC>=0.25，pct.1>=0.1）的基因，使用clusterProfiler包进行富集分析。报告中以Cluster1为例（可在 5_Enrich 目录查看全集）：",
            "",
            "（1）GO富集分析结果：",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"cluster_characterization/5_Enrich/1.up.GO_bar.png\" alt=\"图片1\" width=\"400\"></td>",
            "    <td><img src=\"cluster_characterization/5_Enrich/1.up.GO_dot.png\" alt=\"图片2\" width=\"400\"></td>",
            "  </tr>",
            "</table>"
        ]))

    # --- 模块 2: 参考Marker验证 (do_refmarker) ---
    has_refmarker = args.do_refmarker or os.path.exists(os.path.join(args.outdir, "marker_expression"))
    if has_refmarker:
        cells.append(create_markdown_cell([
            "# 二、参考marker表达情况",
            "",
            "整合客户提供的Marker基因与聚类数据，通过‌跨亚群表达热图展示其在特定Cluster的特异性表达（如log2FC>2，p<0.001）情况。和临近亚群相比，若Marker在目标Cluster呈现显著表达断层，则它可能优先作为该细胞类型的注释依据，同步纳入终版注释结果。",
            "",
            "（1）绘制参考Marker在各细胞群的表达水平气泡图，可直观得看到各细胞类型Marker在细胞群的表达情况："
        ]))
        cells.append(create_markdown_cell([
            "<img src=\"marker_expression/1_DotPlot-RefMarkers.png\" alt=\"RefMarker\" width=\"800\" >",
            "<p style=\"font-size: 12px;\">注：客户提供的Marker基因在不同Cluster中表达量的dotplot热图，横轴为聚类（Cluster），纵轴右侧标签为基因，纵轴左侧标签为Marker基因所属的细胞类型。点/方块的颜色表示表达量高低，从蓝到红表示表达量从低到高，即越红表示表达量越高。点的大小表示某Cluster中有该基因表达的细胞占比，点越大，说明有表达的细胞占比越高。</p>",
            "",
            "（2）绘制Marker基因表达UMAP图，可直观得看到基因整体表达情况：",
            "",
            "<img src=\"marker_expression/2_FeaturePlot.example.png\" alt=\"FeaturePlot\" width=\"500\">",
            "<p style=\"font-size: 12px;\">注：客户提供的Marker基因在不同Cluster中表达量的feature plot图。红色代表高表达，蓝色代表低表达。</p>"
        ]))
        
    # --- 模块 3: 细胞类型注释 (do_annotation) ---
    has_annotation = args.do_annotation or os.path.exists(os.path.join(args.outdir, "annotation"))
    if has_annotation:
        cells.append(create_markdown_cell([
            "# 三、细胞类型注释",
            "",
            "基于降维聚类各亚群特征基因，跨库比对公共数据库（CellMarker/PanglaoDB）与安诺内部人工注释库‌，同步整合客户Marker基因表达验证，经人工校验各细胞群特征，输出细胞分群清晰和注释准确的终版结果。",
            "",
            "## 3.1 所有细胞类型注释结果",
            "细胞Cluster注释结果及参考marker信息如下表（仅示例图）：",
            "",
            "<img src=\"annotation/CellType_All/0_Cluster-CellType.anno.png\" alt=\"Anno Table\" width=\"600\">"
        ]))
        cells.append(create_markdown_cell([
            "（1）将细胞Cluster加上细胞类型注释信息，获得所有细胞类型注释结果：",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"annotation/CellType_All/1_CellType.UMAP-blank.png\" alt=\"全覆盖注释\" width=\"400\"></td>",
            "    <td><img src=\"annotation/CellType_All/1_CellType.UMAP-split.png\" alt=\"切分\" width=\"400\"></td>",
            "  </tr>",
            "</table>",
            "（2）将样本/分组信息映射到UMAP上，展示其分布情况：",
            "",
            "<img src=\"annotation/CellType_All/2_Sample.UMAP-split.png\" alt=\"Sample\" width=\"1000\">",
            "<img src=\"annotation/CellType_All/2_Group.UMAP-split.png\" alt=\"Group\" width=\"1000\">",
            "<p style=\"font-size: 12px;\">注：注释后的细胞UMAP图谱，不同颜色代表不同的细胞/样本/分组。</p>"
        ]))
        cells.append(create_markdown_cell([
            "## 3.2 所有细胞类型注释 marker 基因表达展示",
            "分析细胞类型定义所使用的主要Marker基因在所有细胞群中的表达情况，用于验证注释的准确性，以及Marker的特异性。",
            "",
            "绘制气泡图，展示Marker基因的表达水平及在细胞类型中表达比例：",
            "",
            "<img src=\"annotation/CellType_All/3_CellType.dotplot.clusterAnno.png\" alt=\"AllAnnoDot\" width=\"600\">",
            "<p style=\"font-size: 12px;\">注：细胞类型定义使用的Marker基因在不同细胞类型中表达量的dotplot图，纵轴为定义的细胞类型及聚类号，横轴为基因名。</p>"
        ]))
        cells.append(create_markdown_cell([
            "## 3.3 高质量细胞类型注释结果",
            "去除部分不展示的细胞，比如双胞（Doublets）、低质量细胞（Low quality cells）等，得到高质量细胞类型数据。",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"annotation/CellType_Keep/1_CellType.UMAP-blank.png\" alt=\"全覆盖注释\" width=\"400\"></td>",
            "    <td><img src=\"annotation/CellType_Keep/1_CellType.UMAP-split.png\" alt=\"切分\" width=\"400\"></td>",
            "  </tr>",
            "</table>",
            "<p style=\"font-size: 12px;\">注：去除双细胞和低质量细胞的umap图谱，不同颜色代表不同的细胞类型。</p>"
        ]))
        cells.append(create_markdown_cell([
            "## 3.4 细胞类型定义使用的marker基因表达分布",
            "分析细胞类型定义所使用的主要Marker基因在注释好的细胞群中的表达情况，用于验证注释的准确性，以及Marker的特异性。",
            "",
            "（1）绘制气泡图，展示Marker基因的表达水平及在细胞类型中表达比例：",
            "",
            "<img src=\"annotation/CellType_Keep/4_CellType.dotplot.png\" alt=\"AnnoDot\" width=\"600\">",
            "<p style=\"font-size: 12px;\">注：细胞类型定义使用的Marker基因在不同细胞类型中表达量的dotplot图，纵轴为定义的细胞类型，横轴为基因名。</p>",
            "",
            "（2）绘制热图，展示细胞类型中Marker基因的表达水平：",
            "",
            "<img src=\"annotation/CellType_Keep/3_CellType.heatmap.png\" alt=\"AnnoHeat\" width=\"600\">"
        ]))
        
    # --- 模块 4: 分类比例推算 (do_celltype) ---
    has_celltype = args.do_celltype or os.path.exists(os.path.join(args.outdir, "celltype_fraction"))
    if has_celltype:
        cells.append(create_markdown_cell([
            "# 四、细胞类型比例分析",
            "",
            "基于完整细胞注释，量化跨组细胞丰度差异，明确组间差异模式并锁定显著变化的细胞亚群，为后续机制解析或标志物挖掘提供定向数据支撑。",
            "## 4.1 细胞类型数量统计",
            "为进一步研究样本在细胞类型分布上的差异，统计了细胞类型在各样本中的数量和占比情况，绘制热图："
        ]))
        cells.append(create_markdown_cell([
            "<table>",
            "  <tr>",
            "    <td><img src=\"celltype_fraction/2_CellType_Keep.in.Sample-count-heatmap.png\" alt=\"Count\" width=\"400\"></td>",
            "    <td><img src=\"celltype_fraction/2_CellType_Keep.in.Sample-frac-heatmap.png\" alt=\"Frac\" width=\"400\"></td>",
            "  </tr>",
            "</table>"
        ]))

        if args.no_group:
            cells.append(create_markdown_cell([
                "统计样本中各细胞类型的占比绘制比例分类图：",
                "",
                "<img src=\"celltype_fraction/2_CellType_Keep.in.Sample-bar.png\" alt=\"样本组成\" width=\"600\">"
            ]))
        else:
            cells.append(create_markdown_cell([
                "统计样本/分组中各细胞类型的占比绘制比例分类图及趋势图：",
                "",
                "<table>",
                "  <tr>",
                "    <td><img src=\"celltype_fraction/2_CellType_Keep.in.Sample-bar.png\" alt=\"样本组成\" width=\"400\"></td>",
                "    <td><img src=\"celltype_fraction/2_CellType_Keep.in.Group-bar.png\" alt=\"分组宏观\" width=\"400\"></td>",
                "  </tr>",
                "</table>",
                "## 4.2 细胞类型偏好性分析（Ro/e）",
                "Ro/e（the ratio of observed over expected cell numbers）是一种用于评估细胞亚群在组织或分组中分布偏好性的统计学方法（大于1为表现出富集关联）：",
                "",
                "<img src=\"celltype_fraction/3_CellType_Keep.in.Group-Roe.png\" alt=\"Roe\" width=\"400\">",
                "<p style=\"font-size: 12px;\">注：横轴为细胞群名称，纵轴为组织类型，颜色代表Ro/e值。</p>"
            ]))
        
    # --- 模块 5: 分类间功能富集 (do_celltype_de) ---
    has_celltype_de = args.do_celltype_de or os.path.exists(os.path.join(args.outdir, "celltype_characterization"))
    if has_celltype_de:
        cells.append(create_markdown_cell([
            "# 五、细胞类型特征基因表达分析",
            "",
            "基于已定义的细胞亚群及其特异性标记基因，通过GO/KEGG富集分析检测亚群内显著富集的信号通路，将基因表达特征关联至分子功能与表型变化（如免疫应答或代谢重塑），验证其与疾病进程或药物干预的潜在关联，为后续机制研究提供分子层面的靶点与理论依据。",
            "",
            "## 5.1 细胞类型特征基因计算",
            "经过注释得到了多种细胞类型，为了进一步研究各细胞类型的特征，使用Seurat包的FindAllMarkers函数进行了细胞类型差异基因分析。显著差异基因筛选的标准为：min.pct>=0.10，avg_log2FC>=0.25，p_val_adj<0.05。"
        ]))
        cells.append(create_markdown_cell([
            "展示各细胞群Top5差异基因与热图：",
            "",
            "<img src=\"celltype_characterization/1_CellType_DE/3_CellType.Top5DEGs.png\" alt=\"大类Top\" width=\"800\">",
            "<img src=\"celltype_characterization/1_CellType_DE/2_CellType.Top5DEGs-heatmap.png\" alt=\"TopHeatmap\" width=\"500\">",
            "## 5.2 细胞类型特征基因富集分析",
            "取各细胞类型中avg_log2FC数值Top100的基因，使用clusterProfiler包进行GO与KEGG富集分析：",
            "",
            "<table>",
            "  <tr>",
            "    <td><img src=\"celltype_characterization/2_CellType_Enrich/2_CellType.GO.png\" alt=\"GO\" width=\"450\"></td>",
            "    <td><img src=\"celltype_characterization/2_CellType_Enrich/4_CellType.KEGG.png\" alt=\"KEGG\" width=\"450\"></td>",
            "  </tr>",
            "</table>"
        ]))
        
    # --- 结语 ---
    cells.append(create_markdown_cell([
        "# 结果文件目录结构与软件声明",
        "",
        "```text",
        "├── cluster_characterization         cluster特征：包括cluster分布，样本/分组分布，cluster占比等",
        "├── marker_expression                参考marker表达绘图",
        "├── annotation                       细胞类型注释结果",
        "│   ├── CellType_All                   所有细胞注释结果",
        "│   └── CellType_Keep                  去除低质量、doublet等不感兴趣的细胞后，细胞注释结果",
        "├── celltype_fraction                细胞占比分析",
        "├── celltype_characterization        细胞类型特征，包括差异分析和富集分析",
        "└── CellType.Annotation_report.html    细胞类型注释报告",
        "```",
        "",
        "本次分析所使用的核心软件：",
        "",
        "| 软件名称            | 版本号    | 功能描述                              |",
        "|---------------------|-----------|---------------------------------------|",
        "| R                   | 4.3.1     | R语言软件                             |",
        "| R Seurat            | 4.0.0       | 单细胞数据分析工具                    |",
        "| R clusterProfiler   | 4.10.0       | GO和KEGG富集分析工具                  |"
    ]))
    
    cells.append(create_markdown_cell([
        "# 参考文献",
        "",
        "Abdi H, Williams L J. Principal component analysis[J]. Wiley interdisciplinary reviews: computational statistics, 2010, 2(4): 433-459.  ",
        "",
        "Butler A, Hoffman P, Smibert P, et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species[J]. Nature biotechnology, 2018, 36(5): 411.  ",
        "",
        "Camp J G, Sekine K, Gerber T, et al. Multilineage communication regulates human liver bud development from pluripotency[J]. Nature, 2017, 546(7659): 533.  ",
        "",
        "Ding C, He X. K-means clustering via principal component analysis[C]//Proceedings of the twenty-first international conference on Machine learning. ACM, 2004: 29.  ",
        "",
        "Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21.  ",
        "",
        "Maaten L, Hinton G. Visualizing data using t-SNE[J]. Journal of Machine Learning Research, 2008, 9(Nov): 2579-2605.  ",
        "",
        "Macosko E Z, Basu A, Satija R, et al. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets.[J]. Cell, 2015, 161(5):1202.  ",
        "",
        "Hinton G, Roweis S. Stochastic Neighbor Embedding[J]. Advances in Neural Information Processing Systems, 2010, 41(4):833--840.  ",
        "",
        "Shannon, P. et al. Cytoscape: A software environment for integrated models of biomolecular interaction networks. Genome Research 13, 2498-2504 (2003).  ",
        "",
        "Van der Maaten, L.J.P. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15, 3221-3245 (2014).  ",
        "",
        "Zhong S, Zhang S, Fan X, et al. A single-cell RNA-seq survey of the developmental landscape of the human prefrontal cortex[J]. Nature, 2018, 555(7697): 524.  ",
        "",
        "Wang Y , Wang R , Zhang S , et al. iTALK: an R Package to Characterize and Illustrate Intercellular Communication. 2019.  ",
        "",
        "Cabello-Aguilar S , Alame M , Kon-Sun-Tack F , et al. SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics[J]. Nucleic Acids Research, 2020.  ",
        "",
        "Patel AP, Tirosh I, Trombetta JJ, et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science.2014;344(6190):1396-1401. doi:10.1126/science.1254257.  ",
        "",
        "La Manno G, Soldatov R, Zeisel A, et al.RNA velocity of single cells.Nature. 2018;560(7719):494-498. doi:10.1038/s41586-018-0414-6.  ",
        "",
        "Gao R , Bai S , Ying C H , et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes[J]. Nature Biotechnology, 2021:1-10.  ",
        "",
        "Zhang, L., Yu, X., Zheng, L. et al. Lineage tracking reveals dynamic relationships of T cells in colorectal cancer. Nature.2018; 564: 268–272.  "
    ]))

    # 写入文件
    target_ipynb = os.path.join(args.outdir, "report_custom.ipynb")
    os.makedirs(args.outdir, exist_ok=True)
    with open(target_ipynb, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, ensure_ascii=False, indent=1)
        
    print(f"[*] 依据 5 组开关状态成功自适应拼装 Notebook：{target_ipynb}")

if __name__ == "__main__":
    main()
