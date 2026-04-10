# scMax v2.0: 高级单细胞转录组全栈分析与管理平台 🧬

**scMax** 是一套面向生产环境设计的端到端单细胞分析系统。它将生信底层计算（Seurat/SingleCellExperiment）、自动化报告拼装（Jupyter/Python）、多维数据资产入库（SQLite）以及 Web 可视化看板（Flask/Bootstrap）进行了深度垂直整合，旨在解决单细胞分析中“交付慢、数据碎、回溯难”的核心痛点。

---

## 🌟 核心模块架构 (Standard Modules)

项目采用严格的模块化层级体系，确保代码的高度解耦与可维护性：

- **`01_scQC`**: 数据整合与前置质控。支持多样本合并、双细胞判定（DoubletFinder）及环境污染去除（decontX）。
- **`02_filter`**: 深度过滤裁切。基于物种敏感阈值去除死细胞、低质量细胞及物理双细胞。
- **`03_cluster`**: 降维聚类核心。集成 Harmony、RPCA 等主流去批次算法并完成全局 Marker 挖掘。
- **`04_subcluster`**: 亚群精细分化。针对特定簇进行二次提取与特异性重聚类。
- **`05_celltype`**: 注释、分发与报告。实现多级人工标注映射，生成交互式 HTML 结题报告并自动注册中央数据库。
- **`06_merge_sub` [NEW]**: 亚群注释回填整合。将亚群注释结果合并回大群对象，输出层级 UMAP 与样本/分组统计。
- **`07_differential` [NEW]**: 深度差异挖掘。支持全局标签差异分析（Mode A）以及精细亚群内部的组间两两比较（Mode B），可直接基于合并后的对象继续分析。

---

## 🛠️ 快速开始

### 核心工具入口
- **`scMax.py`**: **[推荐]** 现代化的 CLI 工具入口，支持子命令驱动，适合日常测试与精细调节。
- **`scMax_pipeline.py`**: 底层流水线引擎，负责解析 YAML 配置并生成高性能 Bash 执行脚本。

### 场景一：全链路自动化生产 (Run All)
如果您已经配置好了 `base_config.yaml` 或项目专属 YAML，通过以下指令即可完成从 raw 矩阵到 HTML 报告及数据库入库的全流程：
```bash
python3 scMax.py run_all -c project_config.yaml -o ./analysis_scripts
```
> **Tip**: 生成的脚本会自动附带 `_YYYYMMDD_HHMMSS` 时间戳，防止多次运行相互覆盖。

### 场景二：分步精细化调试 (Step-by-Step)
您可以随时利用 CLI 系统对特定环节进行针对性重跑或参数微调：
```bash
# 例子：只重跑聚类阶段，调整分辨率和去批次方法
python3 scMax.py cluster --rds output_02/filtered_data.rds --methods harmony,rpca --resolutions 0.2,0.4,0.6 -o ./res_test

# 例子：进入专家注释环节
python3 scMax.py celltype --rds cluster_out.rds --annofile final_mapping.xlsx -o ./delivery
```

### 场景三：深度生物学挖掘 (Differential)
利用新增的 06 模块进行精细化的 Case vs Control 分析：
```bash
# 模式 B：在每个 CellType 内部比较组 1 vs 组 2
python3 scMax.py differential --rds celltype.rds --method group_comparison --groupby Group --split_by CellType --cmp_file pairs.txt
```

---

## 📊 数据资产管理与看板

scMax 不仅仅产出静态图表，更关注分析经验的“原子化”沉淀。

### 1. 中央数据库入库 (`scMax_db_manager.py`)
每次分析流程结束，系统会自动将分析元数据（项目编号、分析人、物种、核心图谱、统计表）注册至 `WebUI/scMax_projects.db`。同时，所有核心图表会同步软链/拷贝至资产隔离区 `WebUI/assets/`，确保数据的永久可回溯性。

### 2. 网页端可视化大盘 (Web看板)
启动本地服务，在 UI 界面中全局搜索、筛选并对比历次项目的分析结论：
```bash
python3 WebUI/app.py
```
默认访问地址：`http://127.0.0.1:5050`

---

## 📂 项目目录规约

```text
scMax/
├── scMax.py               # 🚀 CLI 统一入口
├── scMax_pipeline.py      # ⚙️ 管道引擎核心
├── scMax_db_manager.py    # 📦 数据库入库主管
├── config_template.yaml   # 唯一参数蓝图
│
├── 01_scQC/               # 数据质控
├── 02_filter/             # 数据过滤
├── 03_cluster/            # 全局聚类
├── 04_subcluster/         # 亚群细分
├── 05_celltype/           # 注释与结题报告
├── 06_differential/       # 差异分析脚本目录（当前实现仍复用此目录名）
├── config_template_06_merge_sub.yaml  # Step 06 配置模板
├── config_template_07_diff.yaml       # Step 07 配置模板
│
├── WebUI/                 # 可视化看板系统 (Flask)
│   ├── scMax_projects.db  # 项目元数据库
│   └── assets/           # 资产隔离物理留底盘
└── jupyter_templates/     # 构成动态报告的 IPYNB 插件库
```

---

## 💡 进阶 Tips
- **智能软链机制**: 在 `05_celltype` 阶段，如果上游已存在预计算好的绘图数据，系统会自动建立软链以节约存储与算力。
- **环境隔离**: 建议使用 `conda` 管理 R 4.2+ 和 Python 3.8+ 的交叉依赖。
- **经验反哺**: 利用终端日志提取出的显著基因，结合数据库历史数据，可以极大加速您的细胞类型定义过程。

---
_Built with ❤️ for Single-Cell Researchers._
