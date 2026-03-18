# scMax 单细胞转录组分析流水线与 WebUI 智能管理平台 🧬

scMax 是一套端到端的自动化单细胞转录组分析、可视化报告生成及全周期数据沉淀检索工具。
不仅整合了标准的生信降维注释处理，更是首创性地结合了动态图文全自动 HTML 出片引擎以及轻量级的 SQLite + Flask 本地多维经验大盘展示平台，彻底打通了生信跑批到项目信息结构化维护的“最后一公里”。

## ✨ 核心特性

- 🧪 **自动化流水线 (Pipeline)**
  基于 Python `scMax_pipeline.py` 调度各项子流程（包括基础统计、质控、双胞打分、分群降维以及自动注释），通过一份集中的 `config_template.yaml` 灵活接管所有子任务参数。
- 📊 **动态报告组件库 (`generate_dynamic_report.py`)**
  告别手工整理结果页。流程末尾自动接管图表产出，通过组件化拼装架构以及 `jupyter nbconvert` 的无缝转译，一键输出附带排版与动态关联的 `report_custom.html` 交互全景图文报告。
- 📦 **全自动的 SQLite 资产双轨入库 (`scMax_db_manager.py`)**
  内置防重入防护库管理引擎。
  每次执行结尾均会自动识别并更新：
  - **主项目大盘库 (`scMax_projects.db`)**：涵盖项目 ID、CRM 单号、分群层次、分析人员等全要素，自动追踪记录当前靶向产生的 UMAP 图与细胞占比统计表。
  - **知识图谱经验字典 (`scMax_experience.db`)**：精细切粒度提取单细胞转录组生成的 Marker 统计表，细化至种属、脏器以及每群特异性标识性高表达基因，成为长尾数据的算法沉淀集。
- 🤖 **智能辅助注释草图系统 (`generate_anno_draft.py`)**。
  反哺机制的终极体现。内置该脚本可拦截 Seurat 生成的 `markers.csv`，通过筛选（`p_val_adj < 0.05` 且降序提取 Top 30 `avg_log2FC`）最严苛显著的差异基因，连线 SQLite `scMax_experience.db` 经验大库。
  系统能够依据历次录入的真实经验自动推断高分的细胞亚型，将画点图 Markers、参考描述(`description`)与文献(`reference`)秒出成规范合规的 `annofile` 草图（Excel）。专家只需“批改作业”，免去手动填表苦差事。
- 🛡️ **隔离图床缓存机制 (Asset Quarantine Backup)**
  完美解决传统出图目录 `outdir` 被人为误删或腾挪后的数据“失独断链”问题。每次提取入库时皆会通过脚本同步硬拷贝产出的核心图表及报告文本入 `WebUI/assets/` 物理归档隔离舱区，建立最深防线的“双库挂载”。
- 🌐 **高可用 Flask 可视化大盘检索平台 (`WebUI`)**
  通过轻量化的 Python Flask 在本地开启看板，搭载 Bootstrap 四维实时仪表盘（涵盖统计量）、高阶 DataTables 前端筛选搜索组件，以及基于备份路径回退技术的降级灾备代理机制，实现本地长年不断供的单细胞经验数据中心查询门户。

## 📍 架构层级

```text
scMax/
├── README.md                          # 项目说明文件
├── config_template.yaml               # 唯一配置中心
├── scMax_pipeline.py                  # 流水线与分支路逻辑主车间
├── generate_anno_draft.py             # 🌟 经验反哺：基于数据库倒推打标初稿的 AI 脚本
├── generate_dynamic_report.py         # 智能拼装 HTML 单细胞交互报告
├── scMax_db_manager.py                # 跨平台 SQLite 数据库防重落盘总管与资产隔离器
├── .gemini/                           # AI 助理的工作台规划缓存
│
├── jupyter_templates/                 # 构成最终分析报告的可插拔零件模型组件
│   └── *.ipynb
│
├── WebUI/                             # 网页看板主程序包
│   ├── app.py                         # 基于 Flask 搭载的数据中台 Server 与图床灾备处理代理路由
│   ├── scMax_projects.db              # 运行时自动生成的单细胞数据主从 SQLite 关联实体
│   ├── assets/                        # 受防删灾备隔离机制保护的流水线 UMAP / Bar 等图谱硬留底归档盘
│   └── templates/                     # Bootstrap 前端静态交互大盘
│       ├── layout.html
│       ├── index.html
│       └── experience.html
│
└── tmp_code/                          # 流水线跑批底层生信调度引擎及子函数库 (R & Shell 杂烩)
```

## 🚀 启动指南

### 1. 运行核心分析与打库流水线

由于脚本本身包含完整的模块钩子化以及前置防护层，您可以在修改好 `config_template.yaml` 内的项目参数（如输入输出目录和组织物种属性）后直接通过执行管道主入口触发端到端产出：
```bash
python3 scMax_pipeline.py -c config_template.yaml
```
程序将会一次性走完：生信底层解析计算 -> 组件拼装图文 HTML 报告产出 -> 图片分离开垦拷贝 -> SQLite 防重注册归档 四大流程。

### 2. 长青开启网页前端经验面板

如果你需要调阅分析结论或者是通过以往的标记历史经验查询新的细胞集群走向，在主根目录启动独立中台进程：
```bash
python3 WebUI/app.py
```
默认会在本地挂载服务监听端口 `5050`。打开任意浏览器直塞 `http://127.0.0.1:5050` 见证单细胞图谱大数据！

### 3. (高阶) 智能生成注释草图 (Draft Annotation)

如果在第 5 步你苦于不想手工创建专家鉴别表（`annofile`），那可以利用前序产出的集落差异 Marker 表格，使用经验库反打标签：

```bash
python3 generate_anno_draft.py -i <你的_03或04_产物/markers_res0.4.csv> \
                               -c config_template.yaml \
                               -o my_draft_anno.xlsx
```
脚本提取的高频高显著性基因（`p_adj < 0.05` & Top 30 `avg_log2FC`）会和你的历史经验进行比对，光速为你生成填好推断名字和描述文献的雏形 Excel！你核对无误后将其路径放进 YAML 跑第 5 步即可。

---

_Powered By Angela & Antigravity Assistant._
