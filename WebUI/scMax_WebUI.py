import os
import sqlite3
from flask import Flask, render_template_string, send_file, request, abort

app = Flask(__name__)

DB_PATH = os.environ.get("SCMAX_DB_PATH", os.path.join(os.path.dirname(__file__), "scMax_projects.db"))
FILE_BASE = os.environ.get("SCMAX_FILE_BASE", "")

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="zh">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>scMax 单细胞项目数据库总屏</title>
    <!-- 使用高度兼容的 Bootstrap 和 DataTables CDN -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/dataTables.bootstrap5.min.css">
    <style>
        body { background-color: #f4f6f9; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
        .navbar { background: linear-gradient(135deg, #1d2b64 0%, #f8cdda 100%); }
        .card { border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border: none; }
        .table-hover tbody tr:hover { background-color: #e9ecef; }
    </style>
</head>
<body>
    <nav class="navbar navbar-dark shadow-sm">
        <div class="container-fluid">
            <span class="navbar-brand mb-0 h1 fw-bold text-white">✨ scMax 智能聚类分析与数据库中央大屏 ✨</span>
        </div>
    </nav>
    
    <div class="container-fluid mt-4 px-4">
        <div class="card p-4">
            <table id="projectsTable" class="table table-striped table-hover align-middle">
                <thead class="table-dark">
                    <tr>
                        <th>入库时间</th>
                        <th>标分编号</th>
                        <th>自定义编号</th>
                        <th>物种</th>
        <th>注释 UMAP</th>
        <th>细胞比例分布图</th>
        <th>细胞注释表格</th>
        <th>质控与大屏分析报告 (HTML)</th>
                        <th>项目本地存放路径</th>
                        <th>备注</th>
                    </tr>
                </thead>
                <tbody>
                    {% for row in rows %}
                    <tr>
                        <td>{{ row.created_at }}</td>
                        <td><span class="badge bg-primary fs-6">{{ row.std_project_id }}</span></td>
                        <td><span class="badge bg-secondary fs-6">{{ row.custom_project_id }}</span></td>
                        <td>{{ row.species }}</td>
                        <td>
                            {% if row.umap_path and row.umap_path != "NA" %}
                                <a href="{{ url_for('serve_file', path=row.umap_path) }}" target="_blank">
                                    <img src="{{ url_for('serve_file', path=row.umap_path) }}" alt="UMAP" style="max-width:120px; height:auto;">
                                </a>
                            {% else %}
                                <span class="text-muted">无</span>
                            {% endif %}
                        </td>
                        <td>
                            {% if row.fraction_path and row.fraction_path != "NA" %}
                                <a href="{{ url_for('serve_file', path=row.fraction_path) }}" target="_blank">
                                    <img src="{{ url_for('serve_file', path=row.fraction_path) }}" alt="Fraction" style="max-width:120px; height:auto;">
                                </a>
                            {% else %}
                                <span class="text-muted">无</span>
                            {% endif %}
                        </td>
                        <td>
                            {% if row.anno_table_path and row.anno_table_path != "NA" %}
                                <a href="{{ url_for('serve_file', path=row.anno_table_path) }}" target="_blank" class="btn btn-sm btn-outline-primary">
                                    下载/查看
                                </a>
                            {% else %}
                                <span class="text-muted">无</span>
                            {% endif %}
                        </td>
                        <td>
                            {% if row.anno_html_path and row.anno_html_path != "NA" %}
                                <a href="{{ url_for('serve_file', path=row.anno_html_path) }}" target="_blank" class="btn btn-sm btn-outline-success">
                                    在浏览器展开详细报告
                                </a>
                            {% else %}
                                <span class="badge bg-secondary">未制备报告</span>
                            {% endif %}
                        </td>
                        <td class="text-break" style="max-width:250px;"><code>{{ row.result_path }}</code></td>
                        <td>{{ row.notes }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
    </div>

    <!-- 必要的交互 JS -->
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/dataTables.bootstrap5.min.js"></script>
    <script>
        $(document).ready(function() {
            $('#projectsTable').DataTable({
                "language": {
                    "url": "https://cdn.datatables.net/plug-ins/1.13.6/i18n/zh.json"
                },
                "order": [[ 0, "desc" ]],
                "pageLength": 10
            });
        });
    </script>
</body>
</html>
"""

def init_db():
    conn = sqlite3.connect(DB_PATH)
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
    conn.commit()
    conn.close()

def _resolve_path(raw_path: str) -> str:
    if not raw_path:
        return ""
    path = os.path.expanduser(raw_path)
    if os.path.isabs(path):
        abs_path = path
    else:
        base = FILE_BASE or os.getcwd()
        abs_path = os.path.abspath(os.path.join(base, path))
    if FILE_BASE:
        base = os.path.abspath(FILE_BASE)
        if not (abs_path == base or abs_path.startswith(base + os.sep)):
            return ""
    return abs_path

@app.route('/file')
def serve_file():
    raw_path = request.args.get("path", "")
    path = _resolve_path(raw_path)
    if not path or not os.path.exists(path):
        abort(404)
    return send_file(path, as_attachment=False)

@app.route('/')
def index():
    init_db()
    
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    c.execute("SELECT * FROM projects ORDER BY id DESC")
    rows = c.fetchall()
    conn.close()
    
    return render_template_string(HTML_TEMPLATE, rows=rows)

if __name__ == '__main__':
    print(f"[*] 正在启动 scMax WebUI 服务器 ...")
    print(f"[*] 监听数据库: {DB_PATH}")
    print(f"[*] 建议使用 Nginx/Gunicorn 包裹在生产环境，本地测试在浏览器打开 http://0.0.0.0:8080")
    init_db()
    app.run(host='0.0.0.0', port=8080, debug=True)
