from flask import Flask, render_template, request, g, send_file
import sqlite3
import os
from urllib.parse import urlparse, parse_qs

app = Flask(__name__)
app.config['SECRET_KEY'] = 'scmax_secret_ui'

# 配置文件中定义的数据库路径
DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scMax_projects.db')

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DB_PATH)
        db.row_factory = sqlite3.Row
        
        # 兼容旧数据库，补充新字段
        cursor = db.cursor()
        for col in ['disease_type', 'seq_tech', 'anno_png_path', 'bk_anno_png']:
            try:
                cursor.execute(f"ALTER TABLE projects ADD COLUMN {col} TEXT")
            except sqlite3.OperationalError:
                pass
        db.commit()
    return db

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

@app.route('/')
def index():
    """主项目记录面板"""
    db = get_db()
    cursor = db.cursor()
    cursor.execute('''
        SELECT * FROM projects 
        ORDER BY created_at DESC 
        LIMIT 100
    ''')
    projects = cursor.fetchall()
    
    # 增加仪表盘大盘概况统计（分布汇总，移除上限，展示全库种类）
    cursor.execute('SELECT species_zh, COUNT(*) as cnt FROM projects WHERE species_zh IS NOT NULL AND species_zh != "" GROUP BY species_zh ORDER BY cnt DESC')
    species_dist = cursor.fetchall()
    
    cursor.execute('SELECT tissue_type_zh, COUNT(*) as cnt FROM projects WHERE tissue_type_zh IS NOT NULL AND tissue_type_zh != "" GROUP BY tissue_type_zh ORDER BY cnt DESC')
    tissue_dist = cursor.fetchall()
    
    # 因为 disease_type 是新加字段，旧数据可能为 NULL
    cursor.execute('SELECT disease_type, COUNT(*) as cnt FROM projects WHERE disease_type IS NOT NULL AND disease_type != "" GROUP BY disease_type ORDER BY cnt DESC')
    disease_dist = cursor.fetchall()

    cursor.execute('SELECT COUNT(DISTINCT custom_project_id), COUNT(DISTINCT species), COUNT(DISTINCT tissue_type) FROM projects')
    proj_row = cursor.fetchone()
    
    cursor.execute('SELECT COUNT(DISTINCT cell_type) FROM marker_experience')
    cell_row = cursor.fetchone()
    
    stats = {
        'projects': proj_row[0] or 0,
        'species': proj_row[1] or 0,
        'tissues': proj_row[2] or 0,
        'cell_types': cell_row[0] or 0,
        'species_dist': [{'name': r[0], 'cnt': r[1]} for r in species_dist],
        'tissue_dist': [{'name': r[0], 'cnt': r[1]} for r in tissue_dist],
        'disease_dist': [{'name': r[0], 'cnt': r[1]} for r in disease_dist]
    }
    
    return render_template('index.html', projects=projects, stats=stats, active_page='projects')

@app.route('/experience')
def experience():
    """经验库查询面板"""
    db = get_db()
    cursor = db.cursor()
    cursor.execute('''
        SELECT * FROM marker_experience 
        ORDER BY created_at DESC 
        LIMIT 500
    ''')
    experiences = cursor.fetchall()
    return render_template('experience.html', experiences=experiences, active_page='experience')

@app.route('/view_image')
def view_image():
    """提供本地绝对路径图片的静态代理 (支持备份图库降级)"""
    path = request.args.get('path')
    bk_path = request.args.get('bk_path')
    
    if path and os.path.exists(path):
        return send_file(path)
    if bk_path and os.path.exists(bk_path):
        return send_file(bk_path)
    return "Image not found (Source and backup may be deleted)", 404

@app.route('/view_report')
def view_report():
    """提供本地绝对路径HTML的代理挂载 (支持备份图库降级)"""
    path = request.args.get('path')
    bk_path = request.args.get('bk_path')
    
    if path and os.path.exists(path):
        return send_file(path)
    if bk_path and os.path.exists(bk_path):
        return send_file(bk_path)
    return "Report not found (Source and backup may be deleted)", 404

@app.route('/<path:filename>')
def serve_report_assets(filename):
    """
    用于拦截网页报告（HTML）内部加载本地相对路径图片时引发的 404 请求。
    通过读取 HTTP Referer 头来溯源当前正在打开的报告文件位置，并拼接触发该请求的静态图片路径。
    """
    referer = request.headers.get('Referer')
    if referer:
        parsed_url = urlparse(referer)
        query_params = parse_qs(parsed_url.query)
        for arg in ['path', 'bk_path']:
            if arg in query_params:
                report_path = query_params[arg][0]
                if os.path.exists(report_path):
                    report_dir = os.path.dirname(report_path)
                    asset_path = os.path.join(report_dir, filename)
                    if os.path.exists(asset_path):
                        return send_file(asset_path)
    return "Asset not found", 404

if __name__ == '__main__':
    # 供本地直接调试用，生产可由 wsgi 启动
    app.run(host='0.0.0.0', port=5050, debug=True)
