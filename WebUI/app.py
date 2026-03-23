from flask import Flask, render_template, request, g, send_file
import sqlite3
import os

app = Flask(__name__)
app.config['SECRET_KEY'] = 'scmax_secret_ui'

# 配置文件中定义的数据库路径
DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scMax_projects.db')

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DB_PATH)
        db.row_factory = sqlite3.Row
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
    
    # 增加仪表盘大盘概况统计
    cursor.execute('SELECT COUNT(DISTINCT custom_project_id), COUNT(DISTINCT species), COUNT(DISTINCT tissue_type) FROM projects')
    proj_row = cursor.fetchone()
    
    cursor.execute('SELECT COUNT(DISTINCT cell_type) FROM marker_experience')
    cell_row = cursor.fetchone()
    
    stats = {
        'projects': proj_row[0] or 0,
        'species': proj_row[1] or 0,
        'tissues': proj_row[2] or 0,
        'cell_types': cell_row[0] or 0
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
    
    if bk_path and os.path.exists(bk_path):
        return send_file(bk_path)
    if path and os.path.exists(path):
        return send_file(path)
    return "Image not found (Source and backup may be deleted)", 404

@app.route('/view_report')
def view_report():
    """提供本地绝对路径HTML的代理挂载 (支持备份图库降级)"""
    path = request.args.get('path')
    bk_path = request.args.get('bk_path')
    
    if bk_path and os.path.exists(bk_path):
        return send_file(bk_path)
    if path and os.path.exists(path):
        return send_file(path)
    return "Report not found (Source and backup may be deleted)", 404

if __name__ == '__main__':
    # 供本地直接调试用，生产可由 wsgi 启动
    app.run(host='0.0.0.0', port=5050, debug=True)
