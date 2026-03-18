# img2singlehtml.py
import base64
import re
from pathlib import Path

def convert_images_to_base64(html_file):
    html_path = Path(html_file).resolve()
    html_content = html_path.read_text()
    
    # 匹配所有图片标签
    pattern = re.compile(r'<img\s+[^>]*src="([^"]+\.(?:png|jpg|jpeg|gif))"[^>]*>', re.IGNORECASE)
    
    for match in set(pattern.findall(html_content)):  # 去重处理
        img_path = html_path.parent / match
        if img_path.exists():
            mime_type = f"image/{img_path.suffix[1:].lower()}"
            b64_data = base64.b64encode(img_path.read_bytes()).decode('utf-8')
            html_content = html_content.replace(
                f'src="{match}"', 
                f'src="data:{mime_type};base64,{b64_data}"'
            )
    
    output_path = html_path.parent / f"CellType.Annotation_{html_path.name}"
    output_path.write_text(html_content)
    print(f"生成自包含文件: {output_path}")

if __name__ == "__main__":
    import sys
    convert_images_to_base64(sys.argv[1])
