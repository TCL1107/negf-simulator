# preview_readme.py
# 將 README.md 轉成 HTML，並在瀏覽器中預覽

import pathlib
import webbrowser

try:
    import markdown
except ImportError:
    raise SystemExit(
        "請先安裝 markdown 套件：\n\n    pip install markdown\n"
    )

def main():
    root = pathlib.Path(__file__).resolve().parent
    md_path = root / "README.md"
    html_path = root / "README_preview.html"

    if not md_path.exists():
        raise SystemExit(f"找不到 README.md：{md_path}")

    # 讀取 README.md
    md_text = md_path.read_text(encoding="utf-8")

    # 轉成 HTML（啟用常用延伸）
    html_body = markdown.markdown(
        md_text,
        extensions=[
            "extra",        # 支援表格等
            "toc",          # 目錄
            "codehilite",   # 稍微美化 code block
        ]
    )

    # 包一個最基本的 HTML 殼
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>README Preview</title>
<style>
  body {{
    max-width: 900px;
    margin: 40px auto;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    line-height: 1.6;
    padding: 0 16px;
  }}
  pre {{
    background: #f6f8fa;
    padding: 12px;
    overflow-x: auto;
  }}
  code {{
    font-family: "JetBrains Mono", "Fira Code", monospace;
  }}
  img {{
    max-width: 100%;
    height: auto;
  }}
  h1, h2, h3 {{
    border-bottom: 1px solid #eaecef;
    padding-bottom: 0.3em;
  }}
</style>
</head>
<body>
{html_body}
</body>
</html>
"""

    # 寫出 HTML 檔
    html_path.write_text(full_html, encoding="utf-8")
    print(f"已產生：{html_path}")

    # 在瀏覽器打開
    webbrowser.open(html_path.as_uri())


if __name__ == "__main__":
    main()
