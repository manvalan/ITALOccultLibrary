#!/usr/bin/env python3
"""
Script per generare il PDF del report di validazione
"""

import markdown
from weasyprint import HTML, CSS
from pathlib import Path

# Percorsi
script_dir = Path(__file__).parent
md_file = script_dir / "VALIDATION_REPORT.md"
pdf_file = script_dir / "VALIDATION_REPORT.pdf"

# Leggi il markdown
with open(md_file, 'r', encoding='utf-8') as f:
    md_content = f.read()

# Converti in HTML
md = markdown.Markdown(extensions=['tables', 'fenced_code', 'codehilite', 'toc'])
html_body = md.convert(md_content)

# CSS per lo styling
css = """
@page {
    size: A4;
    margin: 2cm;
    @top-center {
        content: "AstDyn Library - Validation Report";
        font-size: 9pt;
        color: #666;
    }
    @bottom-center {
        content: counter(page) " / " counter(pages);
        font-size: 9pt;
        color: #666;
    }
}

body {
    font-family: 'Helvetica Neue', Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.5;
    color: #333;
}

h1 {
    color: #1a365d;
    font-size: 24pt;
    border-bottom: 3px solid #1a365d;
    padding-bottom: 10px;
    margin-top: 30px;
}

h2 {
    color: #2c5282;
    font-size: 18pt;
    border-bottom: 2px solid #e2e8f0;
    padding-bottom: 8px;
    margin-top: 25px;
    page-break-after: avoid;
}

h3 {
    color: #2d3748;
    font-size: 14pt;
    margin-top: 20px;
    page-break-after: avoid;
}

h4 {
    color: #4a5568;
    font-size: 12pt;
    margin-top: 15px;
}

table {
    border-collapse: collapse;
    width: 100%;
    margin: 15px 0;
    font-size: 10pt;
    page-break-inside: avoid;
}

th {
    background-color: #2c5282;
    color: white;
    padding: 10px 8px;
    text-align: left;
    font-weight: bold;
}

td {
    padding: 8px;
    border-bottom: 1px solid #e2e8f0;
}

tr:nth-child(even) {
    background-color: #f7fafc;
}

tr:hover {
    background-color: #edf2f7;
}

code {
    font-family: 'Monaco', 'Consolas', monospace;
    font-size: 9pt;
    background-color: #edf2f7;
    padding: 2px 6px;
    border-radius: 3px;
}

pre {
    background-color: #2d3748;
    color: #e2e8f0;
    padding: 15px;
    border-radius: 5px;
    overflow-x: auto;
    font-size: 9pt;
    page-break-inside: avoid;
}

pre code {
    background-color: transparent;
    padding: 0;
    color: inherit;
}

blockquote {
    border-left: 4px solid #4299e1;
    margin: 15px 0;
    padding: 10px 20px;
    background-color: #ebf8ff;
    color: #2c5282;
}

ul, ol {
    margin: 10px 0;
    padding-left: 25px;
}

li {
    margin: 5px 0;
}

strong {
    color: #1a365d;
}

/* Classi speciali per i risultati */
.pass, .success {
    color: #22543d;
    background-color: #c6f6d5;
    padding: 2px 6px;
    border-radius: 3px;
}

.fail, .error {
    color: #742a2a;
    background-color: #fed7d7;
    padding: 2px 6px;
    border-radius: 3px;
}

/* Prima pagina */
.title-page {
    text-align: center;
    padding-top: 100px;
}

/* Footer per ogni sezione */
.section-footer {
    margin-top: 30px;
    border-top: 1px solid #e2e8f0;
    padding-top: 10px;
    font-size: 9pt;
    color: #718096;
}

/* Box per statistiche */
.stats-box {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 20px;
    border-radius: 10px;
    margin: 20px 0;
    text-align: center;
}

/* Emoji support */
.emoji {
    font-size: 14pt;
}
"""

# HTML completo
html_doc = f"""
<!DOCTYPE html>
<html lang="it">
<head>
    <meta charset="UTF-8">
    <title>AstDyn Library - Validation Report</title>
</head>
<body>
{html_body}
</body>
</html>
"""

# Genera PDF
print("Generazione PDF in corso...")
HTML(string=html_doc).write_pdf(
    pdf_file,
    stylesheets=[CSS(string=css)]
)

print(f"✓ PDF generato: {pdf_file}")
print(f"  Dimensione: {pdf_file.stat().st_size / 1024:.1f} KB")
