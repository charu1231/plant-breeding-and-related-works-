#!/usr/bin/env python3
"""05_make_pdf.py — render paper/PAPER.md into paper/PAPER.pdf (A4, DejaVu fonts,
embedded 6 figures, styled tables) via markdown -> HTML -> fpdf2."""
import os
import re
from pathlib import Path
import markdown
import matplotlib
from fpdf import FPDF

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "paper" / "PAPER.md"
OUT = ROOT / "paper" / "PAPER.pdf"

text = SRC.read_text(encoding="utf-8")

# ---- front matter -> title block
fm = {}
if text.startswith("---"):
    end = text.index("\n---", 3)
    block = text[3:end]
    key = None
    for ln in block.splitlines():
        m = re.match(r"^(\w[\w-]*):\s*(.*)$", ln)
        if m:
            key, val = m.group(1), m.group(2).strip().strip('"')
            fm[key] = val
        elif key and ln.startswith(" "):
            fm[key] += " " + ln.strip().strip('"')
    body_md = text[end + 4:]
else:
    body_md = text

header_html = f"""
<p align="center"><font size="17"><b>{fm.get('title','').replace('|','').strip()}</b></font></p>
<p align="center"><font size="11">{fm.get('authors','')}</font></p>
<p align="center"><font size="10">{fm.get('date','')} &nbsp;|&nbsp; <i>{fm.get('status','')}</i></font></p>
<hr/>
"""

# escape special chars that break the HTML parser (keep markdown syntax intact)
body_md = body_md.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

html = markdown.markdown(body_md, extensions=["tables", "sane_lists"])

# avoid core-only "courier" for <code>/<pre> (no unicode): render in DejaVu
html = html.replace("<code>", '<font size="9">').replace("</code>", "</font>")
html = html.replace("<pre>", '<p><font size="9">').replace("</pre>", "</font></p>")

# table styling: borders + small font
html = html.replace("<table>", '<table border="1" align="center" width="100%">')
html = re.sub(r"<th>", '<th align="center"><font size="8"><b>', html)
html = re.sub(r"</th>", "</b></font></th>", html)
html = re.sub(r"<td>", '<td align="center"><font size="8">', html)
html = re.sub(r"</td>", "</font></td>", html)

# image sizing + centering
def fix_img(m):
    alt, src = m.group(1), m.group(2)
    return f'<p align="center"><img alt="{alt}" src="{src}" width="660"/></p>'
html = re.sub(r'<img alt="([^"]*)" src="([^"]+)"\s*/?>', fix_img, html)

full_html = header_html + html

# ---- build PDF
FDIR = Path(os.path.dirname(matplotlib.__file__)) / "mpl-data" / "fonts" / "ttf"
pdf = FPDF(orientation="P", unit="mm", format="A4")
pdf.set_margins(17, 16, 17)
pdf.set_auto_page_break(True, margin=18)
pdf.add_font("djv", "", str(FDIR / "DejaVuSans.ttf"))
pdf.add_font("djv", "B", str(FDIR / "DejaVuSans-Bold.ttf"))
pdf.add_font("djv", "I", str(FDIR / "DejaVuSans-Oblique.ttf"))
pdf.add_font("djv", "BI", str(FDIR / "DejaVuSans-BoldOblique.ttf"))
pdf.set_font("djv", size=10.5)
pdf.add_page()

os.chdir(ROOT / "paper")   # resolve ../results/figures/*.png
pdf.write_html(full_html,
               font_family="djv",
               heading_sizes={"h1": 15, "h2": 12.5, "h3": 11.5},
               table_line_separators=False)
pdf.output(str(OUT))
print("wrote", OUT, f"({OUT.stat().st_size/1e6:.2f} MB)")
