"""Render the public Markdown submission files as standalone PDFs.

Requires ReportLab (the Codex bundled document runtime provides it), DejaVu
Sans fonts, and optionally Poppler for external visual inspection. No network,
patient-scale input, or external publishing is used. Source Markdown stays
authoritative; run this script again after edits.
"""
from __future__ import annotations

import argparse
from html import escape
import os
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FONT_DIR = Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype'
NAVY = colors.HexColor('#123047')
TEAL = colors.HexColor('#0c656b')
GRAY = colors.HexColor('#56626d')
LIGHT = colors.HexColor('#edf3f5')
LINK = re.compile(r'\[([^\]]+)\]\(([^\s]+)\)')


def fonts(directory: Path) -> None:
    for face, suffix in [('DV', ''), ('DV-Bold', '-Bold'), ('DV-Italic', '-Oblique'), ('DV-BoldItalic', '-BoldOblique')]:
        pdfmetrics.registerFont(TTFont(face, str(directory / f'DejaVuSans{suffix}.ttf')))
    pdfmetrics.registerFontFamily('DV', normal='DV', bold='DV-Bold', italic='DV-Italic', boldItalic='DV-BoldItalic')


def plain_inline(text: str) -> str:
    # ReportLab handles these glyphs with the embedded font. ASCII hyphens are
    # used for dashes so copy/paste and PDF search remain predictable.
    text = escape(text.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', ' - '))
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<i>\1</i>', text)
    text = re.sub(r'`([^`]+)`', r'<font color="#174a5b">\1</font>', text)
    return text


def inline(text: str, source: Path) -> str:
    out, last = [], 0
    for match in LINK.finditer(text):
        out.append(plain_inline(text[last:match.start()]))
        label, target = match.groups()
        if target.startswith(('https://', 'http://')):
            out.append(f'<link href="{escape(target, quote=True)}" color="#0c656b">{plain_inline(label)}</link>')
        else:
            resolved = (source.parent / target.split('#')[0]).resolve()
            try:
                relative = resolved.relative_to(ROOT).as_posix()
            except ValueError:
                relative = target
            # A standalone submission cannot depend on file:// links on the
            # author's computer. Retain companion paths as readable text.
            out.append(f'{plain_inline(label)} <font size="7.1" color="#56626d">[{escape(relative)}]</font>')
        last = match.end()
    out.append(plain_inline(text[last:]))
    return ''.join(out)


def styles(compact: bool):
    body_size, leading = (9.0, 11.25) if compact else (9.35, 13.2)
    base = getSampleStyleSheet()
    return {
        'body': ParagraphStyle('BodyDV', fontName='DV', fontSize=body_size, leading=leading, textColor=NAVY, spaceAfter=5 if compact else 7, alignment=TA_LEFT, splitLongWords=True),
        'h1': ParagraphStyle('TitleDV', fontName='DV-Bold', fontSize=16 if compact else 21, leading=19 if compact else 25, textColor=NAVY, spaceAfter=9 if compact else 13, keepWithNext=True),
        'h2': ParagraphStyle('HeadingDV', fontName='DV-Bold', fontSize=12.8 if compact else 13.2, leading=17, textColor=TEAL, spaceBefore=9, spaceAfter=5, keepWithNext=True),
        'h3': ParagraphStyle('SubheadingDV', fontName='DV-Bold', fontSize=10.4, leading=14.2, textColor=TEAL, spaceBefore=6, spaceAfter=4, keepWithNext=True),
        'cell': ParagraphStyle('CellDV', fontName='DV', fontSize=8.2 if compact else 8.8, leading=10.3 if compact else 12.2, textColor=NAVY, spaceAfter=0, splitLongWords=True),
        'cell_head': ParagraphStyle('CellHeadDV', fontName='DV-Bold', fontSize=8.3 if compact else 9, leading=10.4 if compact else 12.3, textColor=colors.white),
        'caption': ParagraphStyle('CaptionDV', fontName='DV-Italic', fontSize=8, leading=10.8, textColor=GRAY, spaceAfter=9),
    }


def table_rows(lines: list[str], source: Path, st, width: float, compact: bool):
    rows = []
    for line in lines:
        fields = [x.strip() for x in line.strip().strip('|').split('|')]
        if all(re.fullmatch(r':?-+:?', x.replace(' ', '')) for x in fields):
            continue
        style = st['cell_head'] if not rows else st['cell']
        rows.append([Paragraph(inline(x, source), style) for x in fields])
    columns = len(rows[0])
    if columns == 2:
        first = width * (0.185 if compact else 0.28)
        widths = [first, width - first]
    else:
        widths = [width / columns] * columns
    table = Table(rows, colWidths=widths, hAlign='LEFT', repeatRows=1)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), TEAL),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, LIGHT]),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 6),
        ('RIGHTPADDING', (0, 0), (-1, -1), 6),
        ('TOPPADDING', (0, 0), (-1, -1), 4 if compact else 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4 if compact else 5),
        ('LINEBELOW', (0, 0), (-1, 0), 0.5, TEAL),
    ]))
    return table


def story(text: str, source: Path, width: float, compact: bool):
    st = styles(compact)
    lines = text.splitlines()
    flowables, index = [], 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line.startswith('|'):
            table_lines = []
            while index < len(lines) and lines[index].strip().startswith('|'):
                table_lines.append(lines[index])
                index += 1
            flowables.extend([table_rows(table_lines, source, st, width, compact), Spacer(1, 6)])
            continue
        image = re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)', line)
        if image:
            path = (source.parent / image.group(2)).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            picture = Image(str(path))
            scale = min(width / picture.imageWidth, 245 / picture.imageHeight)
            picture.drawWidth = picture.imageWidth * scale
            picture.drawHeight = picture.imageHeight * scale
            picture.hAlign = 'LEFT'
            flowables.extend([picture, Paragraph(plain_inline(image.group(1)), st['caption'])])
            index += 1
            continue
        heading = re.match(r'^(#{1,3})\s+(.+)', line)
        if heading:
            flowables.append(Paragraph(inline(heading.group(2), source), st[f'h{len(heading.group(1))}']))
            index += 1
            continue
        if line.startswith('```'):
            block = []
            index += 1
            while index < len(lines) and not lines[index].strip().startswith('```'):
                block.append(lines[index])
                index += 1
            flowables.append(Paragraph('<br/>'.join(escape(x) for x in block), st['body']))
            index += 1
            continue
        paragraph = [line]
        index += 1
        while index < len(lines) and lines[index].strip() and not re.match(r'^(#{1,3}\s|\||!\[|```)', lines[index].strip()):
            # Markdown list items are separate paragraphs, not one long run.
            if re.match(r'^(?:[-*]|\d+\.)\s', lines[index].strip()):
                break
            paragraph.append(lines[index].strip())
            index += 1
        rendered = inline(' '.join(paragraph), source)
        if rendered.startswith('- '):
            rendered = '• ' + rendered[2:]
        flowables.append(Paragraph(rendered, st['body']))
    return flowables


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#cad7dd'))
    canvas.line(doc.leftMargin, 29, A4[0] - doc.rightMargin, 29)
    canvas.setFont('DV', 7.1)
    canvas.setFillColor(GRAY)
    canvas.drawString(doc.leftMargin, 17, 'carcango | MVA Hackathon 2026 | Research hypothesis')
    canvas.drawRightString(A4[0] - doc.rightMargin, 17, str(doc.page))
    canvas.restoreState()


def render(source: Path, output: Path, compact: bool = False, extra: str = ''):
    margin = 28 if compact else 42
    doc = SimpleDocTemplate(str(output), pagesize=A4, rightMargin=margin, leftMargin=margin,
                            topMargin=24 if compact else 33, bottomMargin=39,
                            title=source.read_text().splitlines()[0].lstrip('# '), author='Carson Olafson / carcango',
                            subject='MVA Hackathon 2026 research submission', pageCompression=1)
    text = source.read_text() + extra
    doc.build(story(text, source, A4[0] - 2 * margin, compact), onFirstPage=footer, onLaterPages=footer)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font-dir', type=Path, default=Path(os.environ.get('MVA_PDF_FONT_DIR', DEFAULT_FONT_DIR)))
    args = parser.parse_args()
    fonts(args.font_dir)
    report = ROOT / 'report/track2_report.md'
    acknowledgement = report.read_text().split('## Acknowledgement\n\n', 1)[1].split('\n## ', 1)[0].strip()
    track1_extra = '\n\n## AI use\n\nOpenAI; Pro; model training and data sharing disabled, as reported by Carson. No genotype-scale input is included in this methods document.\n\n## Acknowledgement\n\n' + acknowledgement + '\n'
    outputs = [
        (report, ROOT / 'report/carcango_track2_report.pdf', False, ''),
        (ROOT / 'results/codex-summary.md', ROOT / 'results/carcango_codex_summary.pdf', True, ''),
        (ROOT / 'results/track1_methods.md', ROOT / 'report/carcango_track1_methods.pdf', False, track1_extra),
    ]
    for source, output, compact, extra in outputs:
        render(source, output, compact, extra)
        print(output.relative_to(ROOT))


if __name__ == '__main__':
    main()
