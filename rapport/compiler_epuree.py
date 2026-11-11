"""Build the compact presentation with short oral prompts."""
from collections import Counter
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

import pymupdf


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NAME = 'presentation_epuree'
KEEP_PANEL_TITLES = {
    'correlation_01_xi', 'correlation_02_real_fourier', 'correlation_03_relaxation',
}
PANEL_ANALYSIS = (
    'Une variante différente peut changer',
    'Grandes valeurs : ordre et mélange',
    'Une approximation à comparer au calcul',
)


def prepare_figure(name, target):
    if name.startswith('plaquette_'):
        origin = ROOT / 'results/thermal_cycle_L6/visualisation/animations' / (name + '.pdf')
    else:
        origin = ROOT / 'results/production/report/figures' / (name + '.pdf')
    removed = []
    with pymupdf.open(origin) as doc:
        page = doc[0]
        blocks = sorted(page.get_text('blocks'), key=lambda b: b[1])
        top = 0
        if name not in KEEP_PANEL_TITLES:
            top = blocks[0][3] + 3
            removed.append(blocks[0][4].strip())
        if name.startswith('plaquette_'):
            footers = [b for b in blocks if b[1] > page.rect.height - 45]
        else:
            footers = blocks[-1:]
        bottom = min(b[1] for b in footers) - 3
        removed.extend(b[4].strip() for b in footers)
        # The thermodynamic figure places its footnote unusually close to
        # the x-axis labels. Remove the footnote with tight glyph bounds
        # instead of cropping through the labels.
        if name == '03_thermodynamique':
            pymupdf.TOOLS.set_small_glyph_heights(True)
            try:
                tight_footer = sorted(page.get_text('blocks'), key=lambda b: b[1])[-1]
                page.add_redact_annot(pymupdf.Rect(tight_footer[:4]), fill=(1, 1, 1))
                page.apply_redactions(images=0, graphics=0)
            finally:
                pymupdf.TOOLS.set_small_glyph_heights(False)
            bottom = min(page.rect.height, max(b[3] for b in blocks[:-1]) + 3)
        for block in blocks:
            if any(phrase in block[4] for phrase in PANEL_ANALYSIS):
                page.add_redact_annot(pymupdf.Rect(block[:4]), fill=(1, 1, 1))
                removed.append(block[4].strip())
        if page.annots():
            page.apply_redactions(images=0, graphics=0)
        page.set_cropbox(pymupdf.Rect(0, top, page.rect.width, bottom))
        doc.save(target, garbage=4, deflate=True)
    return {
        'file': str(target.relative_to(HERE)), 'source': str(origin.relative_to(ROOT)),
        'source_sha256': hashlib.sha256(origin.read_bytes()).hexdigest(),
        'omitted_text': removed,
    }


def main():
    preserved = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                 for name in ('presentation.tex', 'presentation.pdf')}
    source = (HERE / (NAME + '.tex')).read_text(encoding='utf8')
    # Short oral prompts are allowed; the earlier analysis boxes stay absent.
    for forbidden in (r'\takeaway', r'\notecrit', r'\begin{block}',
                      r'\begin{enumerate}'):
        if forbidden in source:
            raise RuntimeError('Unexpected long-form element: ' + forbidden)
    frames = re.findall(r'\\begin\{frame\}(?:\[[^\]]*\])?(.*?)\\end\{frame\}', source, re.S)
    if len(frames) != 24:
        raise RuntimeError('Expected 24 source frames.')
    for index, frame in enumerate(frames):
        expected_macro = r'\unpoint{' if index == 8 else r'\deuxpoints{'
        other_macro = r'\deuxpoints{' if index == 8 else r'\unpoint{'
        if frame.count(expected_macro) != 1 or other_macro in frame:
            raise RuntimeError(f'Unexpected number of bullet points on slide {index + 1}.')
        if index and frame.count(r'\framesubtitle{') != 1:
            raise RuntimeError(f'Slide {index + 1} must have one explanatory subtitle.')
    names = re.findall(r'\\(?:fpic|gifpic)\{([^}]+)\}', source)
    duplicates = [name for name, count in Counter(names).items() if count > 1]
    if duplicates:
        raise RuntimeError('Repeated figures: ' + repr(duplicates))
    figures = HERE / 'figures_epurees'
    figures.mkdir(exist_ok=True)
    assets = [prepare_figure(name, figures / (name + '.pdf')) for name in names]
    animations = HERE / 'animations'
    animations.mkdir(exist_ok=True)
    for name in ('plaquette_cycle_thermique.gif', 'plaquette_relaxation_isotherme.gif'):
        target = animations / name
        if not target.exists():
            shutil.copy2(ROOT / 'report/animations' / name, target)

    build = HERE / 'build/epuree'
    build.mkdir(parents=True, exist_ok=True)
    executable = shutil.which('pdflatex')
    if not executable:
        raise RuntimeError('pdflatex is required.')
    for iteration in range(3):
        run = subprocess.run(
            [executable, '-interaction=nonstopmode', '-halt-on-error',
             '-output-directory=' + str(build), NAME + '.tex'],
            cwd=HERE, capture_output=True, encoding='utf8', errors='replace',
        )
        (build / f'compile_pass_{iteration + 1}.log').write_text(
            run.stdout + run.stderr, encoding='utf8')
        if run.returncode:
            raise RuntimeError(run.stdout[-4000:])
        log = (build / (NAME + '.log')).read_text(encoding='utf8', errors='replace')
        if iteration >= 1 and 'Rerun' not in log and 'undefined references' not in log:
            break
    warnings = [line for line in log.splitlines()
                if re.search(r'Warning|Overfull|Underfull|undefined references|! LaTeX Error', line)]
    if warnings:
        raise RuntimeError('Layout warnings: ' + repr(warnings))
    with pymupdf.open(build / (NAME + '.pdf')) as doc:
        texts = [page.get_text() for page in doc]
        question_pages = {}
        for question in range(1, 10):
            pattern = re.compile(r'Question\s+' + str(question) + r'\s+[\u2013\u2014-]')
            pages = [i + 1 for i, text in enumerate(texts) if pattern.search(text)]
            if len(pages) != 1:
                raise RuntimeError(f'Question {question}: {pages}')
            question_pages[str(question)] = pages[0]
        if len(doc) != 24:
            raise RuntimeError(f'Expected 24 slides, got {len(doc)}')
        pdf_text = '\n'.join(texts)
        for phrase in ('À retenir', 'Lecture critique', 'Organisation du diaporama'):
            if phrase in pdf_text:
                raise RuntimeError('Unexpected analysis: ' + phrase)
        gif_links = []
        for page in doc:
            for link in page.get_links():
                target = link.get('file', '')
                if target.endswith('.gif'):
                    if not (HERE / target).is_file():
                        raise FileNotFoundError(target)
                    gif_links.append(target)
        if len(set(gif_links)) != 2:
            raise RuntimeError('Both animations must remain accessible.')
        (build / 'extracted.txt').write_text(pdf_text, encoding='utf8')
        preview = build / 'previews'
        preview.mkdir(exist_ok=True)
        for start in range(0, len(doc), 6):
            contact = pymupdf.open()
            page = contact.new_page(width=960, height=810)
            for index in range(start, min(start + 6, len(doc))):
                j = index - start
                page.show_pdf_page(pymupdf.Rect(j % 2 * 480, j // 2 * 270,
                                               (j % 2 + 1) * 480, (j // 2 + 1) * 270), doc, index)
            page.get_pixmap().save(preview / f'contact_{start // 6 + 1}.png')
            contact.close()
        for index in (1, 12, 14, 22):
            doc[index].get_pixmap(matrix=pymupdf.Matrix(1.6, 1.6)).save(
                preview / f'slide_{index + 1:02d}.png')
    for name, digest in preserved.items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Detailed presentation was modified: ' + name)
    shutil.copy2(build / (NAME + '.pdf'), HERE / (NAME + '.pdf'))
    audit = {
        'date': datetime.now().astimezone().isoformat(), 'pages': len(texts),
        'question_pages': question_pages, 'unique_figures': len(names),
        'bullets_per_slide': {str(i): (1 if i == 9 else 2) for i in range(1, 25)},
        'subtitles': '23 frame subtitles plus the cover subtitle',
        'latex_warnings': warnings, 'gif_links': sorted(set(gif_links)),
        'detailed_presentation_preserved': True, 'assets': assets,
    }
    (HERE / 'verification_epuree.json').write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    print(f'Published {len(texts)} slides, {len(names)} figures; no LaTeX warnings.')


if __name__ == '__main__':
    main()
