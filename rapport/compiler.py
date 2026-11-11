"""Compile the local Beamer presentation and verify the resulting document.

Run from any directory with the project's system Python (PyMuPDF installed):
    python rapport/compiler.py
Figures and animations are refreshed from the archived report, without MC runs.
"""
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


def main():
    source = (HERE / 'presentation.tex').read_text(encoding='utf8')
    figure_names = re.findall(r'\\fpic\{([^}]+)\}', source)
    figures = HERE / 'figures'
    figures.mkdir(exist_ok=True)
    assets = []
    for name in dict.fromkeys(figure_names):
        filename = name + '.pdf'
        if name.startswith('plaquette_'):
            origin = ROOT / 'results/thermal_cycle_L6/visualisation/animations' / filename
        else:
            origin = ROOT / 'results/production/report/figures' / filename
        destination = figures / filename
        if not origin.is_file():
            raise FileNotFoundError(origin)
        shutil.copy2(origin, destination)
        assets.append({
            'file': str(destination.relative_to(HERE)),
            'source': str(origin.relative_to(ROOT)),
            'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
        })

    animations = HERE / 'animations'
    animations.mkdir(exist_ok=True)
    for filename in ('plaquette_cycle_thermique.gif', 'plaquette_relaxation_isotherme.gif'):
        shutil.copy2(ROOT / 'report/animations' / filename, animations / filename)

    executable = shutil.which('pdflatex')
    if not executable:
        raise RuntimeError('pdflatex is required to compile the presentation.')
    build = HERE / 'build'
    build.mkdir(exist_ok=True)
    for iteration in range(3):
        run = subprocess.run(
            [executable, '-interaction=nonstopmode', '-halt-on-error',
             '-output-directory=' + str(build), 'presentation.tex'],
            cwd=HERE, capture_output=True, encoding='utf8', errors='replace',
        )
        (build / f'compile_pass_{iteration + 1}.log').write_text(
            run.stdout + run.stderr, encoding='utf8')
        if run.returncode:
            raise RuntimeError(run.stdout[-4000:])
        log = (build / 'presentation.log').read_text(encoding='utf8', errors='replace')
        if iteration >= 1 and 'Rerun' not in log and 'undefined references' not in log:
            break

    pdf = build / 'presentation.pdf'
    warnings = [line for line in log.splitlines()
                if re.search(r'Warning|Overfull|Underfull|undefined references|! LaTeX Error', line)]
    with pymupdf.open(pdf) as doc:
        texts = [page.get_text() for page in doc]
        question_pages = {}
        for question in range(1, 10):
            pattern = re.compile(r'Question\s+' + str(question) + r'\s+[\u2013\u2014-]')
            pages = [i + 1 for i, text in enumerate(texts) if pattern.search(text)]
            expected_count = 2 if question == 1 else 1
            if len(pages) != expected_count:
                raise RuntimeError(f'Question {question}: expected {expected_count} pages, found {pages}')
            question_pages[str(question)] = pages
        if len(doc) != 25:
            raise RuntimeError(f'Expected 25 slides, found {len(doc)}')
        if question_pages['1'] != [2, 3] or question_pages['4'] != [25]:
            raise RuntimeError('Expected two opening Q1 slides and Q4 as the final appendix.')
        missing_links = []
        gif_links = []
        for page in doc:
            for link in page.get_links():
                filename = link.get('file', '')
                if filename.endswith('.gif'):
                    gif_links.append(filename)
                    if not (HERE / filename).exists():
                        missing_links.append(filename)
        if missing_links:
            raise RuntimeError(f'Missing animation links: {missing_links}')
        dimensions = [round(value * 25.4 / 72, 2)
                      for value in (doc[0].rect.width, doc[0].rect.height)]
        (build / 'presentation_extracted.txt').write_text('\n\n'.join(texts), encoding='utf8')
        audit = {
            'date': datetime.now().astimezone().isoformat(),
            'pages': len(doc), 'slide_dimensions_mm': dimensions,
            'question_pages': question_pages, 'main_slides_including_cover': 10,
            'bonus_appendix_slides': 14, 'validation_appendix_slides': 1,
            'figure_uses': len(figure_names), 'unique_figures': len(set(figure_names)),
            'repeated_figure_names': [name for name, n in Counter(figure_names).items() if n > 1],
            'latex_warnings': warnings, 'gif_links': gif_links, 'assets': assets,
        }
    (HERE / 'verification.json').write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    if warnings:
        raise RuntimeError('Layout or LaTeX warnings require review: ' + repr(warnings))
    shutil.copy2(pdf, HERE / 'presentation.pdf')
    print(f'Published {len(texts)} slides, {len(set(figure_names))} figures; no LaTeX warnings.')


if __name__ == '__main__':
    main()
