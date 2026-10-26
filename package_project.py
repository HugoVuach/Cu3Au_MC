#!/usr/bin/env python3
"""Archive de livraison reproductible, sans caches ni campagnes exploratoires."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=base.parent/'Cu3Au_MonteCarlo_complet.zip')
    args = parser.parse_args()
    files = []
    for p in sorted(base.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(base)
        if any(part in {'__pycache__', '.venv', '.git'} for part in rel.parts):
            continue
        if p.suffix in {'.pyc', '.nbc', '.nbi', '.aux', '.toc', '.out'}:
            continue
        if rel.parts[0] == 'results' and len(rel.parts) > 1:
            if rel.parts[1] not in {'production', 'production_run.log', 'production_pipeline.log'}:
                continue
        if rel.name == 'MANIFEST_SHA256.json':
            continue
        files.append(p)
    hashes = {str(p.relative_to(base)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in files}
    manifest = base/'MANIFEST_SHA256.json'
    manifest.write_text(json.dumps(hashes, indent=2)+'\n', encoding='utf8')
    files.append(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in files:
            z.write(p, Path(base.name)/p.relative_to(base))
    with zipfile.ZipFile(args.output) as z:
        bad = z.testzip()
        if bad:
            raise RuntimeError(f'Archive invalide : {bad}')
    print(f'{args.output.resolve()} | {len(files)} fichiers | {args.output.stat().st_size/1e6:.2f} Mo')


if __name__ == '__main__':
    main()
