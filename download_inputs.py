"""Download the eight analysis inputs and verify their recorded SHA-256 checksums."""
from pathlib import Path
import argparse
import csv
import hashlib
import urllib.request

ROOT = Path(__file__).resolve().parent

def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    with (ROOT / 'data/input_manifest.csv').open(encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        path = (ROOT / row['path']).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('Input destination is outside the repository')
        if not path.exists() or digest(path) != row['sha256']:
            if args.check_only:
                raise RuntimeError(f"Missing or mismatched input: {row['path']}")
            path.parent.mkdir(parents=True, exist_ok=True)
            partial = path.with_suffix(path.suffix + '.part')
            try:
                request = urllib.request.Request(row['source_url'], headers={'User-Agent': 'PFTeDA-reproducibility/1.2.0'})
                with urllib.request.urlopen(request, timeout=120) as source, partial.open('wb') as target:
                    while block := source.read(1024 * 1024):
                        target.write(block)
                if digest(partial) != row['sha256']:
                    raise RuntimeError(f"Source checksum has changed: {row['path']}")
                partial.replace(path)
            finally:
                if partial.exists():
                    partial.unlink()
        print('Verified:', row['path'])
    print(f'All {len(rows)} inputs verified.')

if __name__ == '__main__':
    main()
