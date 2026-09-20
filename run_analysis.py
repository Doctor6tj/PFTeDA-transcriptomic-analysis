"""Run the statistical analysis, evidence synthesis and supplementary-table export."""
from pathlib import Path
import argparse
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rscript', default='Rscript', help='Rscript executable')
    p.add_argument('--output', default='outputs', help='New output directory')
    p.add_argument('--synthesis-only', action='store_true', help='Rebuild evidence and tables from an existing statistical run')
    args = p.parse_args()
    output = Path(args.output).resolve()
    if not args.synthesis_only and output.exists() and any(output.iterdir()):
        raise RuntimeError('Output directory is not empty. Choose a new directory.')
    output.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((ROOT / 'analysis_config.json').read_text(encoding='utf-8'))
    cfg['project_root'] = ROOT.as_posix()
    config = output / 'analysis_config.json'
    config.write_text(json.dumps(cfg, indent=2), encoding='utf-8')
    if not args.synthesis_only:
        subprocess.run([sys.executable, str(ROOT/'download_inputs.py'), '--check-only'], check=True)
        subprocess.run([args.rscript, '--vanilla', str(ROOT/'code/statistics.R'), '--mode=full', f'--config={config}', f'--run-dir={output}'], check=True)
    subprocess.run([sys.executable, str(ROOT/'code/synthesize.py'), '--config', str(config), '--run-dir', str(output), '--evidence-module', str(ROOT/'code/evidence.py')], check=True)
    subprocess.run([sys.executable, str(ROOT/'code/export_tables.py'), '--run-dir', str(output)], check=True)
    subprocess.run([sys.executable, str(ROOT/'verify_results.py'), '--output', str(output)], check=True)

if __name__ == '__main__':
    main()
