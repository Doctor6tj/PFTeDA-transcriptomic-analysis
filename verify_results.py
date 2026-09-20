"""Compare numerical results and all supplementary tables with the reference outputs."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent

def compare(expected,actual):
    a=pd.read_csv(expected);b=pd.read_csv(actual)
    pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=False,rtol=1e-8,atol=1e-10)
    maximum=0.0
    for c in a.select_dtypes(include='number').columns:
        differences=np.abs(a[c].to_numpy(dtype=float)-b[c].to_numpy(dtype=float))
        if np.isfinite(differences).any():maximum=max(maximum,float(np.nanmax(differences)))
    return {'rows':len(a),'columns':len(a.columns),'maximum_absolute_numeric_difference':maximum}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='outputs')
    output=Path(p.parse_args().output)
    checks={'ALL_PATHWAY_RESULTS.csv':compare(ROOT/'reference/ALL_PATHWAY_RESULTS.csv',output/'pathway_results/ALL_PATHWAY_RESULTS.csv')}
    for table in sorted((ROOT/'reference/tables').glob('Table_S*.csv')):
        checks[table.name]=compare(table,output/'tables'/table.name)
    report={'status':'PASS','relative_tolerance':1e-8,'absolute_tolerance':1e-10,'comparisons':checks,'study_design_table':'S1 is copied from the documented study design; S2-S14 are generated from the analysis.'}
    (output/'result_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(f'PASS: {len(checks)} result tables match; categorical values are unchanged.')

if __name__=='__main__':main()
