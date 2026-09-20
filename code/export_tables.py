"""Export Tables S1-S14 using the manuscript's table names and display labels."""
from pathlib import Path
import argparse
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    'primary_common_five': 'Primary common-five analysis',
    'sensitivity_all_available': 'All-available-dose sensitivity analysis',
    'sensitivity_lowest_four': 'Lowest-four-dose sensitivity analysis',
}

def truth(series):
    return series.astype(str).str.lower().eq('true')

def sensitivity_table(run):
    result = pd.read_csv(run/'FAMILY3_GRADE_COMPARISON.csv')
    # Display order follows the supplementary workbook.
    result = result[['pathway', *LABELS]]
    for arm, label in [('primary_common_five','common_five'), ('sensitivity_all_available','all_available'), ('sensitivity_lowest_four','lowest_four')]:
        family = pd.read_csv(run/arm/'FAMILY_PATHWAY_EVIDENCE.csv')
        family = family[family.family_id.eq('FAMILY_3')].copy()
        keep = ['pathway','strong_contexts','weak_contexts','family_vote','family_weak_support','model_restricted','context_states']
        family = family[keep].rename(columns={c:f'family3_{label}_{c}' for c in keep if c != 'pathway'})
        result = result.merge(family,on='pathway',how='left',validate='one_to_one')
    result['reproducibility_grade_varies'] = result[list(LABELS)].nunique(axis=1).gt(1)
    result['dose_definitions'] = 'primary: 0.06, 0.67, 3.35, 6.7, 17 uM at both durations; all-available sensitivity: those doses plus 33 uM at 24 h; lowest-four sensitivity: 0.06, 0.67, 3.35, 6.7 uM at both durations'
    return result.rename(columns=LABELS)

def main():
    p=argparse.ArgumentParser();p.add_argument('--run-dir', required=True)
    run=Path(p.parse_args().run_dir);out=run/'tables';out.mkdir(exist_ok=True)
    def name(n):
        return next((ROOT/'reference/tables').glob(f'Table_S{n}_*.csv')).name
    def write(n, frame):
        frame.to_csv(out/name(n),index=False)
    primary=run/'primary_common_five'
    mapping={3:'CONTEXT_PATHWAY_EVIDENCE.csv',4:'FAMILY_PATHWAY_EVIDENCE.csv',5:'CROSS_FAMILY_REPRODUCIBILITY.csv',6:'R2_R3_REPRODUCIBLE_MODULES.csv',7:'PFTeDA_PFOA_DIRECT_PROGRAMS_ALL.csv',8:'SUPPORTED_PFTeDA_PFOA_DIFFERENTIAL_PROGRAMS.csv',9:'CONTEXT_HETEROGENEITY_ASSESSMENT.csv',10:'SUPPORTED_CONTEXT_HETEROGENEITY.csv'}
    # S1 describes the study design rather than an estimated statistical quantity.
    shutil.copyfile(ROOT/'reference/tables'/name(1),out/name(1))
    write(2,pd.read_csv(run/'qc/MODEL_QC_SUMMARY.csv'))
    for n, source in mapping.items():write(n,pd.read_csv(primary/source))
    ctx=pd.read_csv(primary/mapping[3]);fam=pd.read_csv(primary/mapping[4]);cross=pd.read_csv(primary/mapping[5]);direct=pd.read_csv(primary/mapping[8])
    restricted=fam[truth(fam.model_restricted)].copy()
    restricted['pathway_display']=restricted.pathway
    write(11,restricted)
    cm=int(truth(ctx[ctx.context.eq('iPSC-CM')].strong).sum())
    liver=int(truth(ctx[ctx.context.eq('liver_spheroid_240h')].strong).sum())
    hep=int(direct[direct.context.eq('iPSC-Hep')].pathway.nunique())
    n_r2=int(truth(cross.R2).sum())
    ox=cross[cross.pathway.eq('Oxidative Phosphorylation')].iloc[0]
    conflict=int(ox.FAMILY_1_vote=='UP' and ox.FAMILY_3_vote=='DOWN')
    # Narrative findings describe the specified analysis; changed inputs must not silently reuse that prose.
    if (cm,liver,hep,n_r2,len(restricted),conflict)!=(0,0,0,4,11,1):
        raise RuntimeError('The narrative findings differ from the manuscript; inspect the numerical results.')
    negative=pd.DataFrame([
      ['N01','iPSC-CM had zero strong primary PFTeDA Hallmark programs under the predefined joint camera/singscore criterion.',cm,'strong programs','Table S3','Does not prove absence of biological effect.'],
      ['N02','The 240 h liver-spheroid context had zero strong primary PFTeDA Hallmark programs.',liver,'strong programs','Table S3','Does not prove absence of biological effect; duration/batch-context differences remain.'],
      ['N03','iPSC-Hep had zero direct-supported PFTeDA-PFOA differential programs.',hep,'direct-supported programs','Table S7','No difference was established under the predefined criterion; equivalence was not tested.'],
      ['N04','Only four of fifty tested Hallmark programs met R2 or R3 reproducibility.',n_r2,'R2/R3 programs of 50 tested','Tables S5 and S6','Non-reproducible programs may remain context-specific or underpowered.'],
      ['N05',f'{len(restricted)} family-program records were model-restricted without formal heterogeneity support.',len(restricted),'family-program records','Table S11','A significant-versus-nonsignificant pattern is not heterogeneity.'],
      ['N06','Oxidative Phosphorylation had a family-level direction conflict: Family 1 voted UP and Family 3 voted DOWN.',conflict,'family-direction-conflict program','Tables S4 and S5','This program has no coherent cross-family reproducible direction and must not be summarized as a shared PFTeDA response.'],
    ],columns=['finding_id','finding','value','unit','source','interpretation_boundary'])
    write(12,negative);write(13,sensitivity_table(run))
    metadata=pd.read_csv(run/'metadata/SAMPLE_METADATA.csv')
    f3=metadata[metadata.dataset.eq('GSE145239')]
    write(14,f3.groupby(['duration_h','treatment','dose_uM','batch'],dropna=False).size().reset_index(name='n'))
    print('Exported Tables S1-S14.')

if __name__=='__main__':main()
