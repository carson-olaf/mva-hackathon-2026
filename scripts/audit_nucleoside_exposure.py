"""Reproduce bounded uridine exposure arithmetic from curated public facts.

This is not a dose calculator. It never accepts patient weight, genotype or PK.
Unknown efficacy and safety thresholds stay unknown; a peak ratio cannot fill
them. No network or private files are used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def illustrative_post_peak_time(cmax, half_life, baseline, threshold):
    """Return a conditional post-peak crossing time, not measured coverage."""
    if half_life <= 0 or not 0 <= baseline < threshold < cmax:
        raise ValueError('Require positive half-life and baseline < threshold < peak')
    return half_life * math.log2((cmax-baseline)/(threshold-baseline))


def calculate(inputs):
    if inputs['pk_analyte']!='uridine' or inputs['concentration_unit']!='micromolar':
        raise ValueError('Comparison requires measured uridine in micromolar units')
    comparator=inputs['comparison']['nominal_comparator_uM']
    if comparator <= 0:
        raise ValueError('Comparator must be positive')
    rows=[]
    for profile in inputs['profiles']:
        rows.append({
            'profile':profile['id'], 'n':profile['n'],
            'cmax_mean_uM':profile['cmax_mean_uM'],
            'peak_to_nominal_comparator_ratio':profile['cmax_mean_uM']/comparator,
            'measured_auc_window_h':8,
            'mean_over_first_8h_uM':profile['auc_0_8_mean_uM_h']/8,
            'mean_over_24h_uM':None,
            'measured_trough_24h_uM':None,
            'measured_hours_above_comparator_per_day':None,
            'illustrative_post_peak_crossing_times':[
                {'assumed_baseline_uM':baseline,
                 'hours_after_peak':illustrative_post_peak_time(
                     profile['cmax_mean_uM'],profile['half_life_mean_h'],baseline,comparator),
                 'status':'Sensitivity calculation only; not measured or patient-specific'}
                for baseline in inputs['comparison']['illustrative_baselines_uM']
            ],
        })
    return {
        'nominal_comparator_uM':comparator,
        'comparator_is_verified_uridine_alone_effective_concentration':False,
        'profiles':rows,
        'uridine_alone_bub1b_efficacy_verified':False,
        'cancer_hazard_reported_concentrations_uM':inputs['cancer_hazard']['reported_uridine_cell_exposures_uM'],
        'human_cancer_safe_threshold_uM':None,
        'therapeutic_index':None,
        'patient_specific_exposure':None,
        'limits':[
            'Label concentrations are plasma uridine, not parent uridine triacetate.',
            'AUC0-8 divided by8 is an eight-hour average, not a daily average.',
            'Means from different small groups do not establish dose proportionality.',
            'SD is not a confidence interval or an individual patient range.',
            'Conditional post-peak decay omits the rising limb, absorption, distribution and intracellular retention; it is not total daily time above the comparator.',
            'The crossing-time range varies assumed baseline only, not peak or half-life variability; it is not an uncertainty interval.',
            '30 micromolar is a reference concentration, not an established uridine efficacy threshold; the peak ratio is not an efficacy margin, adequate target coverage or a therapeutic window.',
            'A plasma peak does not establish the required duration or all four mixture components.',
            'Cellular DNA-damage reports are not a human toxicity threshold or clinical cancer risk ratio.'
        ]
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'results/nucleoside_inputs.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/nucleoside_exposure_audit.json')
    args=parser.parse_args()
    raw=args.input.read_bytes()
    result=calculate(json.loads(raw))
    result['input_sha256']=hashlib.sha256(raw).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
