"""Public arithmetic boundaries; no patient data or network."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts/audit_nucleoside_exposure.py'
SPEC=importlib.util.spec_from_file_location('nucleoside_exposure',SCRIPT)
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_cli_distinguishes_peak_eight_hour_average_and_unknown_daily_exposure(tmp_path):
    output=tmp_path/'audit.json'
    subprocess.run([sys.executable,str(SCRIPT),'--output',str(output)],check=True,capture_output=True)
    result=json.loads(output.read_text())
    day28,day160=result['profiles']
    assert day28['peak_to_nominal_comparator_ratio']==pytest.approx(88.7/30)
    assert day28['mean_over_first_8h_uM']==pytest.approx(278.7/8)
    assert day160['mean_over_first_8h_uM']==pytest.approx(465.6/8)
    assert all(r['mean_over_24h_uM'] is None and r['measured_trough_24h_uM'] is None for r in result['profiles'])
    assert result['therapeutic_index'] is None
    assert result['human_cancer_safe_threshold_uM'] is None
    assert result['patient_specific_exposure'] is None
    assert not result['comparator_is_verified_uridine_alone_effective_concentration']


def test_post_peak_sensitivity_retains_endogenous_baseline():
    time=AUDIT.illustrative_post_peak_time(85,2,5,25)
    assert time==pytest.approx(4)
    with pytest.raises(ValueError):
        AUDIT.illustrative_post_peak_time(85,2,5,3)


def test_parent_drug_or_wrong_units_cannot_masquerade_as_uridine():
    inputs=json.loads((ROOT/'results/nucleoside_inputs.json').read_text())
    inputs['pk_analyte']='uridine triacetate'
    with pytest.raises(ValueError,match='measured uridine'):
        AUDIT.calculate(inputs)
    inputs['pk_analyte']='uridine'
    inputs['concentration_unit']='ng/mL'
    with pytest.raises(ValueError,match='micromolar'):
        AUDIT.calculate(inputs)
