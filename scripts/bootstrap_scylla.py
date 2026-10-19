#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Copy a frozen FABM YAML into this stand-alone PAPER repository (no source edits)."""
from pathlib import Path
import hashlib, json, os, platform, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get('ERSEM_ALLOC_ROOT', '/work/dba/ERSEM-ALLOC-0D'))
ORIGINAL = SOURCE/'experiments/05_DYNAMIC_DT300/fabm.yaml'
DEST = ROOT/'configs/fabm.yaml'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(2**20),b''):
            h.update(part)
    return h.hexdigest()

def git_commit(path):
    if not path.exists():
        return None
    p=subprocess.run(['git','-C',str(path),'rev-parse','HEAD'],capture_output=True,text=True)
    return p.stdout.strip() if p.returncode == 0 else None

def main():
    if not ORIGINAL.is_file():
        raise SystemExit('Missing source YAML (no copy made): '+str(ORIGINAL))
    if DEST.exists():
        if sha(DEST)!=sha(ORIGINAL):
            raise SystemExit('STOP: configs/fabm.yaml exists but differs from source YAML. No overwrite. Inspect manually.')
        print('Configuration already staged; SHA256 matches original.')
    else:
        shutil.copy2(ORIGINAL,DEST)
        print('Copied original YAML to:',DEST)
    b4_yaml=ROOT/'experiments/B4_dynamic_vs_fixed/inputs/fabm.yaml'
    b4_yaml.parent.mkdir(parents=True, exist_ok=True)
    if b4_yaml.exists():
        if sha(b4_yaml)!=sha(DEST):
            raise SystemExit('STOP: staged B4 configuration differs from frozen YAML')
    else:
        shutil.copy2(DEST,b4_yaml)
        print('Copied frozen YAML to B4 experiment:',b4_yaml)
    srcfile=SOURCE/'src/ersem/src/primary_producer.F90'
    pyfabm_path=SOURCE/'paper/pyfabm_tests/python/pyfabm/__init__.py'
    audit={
        'date_recorded_utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
        'python_version': sys.version,
        'platform': platform.platform(),
        'external_ersem_alloc_project_root':str(SOURCE),
        'frozen_yaml_source':str(ORIGINAL),
        'frozen_yaml_sha256':sha(DEST),
        'model_fortran_path':str(srcfile) if srcfile.exists() else None,
        'model_fortran_sha256':sha(srcfile) if srcfile.exists() else None,
        'model_source_git_commit':git_commit(SOURCE/'src/ersem'),
        'fabm_source_git_commit':git_commit(SOURCE/'src/fabm'),
        'external_pyfabm_module_path':str(pyfabm_path) if pyfabm_path.exists() else None,
        'external_pyfabm_module_sha256':sha(pyfabm_path) if pyfabm_path.exists() else None,
        'redistribution_note':'Source, build flags, compiled binary and licence must be verified before a public release.',
    }
    report=ROOT/'provenance/scylla_environment.json'
    report.write_text(json.dumps(audit,indent=2)+'\n')
    print('Wrote:',report)

if __name__=='__main__':main()
