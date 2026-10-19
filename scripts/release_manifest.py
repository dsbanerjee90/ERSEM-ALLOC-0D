#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Write SHA256 manifest of selected scientific deliverables; exclude itself."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
items=[]
for section in ('configs','observations','scripts','experiments','results','figures','reference','candidates'):
    for p in sorted((ROOT/section).rglob('*')):
        if p.is_file() and not p.name.endswith(('.log','.pyc')) and '__pycache__' not in p.parts:
            items.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}')
path=ROOT/'provenance/RELEASE_SHA256SUMS.txt'
path.write_text('\n'.join(items)+'\n')
print('Saved',path,'with',len(items),'files')
