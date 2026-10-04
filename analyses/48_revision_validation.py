"""Recompute the idealized proxy experiment after fixing the valley score (numbers only; Figure S1 is drawn by make_figures.fig4_validation).

This evaluates the four specified profile generators, not all possible mimics
and not the deterministic IDSS estimator used in the revised recovery analysis.
"""
import importlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
v = importlib.import_module("04_validation_report")


def main():
    n=500; hits={k:0 for k in v.GENERATORS}
    for seed in range(n):
        result=v.discriminates(v.run_blind(v.GENERATORS,seed),deriv_threshold=v.DERIV_THRESHOLD)
        for key,r in result.items(): hits[key]+=int(r['convergent'])
    summary={key:dict(hits=h,total=n,interval=[float(beta.ppf(.025,h,n-h+1)) if h else 0.,
              float(beta.ppf(.975,h+1,n-h)) if h<n else 1.]) for key,h in hits.items()}
    (ROOT/'output/revision_idealized_validation.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
