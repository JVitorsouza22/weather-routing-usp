# -*- coding: utf-8 -*-
"""Caminhos e configuração centralizados (nada de caminhos absolutos)."""
import json, pathlib, pickle
import sys
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
OUTDIR = ROOT / 'out'
MESHDIR = OUTDIR / 'mesh'
FIGDIR = OUTDIR / 'fig'
for _d in (OUTDIR, MESHDIR, FIGDIR):
    _d.mkdir(parents=True, exist_ok=True)
CFG = json.load(open(ROOT / 'config.json', encoding='utf-8'))


def get_grid():
    """Grade de folga da costa (nm) a partir da máscara GLOBE 1 km; construída uma vez e guardada em out/."""
    f = OUTDIR / 'clearance_grid.pkl'
    if f.exists():
        return pickle.load(open(f, 'rb'))
    from clearance import ClearGrid
    G = ClearGrid(res=CFG['RES_GRAU'])
    pickle.dump(G, open(f, 'wb'))
    return G
