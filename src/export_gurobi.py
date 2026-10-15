# -*- coding: utf-8 -*-
"""Exporta o modelo em rede (arcos com t e h) e os resultados esperados (rotulagem) para verificacao/."""
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from load_ships import *
RES = json.load(open(OUTDIR / 'results.json', encoding='utf-8'))
V = ROOT / 'verificacao'
esp = {}
for k in ('E1', 'E2', 'E3'):
    S = load(k)
    A = S.A[['i', 'j', 'd', 't', 'h']].rename(columns={'d': 'd_nm', 't': 't_h', 'h': 'h_riskh'})
    A.to_csv(V / f'modelo_{k}.csv', index=False, float_format='%.10f')
    r = RES[k]
    esp[k] = dict(origem=int(S.anchors[0]), destino=int(S.anchors[-1]), n_arcos=int(len(A)),
                  fronteira=[dict(T=f['T'], H=f['H'], suportado=bool(f['supported'])) for f in r['frontier']],
                  ideal=dict(T=r['payoff']['ideal']['T'], H=r['payoff']['ideal']['H']),
                  nadir=dict(T=r['payoff']['nadir']['T'], H=r['payoff']['nadir']['H']),
                  ponderadas={n: dict(w=list(v['w']), T=v['T'], H=v['H']) for n, v in r['weighted'].items()},
                  n_caminhos=int(r['enumeration']['n_paths']))
json.dump(esp, open(V / 'esperado.json', 'w', encoding='utf-8'), indent=1)
print('exportado para', V)
