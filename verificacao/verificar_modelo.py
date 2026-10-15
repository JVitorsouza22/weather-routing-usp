# -*- coding: utf-8 -*-
"""Conferências extras do modelo (não precisa do Gurobi). Rode depois de  python rodar_tudo.py

  1. Varredura de pesos. A soma ponderada, com 4.000 pesos entre 0 e 1, deve achar exatamente os pontos
     "suportados" da fronteira de Pareto (os que ficam na envoltória convexa inferior) e nenhum outro.
  2. Piso de velocidade. O código impede que a velocidade sobre o fundo fique abaixo de 0,4 V (proteção numérica).
     Confere que esse piso nunca é usado nestes dados.
  3. Números citados no texto do artigo (tamanho dos arcos, corrente média ao norte de 8°S, deslocamentos laterais).

Uso:   python verificacao/verificar_modelo.py
"""
import pathlib, sys
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / 'src'))
from load_ships import *          # noqa: F401,F403  (load, OUTDIR, lower_hull_mask, ...)

NAVIOS = ('E1', 'E2', 'E3')
ok_geral = True

print('=== 1) Varredura de pesos: a soma ponderada acha só os pontos suportados ===')
for k in NAVIOS:
    S = load(k)
    fr = S.pareto(); suportado = lower_hull_mask(fr)
    I_T, N_T = fr[0]['T'], fr[-1]['T']; I_H, N_H = fr[-1]['H'], fr[0]['H']
    achados = set()
    for wT in np.concatenate([np.linspace(0.0005, 0.9995, 4000), [0.5, 0.8, 0.2]]):
        caminho = S.scalar_dp((wT) / (N_T - I_T), (1 - wT) / (N_H - I_H))
        m = S.path_metrics(caminho)
        achados.add((round(m['time'], 9), round(m['rh'], 9)))
    sup = {(round(f['T'], 9), round(f['H'], 9)) for f, s in zip(fr, suportado) if s}
    todos = {(round(f['T'], 9), round(f['H'], 9)) for f in fr}
    ok = (achados == sup) and (achados <= todos)
    ok_geral &= ok
    print(f'{k}: pontos achados pela soma ponderada = {len(achados)}; suportados na envoltória = {len(sup)}; '
          f'de {len(todos)} eficientes  ->', 'OK' if ok else 'FALHOU')

print('\n=== 2) Piso de 0,4 V nunca é usado ===')
for k in NAVIOS:
    for com_corrente in (True, False):
        S = load(k, use_current=com_corrente)
        v_fundo = np.array([S.V + r.u for r in S.A.itertuples()])
        usado = bool((v_fundo < 0.4 * S.V).any())
        ok_geral &= not usado
        print(f'{k} {"com corrente" if com_corrente else "sem corrente"}: V = {S.V:.3f} nós; menor velocidade sobre o fundo em um arco = '
              f'{v_fundo.min():.3f} nós; piso (0,4 V) = {0.4 * S.V:.3f} nós  ->', 'piso não usado, OK' if not usado else 'PISO USADO')

print('\n=== 3) Números citados no texto ===')
cordas = []
for k in NAVIOS:
    S = load(k)
    ds = [S.chord(i, j)[0] for i, j in zip(S.anchors[:-1], S.anchors[1:])]
    print(f'{k}: cordas entre pontos de referência consecutivos (nm) = {[round(x) for x in ds]}; soma = {sum(ds):.1f}')
    cordas += ds
print(f'Todas as cordas: mediana {np.median(cordas):.1f} nm, mínimo {min(cordas):.1f}, máximo {max(cordas):.1f} (n = {len(cordas)})')
for k in NAVIOS:
    S = load(k)
    print(f'{k}: maior arco do grafo = {S.A.d.max():.1f} nm; arcos no grafo = {len(S.A)}')
cnt = ce = cn = 0
for k in NAVIOS:
    n = load(k).nodes
    q = n[n.queried & (n.Latitude > -8.0)]
    cnt += len(q); ce += q.cu.sum(); cn += q.cv.sum()
print(f'Nós consultados ao norte de 8°S: {cnt}; corrente média (km/h): leste = {ce / cnt:.2f}, norte = {cn / cnt:.2f}')
for k in NAVIOS:
    S = load(k); n = S.nodes; fr = S.pareto()
    for nome, f in (('tempo mínimo', fr[0]), ('risco-hora mínimo', fr[-1])):
        desl = [float(n.at[i, 'offset_nm']) for i in f['path'] if n.at[i, 'queried']]
        interior = desl[1:-1]
        print(f'{k} {nome}: colunas interiores = {len(interior)}; com |deslocamento| = 30 nm: {sum(abs(x) == 30 for x in interior)}; deslocamentos = {desl}')

print('\nRESULTADO:', 'tudo OK' if ok_geral else 'HÁ FALHAS (veja as linhas marcadas)')
