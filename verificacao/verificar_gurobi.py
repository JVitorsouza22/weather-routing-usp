# -*- coding: utf-8 -*-
"""
Verificação independente com o Gurobi dos resultados obtidos por programação dinâmica.
Uso:   python rodar_tudo.py --gurobi     (jeito mais simples; roda tudo e depois esta conferência)
       ou, só esta conferência: entre na pasta verificacao/ e rode  python verificar_gurobi.py
       (precisa dos arquivos modelo_E*.csv e esperado.json, gravados por src/export_gurobi.py)
Saída: tabela no terminal + gurobi_resultado.json (registro do resultado)

O que é verificado, para cada embarcação:
  (1) tempo mínimo e risco-hora mínimo (otimização lexicográfica, garante ponto de Pareto);
  (2) soma ponderada normalizada (ideal-nadir) resolvida como PL (x contínuo) -> confere integralidade;
  (3) método ε-restrição: para cada ponto da fronteira esperada, min F_T s.a. F_H <= ε (PLI, x binário).
"""
import csv, json, time, sys
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
import gurobipy as gp
from gurobipy import GRB

TOL = 1e-6
esperado = json.load(open('esperado.json', encoding='utf-8'))

def ler(k):
    with open(f'modelo_{k}.csv', newline='', encoding='utf-8') as f:
        return [(int(r['i']), int(r['j']), float(r['t_h']), float(r['h_riskh'])) for r in csv.DictReader(f)]

def construir(arcos, o, s, binario):
    m = gp.Model(); m.Params.OutputFlag = 0
    vt = GRB.BINARY if binario else GRB.CONTINUOUS
    x = {(i, j): m.addVar(vtype=vt, lb=0.0, ub=1.0) for (i, j, t, h) in arcos}
    nos = sorted({a[0] for a in arcos} | {a[1] for a in arcos})
    for k in nos:
        saida = gp.quicksum(x[i, j] for (i, j, _, _) in arcos if i == k)
        entra = gp.quicksum(x[i, j] for (i, j, _, _) in arcos if j == k)
        m.addConstr(saida - entra == (1 if k == o else (-1 if k == s else 0)))
    T = gp.quicksum(t * x[i, j] for (i, j, t, h) in arcos)
    H = gp.quicksum(h * x[i, j] for (i, j, t, h) in arcos)
    return m, x, T, H

def valores(x, arcos):
    T = sum(t * x[i, j].X for (i, j, t, h) in arcos); H = sum(h * x[i, j].X for (i, j, t, h) in arcos)
    return T, H

def lexico(arcos, o, s, primeiro):
    """minimiza 'primeiro' ('T' ou 'H') e, entre os ótimos, o outro critério"""
    m, x, T, H = construir(arcos, o, s, True)
    a, b = (T, H) if primeiro == 'T' else (H, T)
    m.setObjective(a, GRB.MINIMIZE); m.optimize(); va = m.ObjVal
    m.addConstr(a <= va + 1e-7); m.setObjective(b, GRB.MINIMIZE); m.optimize()
    return valores(x, arcos), m.NodeCount

res = {'gurobi_version': '.'.join(map(str, gp.gurobi.version())), 'navios': {}}
ok_global = True
for k in ('E1', 'E2', 'E3'):
    e = esperado[k]; arcos = ler(k); o, s = e['origem'], e['destino']
    r = {'n_arcos': len(arcos)}
    print(f'\n===== {k}: {len(arcos)} arcos | Gurobi {res["gurobi_version"]}')
    # (1) extremos lexicográficos
    (Tt, Ht), nn1 = lexico(arcos, o, s, 'T'); (Th, Hh), nn2 = lexico(arcos, o, s, 'H')
    fT, fH = e['fronteira'][0], e['fronteira'][-1]
    d1 = max(abs(Tt - fT['T']), abs(Ht - fT['H'])); d2 = max(abs(Th - fH['T']), abs(Hh - fH['H']))
    print(f'  tempo mínimo   : T={Tt:.4f} H={Ht:.4f}  (esperado T={fT["T"]:.4f} H={fT["H"]:.4f})  dif={d1:.2e}')
    print(f'  risco-hora mín.: T={Th:.4f} H={Hh:.4f}  (esperado T={fH["T"]:.4f} H={fH["H"]:.4f})  dif={d2:.2e}')
    r['extremos_dif_max'] = max(d1, d2); ok = (d1 < TOL and d2 < TOL)
    # (2) soma ponderada normalizada, PL contínuo
    IT, NT, IH, NH = e['ideal']['T'], e['nadir']['T'], e['ideal']['H'], e['nadir']['H']
    r['ponderadas'] = {}
    for nome, w in e['ponderadas'].items():
        wT, wH = w['w']
        m, x, T, H = construir(arcos, o, s, False)
        m.setObjective(wT * (T - IT) / (NT - IT) + wH * (H - IH) / (NH - IH), GRB.MINIMIZE); m.optimize()
        Tv, Hv = valores(x, arcos); frac = max(min(v.X, 1 - v.X) for v in x.values())
        dif = max(abs(Tv - w['T']), abs(Hv - w['H']))
        print(f'  soma ponderada {nome:11s} w={tuple(w["w"])}: T={Tv:.4f} H={Hv:.4f} dif={dif:.2e}  maior fração={frac:.1e} (integral se ~0)')
        r['ponderadas'][nome] = dict(dif=dif, fracao_max=frac); ok &= (dif < TOL and frac < 1e-7)
    # (3) ε-restrição sobre cada ponto da fronteira esperada
    difs = []; nos_bb = 0; t0 = time.time()
    for f in e['fronteira']:
        m, x, T, H = construir(arcos, o, s, True)
        m.addConstr(H <= f['H'] + 1e-7); m.setObjective(T, GRB.MINIMIZE); m.optimize()
        difs.append(abs(m.ObjVal - f['T'])); nos_bb += int(m.NodeCount)
    r['epsilon'] = dict(n_pontos=len(e['fronteira']), dif_max=max(difs), nos_bb_total=nos_bb, tempo_s=time.time() - t0)
    print(f'  ε-restrição: {len(difs)} pontos da fronteira; dif máx={max(difs):.2e}; nós de B&B (soma)={nos_bb}; tempo={time.time()-t0:.1f}s')
    ok &= max(difs) < TOL
    r['ok'] = bool(ok); ok_global &= ok; res['navios'][k] = r
    print('  =>', 'CONFERE com a programação dinâmica' if ok else 'DIVERGENCIA - avisar')
res['ok_global'] = bool(ok_global)
json.dump(res, open('gurobi_resultado.json', 'w', encoding='utf-8'), indent=1)
print('\nResultado global:', 'TUDO CONFERE' if ok_global else 'HÁ DIVERGÊNCIA')
print('Arquivo gerado: verificacao/gurobi_resultado.json')
