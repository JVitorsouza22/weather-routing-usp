# -*- coding: utf-8 -*-
import sys, json, pathlib, numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from load_ships import *
from geo_check import path_samples, route_geo_report
from global_land_mask import globe

G = get_grid()
MESH = json.load(open(OUTDIR / 'mesh_stats.json', encoding='utf-8'))
TRANSITS = {'E1': 18, 'E2': 18, 'E3': 14}          # trânsitos/ano estimados (Seção 5.8)
WEIGHTS = {'tempo': (0.8, 0.2), 'equilibrado': (0.5, 0.5), 'risco': (0.2, 0.8)}


def fl(x):
    return float(x)


def metrics(S, path):
    m = S.path_metrics(path)
    return {k: (fl(v) if not isinstance(v, int) else v) for k, v in m.items()}


def pareto_from_arrays(T, H):
    o = np.lexsort((H, T)); T, H = T[o], H[o]
    keep, best = [], np.inf
    for t, h in zip(T, H):
        if h < best - 1e-9:
            keep.append((t, h)); best = h
    return keep


def analyse(k):
    S = load(k)
    n = S.nodes
    R = {}
    # ---------------- dados e velocidade
    sog = SOG[k]
    dens, dur_all, _ = S.real_speed_stats()
    R['data'] = dict(v_obs=sog['v_obs'], v_med=sog['v_med'], n_nav=sog['n_nav'], n_samples=sog['n'], n_invalid=sog['n_invalid'],
                     lag_med=sog['lag_med'], lag_med_abs=sog['lag_med_abs'], lag_p90=sog['lag_p90'], lag_max=sog['lag_max'], pct_stale=sog['pct_stale'],
                     dense_dist=fl(dens), window_h=sog['dur_h'],
                     u_bar_real=fl(S.u_bar_real), V_water=fl(S.V))
    R['mesh'] = dict(MESH[k])
    R['mesh'].update(nodes_total=int(len(n)), n_queried=int(n.queried.sum()), n_wp=int((~n.queried).sum()), arcs=int(len(S.A)),
                     paths=int(S.count_paths()))
    # ---------------- contrafactual
    chords = metrics(S, S.real_path())
    proj_path = S.real_projected(); proj = metrics(S, proj_path)
    R['real'] = dict(chords=chords, proj=proj, dense_dist=fl(dens),
                     sinuosity=fl(dens / chords['dist']), proj_offsets=[fl(x) for x in S.real_proj_offsets],
                     proj_detour_pct=fl(100 * (proj['dist'] / chords['dist'] - 1)))
    # ---------------- fronteira exata
    fr = S.pareto(); hm = lower_hull_mask(fr)
    F = []
    for f, sup in zip(fr, hm):
        m = metrics(S, f['path'])
        F.append(dict(T=fl(f['T']), H=fl(f['H']), dist=m['dist'], rbar=m['rbar'], supported=bool(sup), path=[int(i) for i in f['path']]))
    R['frontier'] = F
    R['frontier_summary'] = dict(n=len(F), n_supported=int(sum(hm)), n_unsupported=int(len(F) - sum(hm)))
    # verificação por enumeração exaustiva
    Tall, Hall = S.enumerate_all(6_000_000)
    ef = pareto_from_arrays(Tall, Hall)
    same = (len(ef) == len(F)) and all(abs(a[0] - f['T']) < 1e-7 and abs(a[1] - f['H']) < 1e-7 for a, f in zip(ef, F))
    R['enumeration'] = dict(n_paths=int(len(Tall)), frontier_equal=bool(same), n_frontier_enum=len(ef),
                            T_range=[fl(Tall.min()), fl(Tall.max())], H_range=[fl(Hall.min()), fl(Hall.max())])
    np.savez_compressed(OUTDIR / f'cloud_{k}.npz', T=Tall.astype(np.float32), H=Hall.astype(np.float32))
    # ---------------- payoff, ideal e nadir
    a, b = F[0], F[-1]
    I_T, N_T, I_H, N_H = a['T'], b['T'], b['H'], a['H']
    R['payoff'] = dict(minT=dict(T=a['T'], H=a['H'], dist=a['dist'], rbar=a['rbar']),
                       minH=dict(T=b['T'], H=b['H'], dist=b['dist'], rbar=b['rbar']),
                       ideal=dict(T=I_T, H=I_H), nadir=dict(T=N_T, H=N_H),
                       range_T=N_T - I_T, range_H=N_H - I_H)
    # ---------------- soma ponderada normalizada
    ws = {}
    for name, (wT, wH) in WEIGHTS.items():
        if N_T - I_T > 1e-9 and N_H - I_H > 1e-9:
            path = S.scalar_dp(wT / (N_T - I_T), wH / (N_H - I_H))
        else:
            path = F[0]['path']
        m = metrics(S, path)
        # confere por enumeração da fronteira (mesmo valor escalar)
        vals = [wT * (f['T'] - I_T) / max(N_T - I_T, 1e-12) + wH * (f['H'] - I_H) / max(N_H - I_H, 1e-12) for f in F]
        zdp = wT * (S.path_metrics(path)['time'] - I_T) / max(N_T - I_T, 1e-12) + wH * (S.path_metrics(path)['rh'] - I_H) / max(N_H - I_H, 1e-12)
        ws[name] = dict(w=(wT, wH), T=m['time'], H=m['rh'], dist=m['dist'], rbar=m['rbar'], Z=fl(zdp), Z_frontier_min=fl(min(vals)),
                        match=bool(abs(zdp - min(vals)) < 1e-9), path=[int(i) for i in path])
    R['weighted'] = ws
    # ---------------- ε-restrição (ganhos contrafactuais)
    Tp, Hp = proj['time'], proj['rh']
    fe1 = [f for f in F if f['H'] <= Hp + 1e-9]
    e1 = fe1[0] if fe1 else None                       # min T s.a. H <= H_real
    fe2 = [f for f in F if f['T'] <= Tp + 1e-9]
    e2 = fe2[-1] if fe2 else None                      # min H s.a. T <= T_real
    R['epsilon'] = dict(minT_given_H=e1, minH_given_T=e2,
                        dominated_by=int(sum(1 for f in F if f['T'] <= Tp + 1e-9 and f['H'] <= Hp + 1e-9)))
    # ---------------- distância "ingênua" (rastro denso vs otimizada em cordas)
    R['naive'] = dict(dist_min=min(f['dist'] for f in F), dist_max=max(f['dist'] for f in F),
                      gain_vs_dense_pct=[fl(100 * (f['dist'] / dens - 1)) for f in (F[0], F[-1])],
                      gain_vs_chords_pct=[fl(100 * (f['dist'] / chords['dist'] - 1)) for f in (F[0], F[-1])],
                      gain_vs_proj_pct=[fl(100 * (f['dist'] / proj['dist'] - 1)) for f in (F[0], F[-1])])
    # ---------------- geografia das rotas de referência
    refs = {'real_proj': proj_path, 'minT': F[0]['path'], 'minH': F[-1]['path'], 'equilibrado': ws['equilibrado']['path']}
    geo = {}
    for name, p in refs.items():
        g = route_geo_report(n, p)
        # verificação de terra com a máscara exata (1 nm) fora dos 3 nm do pixel portuário
        geo[name] = g
    R['geo'] = geo
    # ---------------- sensibilidade
    sens = {}
    for label, (wp, cur) in {'principal': ('idw', True), 'linear': ('linear', True), 'sem_corrente': ('idw', False)}.items():
        Sx = load(k, wp=wp, use_current=cur); frx = Sx.pareto(); px = Sx.path_metrics(Sx.real_projected())
        sens[label] = dict(V=fl(Sx.V), real_T=px['time'], real_H=px['rh'], n=len(frx),
                           minT=dict(dT=fl(100 * (frx[0]['T'] / px['time'] - 1)), dH=fl(100 * (frx[0]['H'] / px['rh'] - 1))),
                           minH=dict(dT=fl(100 * (frx[-1]['T'] / px['time'] - 1)), dH=fl(100 * (frx[-1]['H'] / px['rh'] - 1))))
    R['sens'] = sens
    # ---------------- risco médio por arco (métrica antiga) das rotas de referência
    R['r_arc_mean'] = dict(real_proj=proj['r_arc_mean'], minT=metrics(S, F[0]['path'])['r_arc_mean'], minH=metrics(S, F[-1]['path'])['r_arc_mean'],
                           n_arcs_real=proj['n_arcs'], n_arcs_minT=metrics(S, F[0]['path'])['n_arcs'])
    # ---------------- projeção (ilustrativa)
    tr = TRANSITS[k]
    dT_minT = proj['time'] - F[0]['T']; dT_eq = proj['time'] - ws['equilibrado']['T']
    R['projection'] = dict(transits_year=tr, dT_trip_minT=fl(dT_minT), dT_trip_eq=fl(dT_eq),
                           per_period={str(m): dict(trips=tr * m / 12, hours_minT=fl(dT_minT * tr * m / 12), hours_eq=fl(dT_eq * tr * m / 12)) for m in (3, 6, 12)})
    R['V'] = fl(S.V)
    return R


if __name__ == '__main__':
    RES = {k: analyse(k) for k in ('E1', 'E2', 'E3')}
    json.dump(RES, open(OUTDIR / 'results.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    for k, r in RES.items():
        print(k, 'fronteira', r['frontier_summary'], '| enum igual:', r['enumeration']['frontier_equal'], '| dominada por', r['epsilon']['dominated_by'])
        print('   ponderadas ok:', {n: v['match'] for n, v in r['weighted'].items()})
        print('   geo minT:', {a: round(b, 2) for a, b in r['geo']['minT'].items()})
