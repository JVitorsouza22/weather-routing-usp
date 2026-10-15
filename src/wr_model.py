# -*- coding: utf-8 -*-
"""
Modelo bicritério de weather routing em rede acíclica (nós por coluna).
Critérios efetivos: tempo F_T (proporcional ao consumo a potência constante) e
risco-hora F_H = sum(r_ij * t_ij).  Tudo aditivo sobre os arcos.
"""
import numpy as np, pandas as pd
from collections import defaultdict

R_NM = 3440.065        # raio da Terra em milhas náuticas
KMH_TO_KN = 1.0 / 1.852

# ---------------------------------------------------------------- geometria
def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dl = np.radians(lon2) - np.radians(lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_NM * np.arcsin(np.sqrt(a))

def bearing(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2) - np.radians(lon1)
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.arctan2(y, x)            # radianos, 0 = norte, horário

# ---------------------------------------------------------------- índice de risco
def risk_index(hs, wind_kmh, cur_kmh):
    return (0.41 * np.minimum(hs / 4.0, 1.0)
            + 0.41 * np.minimum(wind_kmh / 62.0, 1.0)
            + 0.18 * np.minimum(cur_kmh / 7.5, 1.0))

# ---------------------------------------------------------------- carga
class Ship:
    def __init__(self, key, params, nodes_path, arcs_path,
                 wp_treatment='idw', K_idw=4):
        self.key = key
        nodes = pd.read_csv(nodes_path)
        arcs = pd.read_csv(arcs_path)
        self.params = params          # escalares derivados do rastro AIS (o rastro bruto não é distribuído)
        nodes = nodes.set_index('id').sort_index()
        nodes['queried'] = (nodes['coluna'] == nodes['coluna'].round())
        self.nodes = nodes
        self.arcs = arcs[['origem', 'destino']].copy()
        self.wp_treatment = wp_treatment
        self._complete_weather(K_idw)
        self._build_topology()

    # -- clima dos waypoints (sem consulta à API) ------------------------
    def _complete_weather(self, K):
        n = self.nodes
        n['cu'] = np.nan; n['cv'] = np.nan            # corrente "para onde flui" (km/h)
        th = np.radians(n['corrente_dir_graus'])
        n.loc[n.queried, 'cu'] = (n['corrente_vel_kmh'] * np.sin(th))[n.queried]
        n.loc[n.queried, 'cv'] = (n['corrente_vel_kmh'] * np.cos(th))[n.queried]
        wp = n[~n.queried]
        if len(wp) == 0:
            return
        if self.wp_treatment == 'linear':
            # tratamento "ingênuo": campos interpolados linearmente ao longo da cadeia entre os dois nós
            # consultados que ela liga (inclui, portanto, nós portuários abrigados)
            for i, r in wp.iterrows():
                a, b, f = int(r.chain_from), int(r.chain_to), float(r.chain_frac)
                for fld in ('vento_vel_kmh', 'onda_altura_m', 'cu', 'cv'):
                    n.at[i, fld] = (1 - f) * n.at[a, fld] + f * n.at[b, fld]
                n.at[i, 'corrente_vel_kmh'] = float(np.hypot(n.at[i, 'cu'], n.at[i, 'cv']))
            n['R_use'] = n['R']
            rr = risk_index(n.loc[~n.queried, 'onda_altura_m'], n.loc[~n.queried, 'vento_vel_kmh'],
                            n.loc[~n.queried, 'corrente_vel_kmh'])
            n.loc[~n.queried, 'R_use'] = rr
            return
        # 'idw': inverso do quadrado da distância, K doadores mais próximos entre nós consultados,
        # excluindo os pontos de referência de origem e destino (águas portuárias abrigadas)
        anc = n[(n.offset_nm == 0) & (n.coluna.isin([n.coluna.min(), 9.0]))].index
        donors = n[n.queried & ~n.index.isin(anc)]
        for i, r in wp.iterrows():
            d = haversine(r.Latitude, r.Longitude, donors.Latitude.values, donors.Longitude.values)
            idx = np.argsort(d)[:K]
            w = 1.0 / np.maximum(d[idx], 1e-6) ** 2
            w /= w.sum()
            dd = donors.iloc[idx]
            n.at[i, 'vento_vel_kmh'] = float((w * dd.vento_vel_kmh.values).sum())
            n.at[i, 'onda_altura_m'] = float((w * dd.onda_altura_m.values).sum())
            n.at[i, 'cu'] = float((w * dd.cu.values).sum())
            n.at[i, 'cv'] = float((w * dd.cv.values).sum())
            n.at[i, 'corrente_vel_kmh'] = float(np.hypot(n.at[i, 'cu'], n.at[i, 'cv']))
        n['R_use'] = n['R']
        rr = risk_index(n.loc[~n.queried, 'onda_altura_m'], n.loc[~n.queried, 'vento_vel_kmh'],
                        n.loc[~n.queried, 'corrente_vel_kmh'])
        n.loc[~n.queried, 'R_use'] = rr

    # -- topologia --------------------------------------------------------
    def _build_topology(self):
        n = self.nodes
        self.o = int(n[(n.coluna == n.coluna.min()) & (n.offset_nm == 0)].index[0])
        self.s = int(n[(n.coluna == 9.0) & (n.offset_nm == 0)].index[0])
        col = n['coluna'].to_dict()
        a = self.arcs
        assert all(col[d] > col[o] for o, d in zip(a.origem, a.destino)), 'arco nao avanca coluna'
        self.order = list(n.sort_values(['coluna']).index)
        # pontos de referência da rota real (chamados de 'anchors' no código) = nós de desvio 0 nas colunas inteiras
        anc = n[(n.offset_nm == 0) & n.queried].sort_values('coluna')
        self.anchors = [int(i) for i in anc.index]

    # -- parâmetros por arco ---------------------------------------------
    def chord(self, i, j):
        """distância, proa, corrente média (kn: leste,norte) e risco médio de um arco (i,j)"""
        n = self.nodes
        li, oi, lj, oj = n.at[i, 'Latitude'], n.at[i, 'Longitude'], n.at[j, 'Latitude'], n.at[j, 'Longitude']
        d = float(haversine(li, oi, lj, oj))
        b = float(bearing(li, oi, lj, oj))
        cu = 0.5 * (n.at[i, 'cu'] + n.at[j, 'cu']) * KMH_TO_KN
        cv = 0.5 * (n.at[i, 'cv'] + n.at[j, 'cv']) * KMH_TO_KN
        u_along = cu * np.sin(b) + cv * np.cos(b)
        r = 0.5 * (n.at[i, 'R_use'] + n.at[j, 'R_use'])
        return d, b, u_along, r

    def real_speed_stats(self):
        """velocidade média sobre o fundo observada (rastro AIS denso) e duração"""
        dens = float(self.params['dense_dist_nm']); dur = float(self.params['dur_h'])
        return dens, dur, dens / dur

    def calibrate(self, v_obs, use_current=True):
        """V_w = V_sog_obs (média do SOG reportado em navegação) - componente média da
        corrente ao longo da rota real (nas cordas entre pontos de referência, ponderada pela distância)"""
        ds, us = [], []
        for i, j in zip(self.anchors[:-1], self.anchors[1:]):
            d, b, u, r = self.chord(i, j)
            ds.append(d); us.append(u)
        ds, us = np.array(ds), np.array(us)
        ubar = float((ds * us).sum() / ds.sum()) if use_current else 0.0
        self.v_sog = v_obs
        self.u_bar_real = ubar
        self.V = v_obs - ubar
        self.use_current = use_current
        return self.V

    def arc_params(self, i, j):
        d, b, u, r = self.chord(i, j)
        u_eff = u if self.use_current else 0.0
        v_ground = max(self.V + u_eff, 0.4 * self.V)          # piso numérico (evita dividir por ~0). Nunca é usado nestes dados:
        # a menor velocidade sobre o fundo em qualquer arco é 9,0 nós (verificacao/verificar_modelo.py), contra um piso de 4,1 nós
        t = d / v_ground
        return dict(d=d, t=t, r=r, h=r * t, u=u_eff)

    def build_arcs(self):
        rows = []
        for i, j in zip(self.arcs.origem, self.arcs.destino):
            p = self.arc_params(int(i), int(j)); p['i'] = int(i); p['j'] = int(j)
            rows.append(p)
        self.A = pd.DataFrame(rows)
        self.out = defaultdict(list)
        for k, r in enumerate(self.A.itertuples()):
            self.out[r.i].append((r.j, r.t, r.h, r.d, r.r, k))
        return self.A

    # -- métricas de um caminho -------------------------------------------
    def path_metrics(self, path):
        d = t = h = 0.0; rs = []
        for i, j in zip(path[:-1], path[1:]):
            p = self.arc_params(i, j)
            d += p['d']; t += p['t']; h += p['h']; rs.append(p['r'])
        return dict(dist=d, time=t, rh=h, rbar=h / t, r_arc_mean=float(np.mean(rs)), n_arcs=len(rs))


    # -- rota real projetada na malha ------------------------------------
    def shortest_dist_path(self, a, b):
        """caminho de menor distância de a até b dentro de A (DAG); None se inalcançável"""
        n = self.nodes
        best = {a: 0.0}; pred = {}
        ca, cb = n.at[a, 'coluna'], n.at[b, 'coluna']
        for u in self.order:
            if u not in best or n.at[u, 'coluna'] >= cb:
                continue
            for (j, t, h, d, r, k) in self.out.get(u, []):
                if n.at[j, 'coluna'] > cb:
                    continue
                v = best[u] + d
                if v < best.get(j, np.inf) - 1e-9:
                    best[j] = v; pred[j] = u
        if b not in best:
            return None
        path, node = [b], b
        while node != a:
            node = pred[node]; path.append(node)
        return path[::-1]

    def real_projected(self):
        """Projeção da rota real em A: caminho que minimiza o desvio lateral total em relação aos
        pontos de referência (soma de |desvio| dos nós visitados; desempate por distância). Reproduz a própria
        sequência de pontos de referência quando todos os arcos de desvio 0 pertencem a A."""
        n = self.nodes
        off = n['offset_nm'].abs().to_dict()
        best = {self.o: 0.0}; pred = {}
        for u in self.order:
            if u not in best:
                continue
            for (j, t, h, d, r, k) in self.out.get(u, []):
                v = best[u] + 1000.0 * off[j] + d
                if v < best.get(j, np.inf) - 1e-9:
                    best[j] = v; pred[j] = u
        path, node = [self.s], self.s
        while node != self.o:
            node = pred[node]; path.append(node)
        path = path[::-1]
        self.real_proj_fallback = 0
        self.real_proj_offsets = [float(n.at[i, 'offset_nm']) for i in path]
        return path

    def chord_arc_in_A(self, i, j):
        return any(jj == j for (jj, *_r) in self.out.get(i, []))

    def real_path(self):
        return list(self.anchors)

    # -- fronteira de Pareto exata (rotulagem em ordem topológica) -------
    def pareto(self):
        labels = defaultdict(list)
        labels[self.o] = [(0.0, 0.0, None, None)]
        for n in self.order:
            L = labels[n]
            if not L:
                continue
            L.sort(key=lambda x: (x[0], x[1]))
            keep, best = [], np.inf
            for lab in L:
                if lab[1] < best - 1e-9:
                    keep.append(lab); best = lab[1]
            labels[n] = keep
            for (j, t, h, d, r, k) in self.out.get(n, []):
                for idx, (T, H, _, _) in enumerate(keep):
                    labels[j].append((T + t, H + h, n, idx))
        self.labels = labels
        front = []
        for idx, (T, H, _, _) in enumerate(labels[self.s]):
            front.append(dict(T=T, H=H, path=self._reconstruct(idx)))
        return front

    def _reconstruct(self, idx):
        node, path = self.s, [self.s]
        while True:
            T, H, pn, pi = self.labels[node][idx]
            if pn is None:
                break
            path.append(pn); node, idx = pn, pi
        return path[::-1]

    # -- contagem de caminhos --------------------------------------------
    def count_paths(self):
        cnt = defaultdict(int); cnt[self.o] = 1
        for n in self.order:
            for (j, t, h, d, r, k) in self.out.get(n, []):
                cnt[j] += cnt[n]
        return cnt[self.s]

    # -- enumeração exaustiva (verificação independente) -----------------
    def enumerate_all(self, limit=6_000_000):
        if self.count_paths() > limit:
            return None
        T = {self.o: np.zeros(1)}; H = {self.o: np.zeros(1)}
        inc = defaultdict(list)
        for n in self.order:
            for (j, t, h, d, r, k) in self.out.get(n, []):
                inc[j].append((n, t, h))
        for n in self.order:
            if n == self.o:
                continue
            if n not in inc:
                continue
            Ts = [T[p] + t for (p, t, h) in inc[n] if p in T]
            Hs = [H[p] + h for (p, t, h) in inc[n] if p in T]
            if Ts:
                T[n] = np.concatenate(Ts); H[n] = np.concatenate(Hs)
        return T[self.s], H[self.s]

    # -- escalarização (soma ponderada normalizada) ----------------------
    def scalar_dp(self, cT, cH):
        """min sum(cT*t + cH*h) sobre caminhos; devolve caminho"""
        best = {self.o: 0.0}; pred = {}
        for n in self.order:
            if n not in best:
                continue
            for (j, t, h, d, r, k) in self.out.get(n, []):
                v = best[n] + cT * t + cH * h
                if v < best.get(j, np.inf) - 1e-12:
                    best[j] = v; pred[j] = n
        path, node = [self.s], self.s
        while node != self.o:
            node = pred[node]; path.append(node)
        return path[::-1]


def lower_hull_mask(front, tol=1e-9):
    """front ordenada por T crescente (H decrescente). True = ponto suportado
    (na envoltória convexa inferior; alcançável por soma ponderada com pesos > 0)"""
    P = [(f['T'], f['H']) for f in front]
    n = len(P)
    if n <= 2:
        return [True] * n
    hull = []
    for p in P:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            cross = (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1)
            if cross < -tol:          # vira para "baixo": mantém côncavo para cima -> remove
                hull.pop()
            else:
                break
        hull.append(p)
    hs = set(hull)
    return [p in hs for p in P]
