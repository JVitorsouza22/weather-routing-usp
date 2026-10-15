# -*- coding: utf-8 -*-
"""
Construção do conjunto de arcos A por regra (Seção 4.1 do artigo).
  - Nós: nós da malha com clima consultado (colunas inteiras) que respeitam folga >= D_MIN da costa
         (nós de origem/destino e nós a <= 15 nm deles são isentos, mas não podem estar em terra);
  - Para cada par (i,j) de colunas consecutivas:
        * corda (segmento reto) válida  -> arco direto;
        * corda inválida                -> menor caminho na grade oceânica com folga >= D_PLAN,
                                           simplificado em waypoints com folga >= D_MIN; aceito se
                                           comprimento <= MAX_RATIO * corda.
  Máscara terra-mar: GLOBE 1 km (global-land-mask). Folga: transformada de distância euclidiana.
"""
import numpy as np, pandas as pd, pickle
from global_land_mask import globe
from scipy import ndimage
from skimage.graph import route_through_array
from wr_model import haversine
from geo_check import sample_segment

from paths import CFG
D_MIN = CFG['D_MIN_NM']        # folga mínima da costa (nm) exigida dos arcos
D_PLAN = CFG['D_PLAN_NM']       # folga usada no planejamento (margem sobre D_MIN, > resolução da máscara)
PORT_TOL = CFG['PORT_TOL_NM']      # nm ao redor dos pontos de referência de porto onde a máscara de 1 km não é confiável
PORT_MARGIN = CFG['PORT_MARGIN_NM'] # nm: aproximações portuárias isentas da regra de folga (não podem ser terra)
MAX_RATIO = CFG['MAX_RATIO']   # desvio máximo aceito em relação à corda
RES = CFG['RES_GRAU']         # grau (~0,6 nm): resolução nativa da máscara GLOBE (1 km)


class Mesh:
    def __init__(self, G):
        self.G = G                      # ClearGrid (folga em nm)

    # ------------------------------------------------------------------ utilidades
    def seg_ok(self, la1, lo1, la2, lo2, near_ports, dmin=D_MIN):
        la, lo, dd = sample_segment(la1, lo1, la2, lo2, 1.0)
        c = self.G.clearance(la, lo)
        land = globe.is_land(la, lo)
        exempt = np.zeros(len(la), bool); tol = np.zeros(len(la), bool)
        for (pl, po) in near_ports:
            dd_ = haversine(la, lo, pl, po)
            exempt |= dd_ <= PORT_MARGIN
            tol |= dd_ <= PORT_TOL                 # pixel do porto: máscara 1 km não resolve estuários
        ok_clear = (c >= dmin) | exempt
        ok_land = ~land | tol
        # dentro da margem portuária a máscara de 1 km pode marcar terra em baías estreitas: tolera terra
        # apenas nos primeiros/últimos 3 nm do próprio nó portuário
        return bool(ok_clear.all() and ok_land.all()), float(c[~exempt].min()) if (~exempt).any() else 99.0

    def plan(self, la1, lo1, la2, lo2, near_ports):
        """menor caminho na grade oceânica com folga >= D_PLAN entre dois pontos; devolve polilinha (lat,lon) ou None"""
        pad = 1.2
        lat0, lat1 = min(la1, la2) - pad, max(la1, la2) + pad
        lon0, lon1 = min(lo1, lo2) - pad, max(lo1, lo2) + pad
        lats = np.arange(lat0, lat1 + 1e-9, RES); lons = np.arange(lon0, lon1 + 1e-9, RES)
        LA, LO = np.meshgrid(lats, lons, indexing='ij')
        land = globe.is_land(LA, LO)
        mid = 0.5 * (lat0 + lat1)
        clr = ndimage.distance_transform_edt(~land, sampling=(RES * 60, RES * 60 * np.cos(np.radians(mid))))
        allowed = (clr >= D_PLAN)
        for (pl, po) in near_ports:                       # nos portos basta ser água (tolerância no pixel do porto)
            dd_ = haversine(LA, LO, pl, po)
            allowed |= ((dd_ <= PORT_MARGIN) & ~land) | (dd_ <= PORT_TOL)
        for (pl, po) in ((la1, lo1), (la2, lo2)):         # pontos extremos: disco de 2 nm se for água
            allowed |= (haversine(LA, LO, pl, po) <= 2.0) & ~land
        cost = np.where(allowed, 1.0, np.inf)
        s = (int(round((la1 - lats[0]) / RES)), int(round((lo1 - lons[0]) / RES)))
        e = (int(round((la2 - lats[0]) / RES)), int(round((lo2 - lons[0]) / RES)))
        if not (np.isfinite(cost[s]) and np.isfinite(cost[e])):
            return None
        try:
            idx, c = route_through_array(cost, s, e, fully_connected=True, geometric=True)
        except Exception:
            return None
        if not np.isfinite(c):
            return None
        pts = np.array([(lats[a], lons[b]) for a, b in idx])
        pts[0] = (la1, lo1); pts[-1] = (la2, lo2)
        return self._simplify(pts, near_ports)

    def _simplify(self, pts, near_ports):
        # Douglas-Peucker em plano local (nm) com verificação de folga dos segmentos simplificados
        def dp(P, tol):
            mid = np.radians(P[:, 0].mean())
            xy = np.c_[(P[:, 1] - P[0, 1]) * 60 * np.cos(mid), (P[:, 0] - P[0, 0]) * 60]
            keep = np.zeros(len(P), bool); keep[0] = keep[-1] = True
            stack = [(0, len(P) - 1)]
            while stack:
                a, b = stack.pop()
                if b <= a + 1:
                    continue
                v = xy[b] - xy[a]; L = np.hypot(*v)
                d = (np.abs(np.cross(v, xy[a + 1:b] - xy[a])) / L) if L > 0 else np.hypot(*(xy[a + 1:b] - xy[a]).T)
                m = int(np.argmax(d)); 
                if d[m] > tol:
                    k = a + 1 + m; keep[k] = True; stack += [(a, k), (k, b)]
            return P[keep]
        for tol in (4.0, 2.5, 1.5, 0.8, 0.4):
            Q = dp(pts, tol)
            ok = all(self.seg_ok(Q[i, 0], Q[i, 1], Q[i + 1, 0], Q[i + 1, 1], near_ports)[0] for i in range(len(Q) - 1))
            if ok:
                return Q
        return pts[::4] if len(pts) > 8 else pts
