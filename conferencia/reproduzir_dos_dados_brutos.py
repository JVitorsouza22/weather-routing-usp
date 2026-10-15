# -*- coding: utf-8 -*-
"""Refaz, a partir dos dados AIS brutos, tudo o que vem ANTES da malha e compara com o que está em data/.

Este script só funciona com os arquivos brutos, que NÃO são distribuídos (restrição contratual com o provedor
dos dados AIS). Quem tem os arquivos pode rodá-lo para ver que os números de data/ saem direto dos dados brutos.

Prepare uma pasta com uma subpasta por viagem (E1, E2, E3), cada uma com 3 arquivos, com estes nomes:

    dados_brutos/E1/linha_do_tempo.csv   exportação do sistema de monitoramento (separador ';'; colunas
                                         'Sample time' e 'Position time')
    dados_brutos/E1/janela.csv           trecho da viagem, de hora em hora, com o clima da posição
                                         (colunas Sample time, Latitude, Longitude, Speed, ...)
    dados_brutos/E1/nos.csv              nós consultados, com o clima (colunas id, Latitude, Longitude,
                                         Sample time, vento_vel_kmh, onda_altura_m, corrente_vel_kmh, ...)

Uso:
    python conferencia/reproduzir_dos_dados_brutos.py                       # lê ./dados_brutos
    python conferencia/reproduzir_dos_dados_brutos.py --pasta C:\\caminho\\dados_brutos

O que é conferido
  1. Números da Tabela 2 (amostras, SOG médio, distância do rastro, defasagem da posição) contra data/parametros.json.
  2. Os 150 nós candidatos (10 colunas x 5 deslocamentos x 3 viagens), refeitos a partir do rastro, contra nos.csv.
  3. O clima e o índice de risco dos nós contra data/E*_nos_consultados.csv.
  4. A contagem dos nós: quantos ficaram fora da malha e por quê (terra, folga menor que 3 nm, mar aberto).
"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from paths import DATA, CFG, get_grid
from wr_model import haversine, risk_index, R_NM
from global_land_mask import globe

NAVIOS = ('E1', 'E2', 'E3')
OFFSETS = (-30, -15, 0, 15, 30)       # nm; positivo = a boreste (direita) do rumo
N_COL = 10


def rumo_graus(lat1, lon1, lat2, lon2):
    p1, p2, dl = map(math.radians, (lat1, lat2, lon2 - lon1))
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def ponto_destino(lat, lon, rumo, dist_nm):
    la, lo, b = map(math.radians, (lat, lon, rumo)); d = dist_nm / R_NM
    la2 = math.asin(math.sin(la) * math.cos(d) + math.cos(la) * math.sin(d) * math.cos(b))
    lo2 = lo + math.atan2(math.sin(b) * math.sin(d) * math.cos(la), math.cos(d) - math.sin(la) * math.sin(la2))
    return math.degrees(la2), math.degrees(lo2)


def tabela1(pasta, P):
    """escalares da Tabela 2, recalculados da janela e da linha do tempo"""
    w = pd.read_csv(pasta / 'janela.csv'); w['Sample time'] = pd.to_datetime(w['Sample time'], utc=True)
    t = pd.read_csv(pasta / 'linha_do_tempo.csv', sep=';', encoding='utf-8-sig')
    t['Sample time'] = pd.to_datetime(t['Sample time'], utc=True); t['Position time'] = pd.to_datetime(t['Position time'], utc=True)
    m = w[['Sample time']].merge(t[['Sample time', 'Position time']], on='Sample time', how='left')
    lag = ((m['Sample time'] - m['Position time']).dt.total_seconds() / 3600)
    nav = w[(w.Speed >= 3) & (w.Speed <= 20)]
    dens = float(haversine(w.Latitude.values[:-1], w.Longitude.values[:-1], w.Latitude.values[1:], w.Longitude.values[1:]).sum())
    novo = dict(v_obs=nav.Speed.mean(), v_med=nav.Speed.median(), n_nav=len(nav), n=len(w),
                n_invalid=int((w.Speed >= 102).sum()),          # 102,3 = velocidade indisponível no AIS
                lag_med=lag.median(), lag_med_abs=lag.abs().median(), lag_p90=lag.abs().quantile(0.9), lag_max=lag.abs().max(),
                pct_stale=100 * (lag.abs() > 3).mean(),
                dur_h=(w['Sample time'].iloc[-1] - w['Sample time'].iloc[0]).total_seconds() / 3600,
                dense_dist_nm=dens)
    ok = True
    for k, v in novo.items():
        d = abs(float(v) - float(P[k]))
        marca = 'ok' if d < 1e-6 else 'DIFERENTE'
        ok &= (d < 1e-6)
        print(f'   {k:14s} bruto={float(v):14.6f}  data/={float(P[k]):14.6f}  {marca}')
    return ok, w


def candidatos(w):
    """os 50 nós candidatos de uma viagem: 10 pontos de referência (amostras igualmente espaçadas na janela)
    x 5 deslocamentos laterais, perpendiculares ao rumo (rumo = do ponto de referência ao seguinte; no último, do anterior)"""
    n = len(w)
    idx = [int(round(i * (n - 1) / (N_COL - 1))) for i in range(N_COL)]
    A = w.iloc[idx][['Latitude', 'Longitude']].values
    linhas = []
    for c in range(N_COL):
        rumo = rumo_graus(*A[c], *A[c + 1]) if c < N_COL - 1 else rumo_graus(*A[c - 1], *A[c])
        for j, off in enumerate(OFFSETS):
            if off == 0:
                lat, lon = A[c]
            else:
                lat, lon = ponto_destino(A[c][0], A[c][1], rumo + (90 if off > 0 else -90), abs(off))
            linhas.append(dict(id=c * 5 + j, coluna=c, offset_nm=off, Latitude=lat, Longitude=lon, t_rel_h=float(idx[c])))
    return pd.DataFrame(linhas)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pasta', default='dados_brutos', help='pasta com E1/, E2/ e E3/ (padrão: ./dados_brutos)')
    a = ap.parse_args()
    base = pathlib.Path(a.pasta)
    if not (base / 'E1').exists():
        raise SystemExit(f'Não achei {base}/E1. Leia as instruções no início deste arquivo.')
    P = json.load(open(DATA / 'parametros.json', encoding='utf-8'))
    G = get_grid(); D_MIN = CFG['D_MIN_NM']
    tudo_ok = True
    resumo = dict(candidatos=0, consultados=0, na_malha=0, terra=0, folga=0, mar_aberto=0)
    for k in NAVIOS:
        print(f'\n===== {k} =====')
        pasta = base / k
        print('1) Tabela 2: escalares recalculados dos dados brutos')
        ok1, w = tabela1(pasta, P[k]); tudo_ok &= ok1

        print('2) nós candidatos refeitos a partir do rastro, contra nos.csv')
        C = candidatos(w)
        nos = pd.read_csv(pasta / 'nos.csv')
        m = nos.merge(C, on='id', suffixes=('', '_novo'))
        erro_nm = float(np.hypot(m.Latitude - m.Latitude_novo, (m.Longitude - m.Longitude_novo) * np.cos(np.radians(m.Latitude))).max() * 60)
        print(f'   {len(nos)} nós consultados em nos.csv; erro máximo de posição = {erro_nm:.2e} nm  ->', 'ok' if erro_nm < 1e-6 else 'DIFERENTE')
        tudo_ok &= erro_nm < 1e-6

        print('3) clima e índice de risco contra data/')
        dat = pd.read_csv(DATA / f'{k}_nos_consultados.csv')
        mm = dat.merge(nos, on='id', suffixes=('', '_bruto'))
        dif_clima = max(float((mm[c] - mm[c + '_bruto']).abs().max()) for c in
                        ('Latitude', 'Longitude', 'vento_vel_kmh', 'onda_altura_m', 'corrente_vel_kmh', 'corrente_dir_graus'))
        R_novo = risk_index(mm.onda_altura_m, mm.vento_vel_kmh, mm.corrente_vel_kmh)
        dif_R = float((R_novo - mm.R).abs().max())
        t_nos = (pd.to_datetime(mm['Sample time'], utc=True) - w['Sample time'].iloc[0]).dt.total_seconds() / 3600
        dif_t = float((t_nos - mm.t_rel_h).abs().max())
        print(f'   {len(dat)} nós em data/ ({len(mm)} achados em nos.csv); diferença máxima: clima={dif_clima:.1e}  risco R={dif_R:.1e}  tempo relativo={dif_t:.1e} h')
        tudo_ok &= (dif_clima < 1e-9 and dif_R < 1e-9 and dif_t < 1e-6)

        print('4) contagem dos nós candidatos')
        C['terra'] = globe.is_land(C.Latitude.values, C.Longitude.values)
        C['folga'] = G.clearance(C.Latitude.values, C.Longitude.values)
        C['consultado'] = C.id.isin(nos.id)
        C['na_malha'] = C.id.isin(dat.id)
        fora = C[~C.na_malha]
        terra = fora[fora.terra]; folga = fora[~fora.terra & (fora.folga < D_MIN)]; mar = fora[~fora.terra & (fora.folga >= D_MIN)]
        print(f'   candidatos={len(C)}  consultados={int(C.consultado.sum())}  na malha={int(C.na_malha.sum())}  fora={len(fora)}'
              f'  (terra={len(terra)}, folga<{D_MIN:g} nm={len(folga)}, mar aberto={len(mar)})')
        if len(mar):
            print('   nós de mar aberto fora da malha (id, coluna, deslocamento, folga em nm):',
                  [(int(r.id), int(r.coluna), int(r.offset_nm), round(float(r.folga), 1)) for r in mar.itertuples()])
        resumo['candidatos'] += len(C); resumo['consultados'] += int(C.consultado.sum()); resumo['na_malha'] += int(C.na_malha.sum())
        resumo['terra'] += len(terra); resumo['folga'] += len(folga); resumo['mar_aberto'] += len(mar)
    print('\n===== TOTAL =====')
    print('   ', resumo)
    print('RESULTADO:', 'tudo confere com data/' if tudo_ok else 'HÁ DIFERENÇAS (veja as linhas marcadas)')


if __name__ == '__main__':
    main()
