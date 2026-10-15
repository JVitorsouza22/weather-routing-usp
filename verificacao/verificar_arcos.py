# -*- coding: utf-8 -*-
"""Confere todos os arcos do grafo contra a máscara terra-mar. Rode depois de  python rodar_tudo.py

A regra que cria os arcos (src/build_mesh.py) já exige folga de 3 nm da costa. Este script confere o resultado
de outro jeito, amostra por amostra, para ver se a regra funcionou:

  1. Cada arco é cortado em amostras a cada 1 nm. Nenhuma amostra pode cair em terra, a não ser perto dos portos
     de origem e de destino (a máscara de 1 km às vezes cobre o cais; tolera-se terra até 15 nm do porto).
  2. Para cada amostra fora das aproximações portuárias, mede-se a folga da costa por outro método: procura-se o
     pixel de terra mais perto, em 48 direções, em passos de 0,25 nm, até 6 nm. Informa-se a menor folga e quantas
     amostras ficam abaixo de 3 nm. Como a máscara tem pixels de 1 km (cerca de 0,5 nm de incerteza), uma ou outra
     amostra um pouco abaixo de 3 nm é esperada.

Uso:   python verificacao/verificar_arcos.py        (cerca de 1 a 3 minutos)
"""
import pathlib, sys
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / 'src'))
from load_ships import *                              # noqa: F401,F403  (load, CFG, haversine, ...)
from geo_check import sample_segment, clearance_nm
from global_land_mask import globe

NAVIOS = ('E1', 'E2', 'E3')
MARGEM_PORTO = CFG['PORT_MARGIN_NM']                    # 15 nm
D_MIN = CFG['D_MIN_NM']                                 # 3 nm

for k in NAVIOS:
    S = load(k); n = S.nodes
    porto_o, porto_d = n.loc[S.o], n.loc[S.s]
    LA, LO = [], []
    for r in S.A.itertuples():
        la, lo, _ = sample_segment(n.at[r.i, 'Latitude'], n.at[r.i, 'Longitude'], n.at[r.j, 'Latitude'], n.at[r.j, 'Longitude'], 1.0)
        LA.append(la); LO.append(lo)
    LA = np.concatenate(LA); LO = np.concatenate(LO)
    d_porto = np.minimum(haversine(LA, LO, porto_o.Latitude, porto_o.Longitude), haversine(LA, LO, porto_d.Latitude, porto_d.Longitude))
    terra = globe.is_land(LA, LO)
    miolo = d_porto >= MARGEM_PORTO                     # fora das aproximações portuárias
    clr = clearance_nm(LA[miolo], LO[miolo], rmax=6.0)
    print(f'{k}: {len(S.A)} arcos, {len(LA)} amostras (a cada 1 nm)')
    print(f'   amostras em terra: {int(terra.sum())} no total; {int((terra & miolo).sum())} fora das aproximações portuárias ({MARGEM_PORTO:g} nm)')
    print(f'   folga da costa nas {int(miolo.sum())} amostras fora das aproximações: mínima = {clr.min():.2f} nm; '
          f'abaixo de {D_MIN:g} nm: {int((clr < D_MIN).sum())}; abaixo de 1 nm: {int((clr < 1.0).sum())}')
