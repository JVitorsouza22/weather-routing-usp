# -*- coding: utf-8 -*-
"""Gera, a partir de out/results.json e de data/, os números das tabelas do artigo.

Cada tabela do artigo sai daqui, formatada como no texto (vírgula decimal, ponto nos milhares,
sinal de menos tipográfico). Assim nenhum número é digitado à mão.

Tabelas: 1 exemplo numérico simples (src/exemplo_didatico.py); 2 dados das viagens; 3 malha e arcos;
4 base de comparação; 5 fronteira de Pareto; 6 soluções de referência; 7 sensibilidade; 8 projeção;
A1 a A3 âncoras da rota real (Anexo A).

Uso (depois de rodar_tudo.py):
    python src/tabelas_artigo.py            # imprime todas as tabelas e grava out/tabelas_artigo.txt
    python src/tabelas_artigo.py 6          # só a Tabela 6   (use A para o Anexo A)
"""
import json, pathlib, sys
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
from exemplo_didatico import tabela_exemplo

ROOT = pathlib.Path(__file__).resolve().parent.parent
NAVIOS = ('E1', 'E2', 'E3')
ROTAS = {'E1': 'São Luís/MA → Rio de Janeiro/RJ', 'E2': 'Santos/SP → Belém/PA', 'E3': 'Rio de Janeiro/RJ → Recife/PE'}
MENOS = '−'


def f(x, nd=1, sinal=False):
    """formata número: vírgula decimal, ponto nos milhares, '−' tipográfico; sinal=True força '+' nos positivos."""
    s = f'{abs(x):,.{nd}f}'.replace(',', '#').replace('.', ',').replace('#', '.')
    if round(x, nd) < 0:
        return MENOS + s
    return ('+' + s) if sinal and round(x, nd) > 0 else s


def carregar():
    R = json.load(open(ROOT / 'out' / 'results.json', encoding='utf-8'))
    P = json.load(open(ROOT / 'data' / 'parametros.json', encoding='utf-8'))
    return R, P


def tabela2(R, P):
    cab = ['Item', 'Embarcação 1', 'Embarcação 2', 'Embarcação 3']
    L = []
    L.append(['Rota'] + [ROTAS[k] for k in NAVIOS])
    L.append(['Distância do rastro AIS (nm)'] + [f(R[k]['data']['dense_dist']) for k in NAVIOS])
    L.append(['Amostras horárias na janela da viagem'] + [str(P[k]['n']) for k in NAVIOS])
    L.append(['Amostras em navegação (3 ≤ SOG ≤ 20 nós)'] + [str(P[k]['n_nav']) for k in NAVIOS])
    L.append(['SOG médio reportado em navegação (nós)'] + [f(P[k]['v_obs'], 2) for k in NAVIOS])
    L.append(['Defasagem da posição AIS (valor absoluto): mediana / p90 / máx. (h)'] +
             [f"{f(P[k]['lag_med_abs'])} / {f(P[k]['lag_p90'])} / {f(P[k]['lag_max'])}" for k in NAVIOS])
    L.append(['Amostras com posição defasada em mais de 3 h (%)'] + [f(P[k]['pct_stale']) for k in NAVIOS])
    L.append(['Corrente média ao longo da rota real (nós; + favorável)'] + [f(R[k]['data']['u_bar_real'], 2, True) for k in NAVIOS])
    L.append(['Velocidade de serviço na água, V (nós)'] + [f(R[k]['data']['V_water'], 2) for k in NAVIOS])
    return cab, L


def tabela3(R, P):
    cab = ['Item', 'Embarcação 1', 'Embarcação 2', 'Embarcação 3']
    M = {k: R[k]['mesh'] for k in NAVIOS}
    L = [
        ['Nós candidatos (10 colunas × 5 posições)'] + ['50'] * 3,
        ['Nós com clima consultado mantidos na triagem'] + [str(M[k]['nodes_q']) for k in NAVIOS],
        ['Nós descartados pela regra de folga'] + [str(M[k]['nodes_dropped']) for k in NAVIOS],
        ['Pares de nós avaliados (colunas consecutivas)'] + [str(M[k]['pairs']) for k in NAVIOS],
        ['Cordas válidas (arcos diretos)'] + [str(M[k]['chords_ok']) for k in NAVIOS],
        ['Desvios planejados (cadeias de pontos de passagem)'] + [str(M[k]['planned']) for k in NAVIOS],
        ['Razão desvio/corda: mediana (máx.)'] + [f"{f(M[k]['ratio_med'], 2)} ({f(M[k]['ratio_max'], 2)})" for k in NAVIOS],
        ['Pontos de passagem criados'] + [str(M[k]['waypoints']) for k in NAVIOS],
        ['Nós no grafo (com arcos)'] + [str(M[k]['nodes_total'] - M[k]['nodes_dropped']) for k in NAVIOS],
        ['Arcos no grafo'] + [str(M[k]['arcs']) for k in NAVIOS],
        ['Caminhos origem–destino'] + [f(M[k]['paths'], 0) for k in NAVIOS],
    ]
    return cab, L


def tabela4(R, P):
    cab = ['Embarcação', 'Real: rastro AIS completo (nm)', 'Real: cordas (nm)', 'Real: projetada em A (nm)',
           'Otimizada, tempo mín. (nm)', 'Δ vs. rastro completo (%)', 'Δ vs. projetada (%)']
    L = []
    for i, k in enumerate(NAVIOS, 1):
        r = R[k]; d = r['real']['dense_dist']; p = r['real']['proj']['dist']; o = r['payoff']['minT']['dist']
        L.append([f'Embarcação {i}', f(d), f(r['real']['chords']['dist']), f(p), f(o), f(100 * (o / d - 1), 1, True), f(100 * (o / p - 1), 1, True)])
    return cab, L


def tabela5(R, P):
    cab = ['Embarcação', 'Caminhos', 'Soluções eficientes', 'Suportadas / não suportadas',
           'Tempo: ideal – nadir (h)', 'Risco-hora: ideal – nadir', 'Dominam a rota real']
    L = []
    for i, k in enumerate(NAVIOS, 1):
        r = R[k]; s = r['frontier_summary']; pay = r['payoff']
        L.append([f'Embarcação {i}', f(r['mesh']['paths'], 0), str(s['n']), f"{s['n_supported']} / {s['n_unsupported']}",
                  f"{f(pay['ideal']['T'])} – {f(pay['nadir']['T'])}", f"{f(pay['ideal']['H'], 2)} – {f(pay['nadir']['H'], 2)}",
                  f"{r['epsilon']['dominated_by']} de {s['n']}"])
    return cab, L


def tabela6(R, P):
    cab = ['Embarcação', 'Solução', 'Dist. (nm)', 'Tempo (h)', 'Risco-hora', 'R méd.', 'Δ tempo (%)', 'Δ risco-hora (%)']
    L = []
    for i, k in enumerate(NAVIOS, 1):
        r = R[k]; real = r['real']['proj']
        L.append([f'Embarcação {i}', 'Real (projetada em A)', f(real['dist']), f(real['time']), f(real['rh'], 2), f(real['rbar'], 3), '–', '–'])
        for nome, x in (('Tempo mínimo', r['payoff']['minT']), ('Compromisso (0,5; 0,5)', r['weighted']['equilibrado']),
                        ('Risco-hora mínimo', r['payoff']['minH'])):
            T = x['T']; H = x['H']
            L.append(['', nome, f(x['dist']), f(T), f(H, 2), f(x['rbar'], 3),
                      f(100 * (T / real['time'] - 1), 1, True), f(100 * (H / real['rh'] - 1), 1, True)])
    return cab, L


def tabela7(R, P):
    cab = ['Embarcação', 'Variante', 'Soluções eficientes', 'Tempo mín.: Δ tempo (%)', 'Tempo mín.: Δ risco-hora (%)',
           'Risco-hora mín.: Δ tempo (%)', 'Risco-hora mín.: Δ risco-hora (%)']
    nomes = [('principal', 'Principal (ponderação espacial, com corrente)'), ('linear', 'Interpolação linear nos pontos de passagem'),
             ('sem_corrente', 'Sem efeito da corrente')]
    extra = None
    fx = ROOT / 'out' / 'teste_nos_retirados.json'
    if fx.exists():
        extra = json.load(open(fx, encoding='utf-8'))['com_nos_de_mar_aberto']
    L = []
    for i, k in enumerate(NAVIOS, 1):
        first = True
        for chave, nome in nomes:
            s = R[k]['sens'][chave]
            L.append([f'Embarcação {i}' if first else '', nome, str(s['n']), f(s['minT']['dT'], 1, True), f(s['minT']['dH'], 1, True),
                      f(s['minH']['dT'], 1, True), f(s['minH']['dH'], 1, True)])
            first = False
        if extra:
            x = extra[k]
            L.append(['', 'Com os nós de mar aberto devolvidos à malha', str(x['n']), f(x['dT_minT'], 1, True), f(x['dH_minT'], 1, True),
                      f(x['dT_minH'], 1, True), f(x['dH_minH'], 1, True)])
    return cab, L


def tabela8(R, P):
    cab = ['Embarcação', 'Trânsitos por ano (estimados)', 'Economia por viagem (h)', '3 meses (h)', '6 meses (h)', '12 meses (h)',
           'Redução de tempo por viagem (%)']
    L = []
    for i, k in enumerate(NAVIOS, 1):
        pr = R[k]['projection']; pp = pr['per_period']; real = R[k]['real']['proj']['time']
        red = 100 * (1 - R[k]['payoff']['minT']['T'] / real)
        L.append([f'Embarcação {i}', str(pr['transits_year']), f(pr['dT_trip_minT']), f(pp['3']['hours_minT'], 0), f(pp['6']['hours_minT'], 0),
                  f(pp['12']['hours_minT'], 0), f(red)])
    return cab, L


def anexo_a(R, P, k):
    cab = ['Col.', 'Lat. (°)', 'Lon. (°)', 't (h)', 'Hs (m)', 'Vento (km/h)', 'Corrente (km/h)', 'Dir. corr. (°)', 'R']
    d = pd.read_csv(ROOT / 'data' / f'{k}_nos_consultados.csv')
    d = d[(d.offset_nm == 0) & (d.coluna == d.coluna.round())].sort_values('coluna')
    L = []
    for r in d.itertuples():
        L.append([str(int(r.coluna)), f(r.Latitude), f(r.Longitude), f(r.t_rel_h, 0), f(r.onda_altura_m, 2), f(r.vento_vel_kmh),
                  f(r.corrente_vel_kmh), f(r.corrente_dir_graus, 0), f(r.R, 3)])
    return cab, L


def todas():
    R, P = carregar()
    T = {'1': tabela_exemplo(), '2': tabela2(R, P), '3': tabela3(R, P), '4': tabela4(R, P), '5': tabela5(R, P),
         '6': tabela6(R, P), '7': tabela7(R, P), '8': tabela8(R, P)}
    for i, k in enumerate(NAVIOS, 1):
        T[f'A{i}'] = anexo_a(R, P, k)
    return T


def formatar(nome, cab, L):
    larg = [max(len(str(r[j])) for r in [cab] + L) for j in range(len(cab))]
    linhas = [f'--- Tabela {nome} ---']
    for r in [cab] + L:
        linhas.append('  '.join(str(c).ljust(larg[j]) for j, c in enumerate(r)).rstrip())
    return '\n'.join(linhas) + '\n'


if __name__ == '__main__':
    T = todas()
    alvo = sys.argv[1].upper() if len(sys.argv) > 1 else None
    blocos = [formatar(nome, cab, L) for nome, (cab, L) in T.items()
              if alvo is None or nome == alvo or (alvo == 'A' and nome.startswith('A'))]
    texto = '\n'.join(blocos)
    print(texto)
    if alvo is None:                      # rodando sem argumento, grava também num arquivo de texto (UTF-8)
        destino = ROOT / 'out' / 'tabelas_artigo.txt'
        destino.write_text(texto, encoding='utf-8')
        print(f'Tabelas gravadas também em {destino.relative_to(ROOT)} (texto UTF-8; abra no Bloco de Notas).')
