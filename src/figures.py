# -*- coding: utf-8 -*-
import sys, pathlib, json, pickle, numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from global_land_mask import globe
from load_ships import *

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.6, 'axes.titlesize': 8,
                     'axes.labelsize': 7.5, 'xtick.labelsize': 6.8, 'ytick.labelsize': 6.8, 'legend.fontsize': 6.6})
RES = json.load(open(OUTDIR / 'results.json', encoding='utf-8'))
G = get_grid()
NAMES = {'E1': 'Embarcação 1', 'E2': 'Embarcação 2', 'E3': 'Embarcação 3'}
C_REAL, C_T, C_H, C_EQ = '#c0392b', '#e0a800', '#1f77b4', '#2e8b57'
OUT = str(FIGDIR) + '/'
import os; os.makedirs(OUT, exist_ok=True)


def fig_pareto():
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 5.0), gridspec_kw={'height_ratios': [1.0, 1.05]})
    for c, k in enumerate(('E1', 'E2', 'E3')):
        R = RES[k]; p = R['real']['proj']; Tp, Hp = p['time'], p['rh']
        z = np.load(OUTDIR / f'cloud_{k}.npz')
        rng = np.random.default_rng(1); idx = rng.choice(len(z['T']), size=min(60000, len(z['T'])), replace=False)
        x = 100 * (z['T'][idx] / Tp - 1); y = 100 * (z['H'][idx] / Hp - 1)
        F = R['frontier']; fx = [100 * (f['T'] / Tp - 1) for f in F]; fy = [100 * (f['H'] / Hp - 1) for f in F]
        ax = axs[0, c]
        ax.scatter(x, y, s=1.0, c='#b9bec4', alpha=0.35, linewidths=0, rasterized=True)
        ax.plot(fx, fy, '-', color=C_H, lw=0.9, zorder=3); ax.scatter(fx, fy, s=7, color=C_H, zorder=4)
        ax.scatter([0], [0], marker='*', s=55, color=C_REAL, zorder=5, edgecolor='k', linewidths=0.4)
        ax.axhline(0, color='k', lw=0.4, ls=':'); ax.axvline(0, color='k', lw=0.4, ls=':')
        ax.set_title(f'{NAMES[k]}', fontweight='bold')
        ax.set_xlabel('Δ tempo, % vs. real'); 
        if c == 0: ax.set_ylabel('Δ risco-hora, % vs. real')
        ax.text(0.03, 0.94, f'{R["enumeration"]["n_paths"]:,} caminhos'.replace(',', '.'), transform=ax.transAxes, fontsize=6.2, color='#555', va='top')
        # zoom
        ax = axs[1, c]
        sup = [f['supported'] for f in F]
        ax.plot(fx, fy, '-', color=C_H, lw=0.8, alpha=0.7, zorder=2)
        ax.scatter([a for a, s in zip(fx, sup) if s], [b for b, s in zip(fy, sup) if s], s=16, color=C_H, zorder=4)
        ax.scatter([a for a, s in zip(fx, sup) if not s], [b for b, s in zip(fy, sup) if not s], s=16, facecolor='white', edgecolor=C_H, linewidths=0.9, zorder=4)
        W = R['weighted']
        ax.scatter([100 * (R['payoff']['minT']['T'] / Tp - 1)], [100 * (R['payoff']['minT']['H'] / Hp - 1)], marker='^', s=38, color=C_T, edgecolor='k', linewidths=0.4, zorder=6)
        ax.scatter([100 * (R['payoff']['minH']['T'] / Tp - 1)], [100 * (R['payoff']['minH']['H'] / Hp - 1)], marker='v', s=38, color=C_H, edgecolor='k', linewidths=0.4, zorder=6)
        ax.scatter([100 * (W['equilibrado']['T'] / Tp - 1)], [100 * (W['equilibrado']['H'] / Hp - 1)], marker='D', s=28, color=C_EQ, edgecolor='k', linewidths=0.4, zorder=6)
        ax.scatter([0], [0], marker='*', s=70, color=C_REAL, zorder=7, edgecolor='k', linewidths=0.4)
        ax.axhline(0, color='k', lw=0.4, ls=':'); ax.axvline(0, color='k', lw=0.4, ls=':')
        xmin, xmax = min(fx) - 0.5, 0.6; ymin, ymax = min(fy) - 0.8, 0.9
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
        ax.set_xlabel('Δ tempo, % vs. real')
        if c == 0: ax.set_ylabel('Δ risco-hora, % vs. real')
        ax.set_title(f'detalhe: {len(F)} pontos ({R["frontier_summary"]["n_supported"]} suportados)', fontsize=7, style='italic')
    handles = [Line2D([0], [0], marker='*', color='w', markerfacecolor=C_REAL, markeredgecolor='k', markersize=9, label='rota real (projetada em A)'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor=C_H, markersize=5, label='fronteira: ponto suportado'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor=C_H, markersize=5, label='fronteira: não suportado'),
               Line2D([0], [0], marker='^', color='w', markerfacecolor=C_T, markeredgecolor='k', markersize=6, label='tempo mínimo'),
               Line2D([0], [0], marker='D', color='w', markerfacecolor=C_EQ, markeredgecolor='k', markersize=5, label='compromisso (0,5/0,5)'),
               Line2D([0], [0], marker='v', color='w', markerfacecolor=C_H, markeredgecolor='k', markersize=6, label='risco-hora mínimo'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor='#b9bec4', markersize=4, label='demais caminhos')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig.savefig(OUT + 'fig_pareto.png', dpi=300); plt.close(fig)


def coast_axes(ax, lon0, lon1, lat0, lat1, step):
    lats = np.arange(lat0, lat1 + 1e-9, step); lons = np.arange(lon0, lon1 + 1e-9, step)
    LA, LO = np.meshgrid(lats, lons, indexing='ij'); land = globe.is_land(LA, LO)
    ax.set_facecolor('#dfeaf3')
    ax.contourf(LO, LA, land.astype(float), levels=[0.5, 1.5], colors=['#c9c2b0'], zorder=1)
    ax.contour(LO, LA, land.astype(float), levels=[0.5], colors=['#7a725f'], linewidths=0.35, zorder=2)
    ax.set_xlim(lon0, lon1); ax.set_ylim(lat0, lat1)


def route_ll(S, path):
    n = S.nodes
    return [n.at[i, 'Longitude'] for i in path], [n.at[i, 'Latitude'] for i in path]


def fig_maps():
    fig = plt.figure(figsize=(7.0, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.93], hspace=0.22, wspace=0.14)
    boxes = {'E1': (-46, -33, -24.5, 0), 'E2': (-49, -34, -25.5, 0.5), 'E3': (-45, -33, -24, -7)}
    labels = {'E1': 'Embarcação 1', 'E2': 'Embarcação 2', 'E3': 'Embarcação 3'}
    for k, sp in zip(('E1', 'E2', 'E3'), (gs[0, 0], gs[0, 1], gs[1, 0])):
        ax = fig.add_subplot(sp); S = load(k); R = RES[k]
        b = boxes[k]; coast_axes(ax, b[0], b[1], b[2], b[3], 0.05)
        for name, col, ls, lw in (('real_proj', C_REAL, '--', 1.0), ('minT', C_T, '-', 1.0), ('minH', C_H, '-', 1.0)):
            path = S.real_projected() if name == 'real_proj' else (R['frontier'][0]['path'] if name == 'minT' else R['frontier'][-1]['path'])
            lo, la = route_ll(S, path); ax.plot(lo, la, ls, color=col, lw=lw, zorder=5)
        o, s = S.o, S.s
        ax.scatter([S.nodes.at[o, 'Longitude']], [S.nodes.at[o, 'Latitude']], marker='^', s=30, color='#27ae60', edgecolor='k', linewidths=0.5, zorder=7)
        ax.scatter([S.nodes.at[s, 'Longitude']], [S.nodes.at[s, 'Latitude']], marker='v', s=30, color='#8e44ad', edgecolor='k', linewidths=0.5, zorder=7)
        ax.set_title(labels[k], fontweight='bold'); ax.set_xlabel('Longitude (°)'); ax.set_ylabel('Latitude (°)')
        ax.set_aspect(1.0 / np.cos(np.radians(0.5 * (b[2] + b[3]))))
        ax.set_xticks(np.arange(np.ceil(b[0] / 5) * 5, b[1] + 0.1, 5)); ax.grid(True, lw=0.25, alpha=0.4)
        if k == 'E3':
            # retângulo do detalhe
            ax.plot([-43.5, -39.8, -39.8, -43.5, -43.5], [-23.3, -23.3, -21.5, -21.5, -23.3], color='k', lw=0.6, zorder=8)
    # detalhe Cabo Frio - Macaé (Embarcação 3)
    ax = fig.add_subplot(gs[1, 1]); S = load('E3'); R = RES['E3']
    lon0, lon1, lat0, lat1 = -43.5, -39.8, -23.3, -21.5
    coast_axes(ax, lon0, lon1, lat0, lat1, 0.01)
    lats = np.arange(lat0, lat1 + 1e-9, 0.02); lons = np.arange(lon0, lon1 + 1e-9, 0.02)
    LA, LO = np.meshgrid(lats, lons, indexing='ij')
    clr = np.array([[G.clearance(la, lo)[0] for lo in lons[::3]] for la in lats[::3]])
    ax.contour(lons[::3], lats[::3], clr, levels=[3.0], colors=['#b03a2e'], linewidths=0.7, linestyles='-', zorder=3)
    for name, col, ls in (('real_proj', C_REAL, '--'), ('minT', C_T, '-'), ('minH', C_H, '-')):
        path = S.real_projected() if name == 'real_proj' else (R['frontier'][0]['path'] if name == 'minT' else R['frontier'][-1]['path'])
        lo, la = route_ll(S, path); ax.plot(lo, la, ls, color=col, lw=1.1, zorder=5)
    n = S.nodes; q = n[n.queried]
    ax.scatter(q.Longitude, q.Latitude, s=5, color='#444', zorder=6)
    ax.scatter([n.at[S.o, 'Longitude']], [n.at[S.o, 'Latitude']], marker='^', s=30, color='#27ae60', edgecolor='k', linewidths=0.5, zorder=7)
    for nm, (la, lo, ha, va) in {'Búzios': (-22.75, -41.86, 'left', 'bottom'), 'Macaé': (-22.37, -41.77, 'left', 'bottom'),
                                 'Cabo Frio': (-22.9, -42.05, 'right', 'top'), 'Rio de Janeiro': (-22.93, -43.0, 'left', 'top')}.items():
        ax.text(lo, la, nm, fontsize=5.6, color='#222', zorder=9, ha=ha, va=va)
    ax.set_title('Detalhe: Rio – Cabo Frio – Macaé (Embarcação 3)', fontsize=7.4, fontweight='bold')
    ax.set_xlabel('Longitude (°)'); ax.set_ylabel('Latitude (°)'); ax.set_aspect(1.0 / np.cos(np.radians(22.4)))
    ax.grid(True, lw=0.25, alpha=0.4)
    handles = [Line2D([0], [0], color=C_REAL, lw=1.2, ls='--', label='rota real projetada em A'),
               Line2D([0], [0], color=C_T, lw=1.4, label='tempo mínimo'), Line2D([0], [0], color=C_H, lw=1.4, label='risco-hora mínimo'),
               Line2D([0], [0], color='#b03a2e', lw=1, label='folga de 3 nm da costa (detalhe)'),
               Line2D([0], [0], marker='^', color='w', markerfacecolor='#27ae60', markeredgecolor='k', markersize=6, label='origem'),
               Line2D([0], [0], marker='v', color='w', markerfacecolor='#8e44ad', markeredgecolor='k', markersize=6, label='destino')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(OUT + 'fig_mapas.png', dpi=300, bbox_inches='tight'); plt.close(fig)


def fig_offsets():
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.3), sharey=True)
    for ax, k in zip(axs, ('E1', 'E2', 'E3')):
        S = load(k); R = RES[k]; n = S.nodes
        for name, col, ls in (('minT', C_T, '-'), ('minH', C_H, '-')):
            path = R['frontier'][0]['path'] if name == 'minT' else R['frontier'][-1]['path']
            pts = [(n.at[i, 'coluna'], n.at[i, 'offset_nm']) for i in path if n.at[i, 'queried']]
            ax.step([p[0] for p in pts], [p[1] for p in pts], where='mid', color=col, lw=1.2, ls=ls, marker='o', ms=3, label='tempo mínimo' if name == 'minT' else 'risco-hora mínimo')
        pj = S.real_projected(); pts = [(n.at[i, 'coluna'], n.at[i, 'offset_nm']) for i in pj if n.at[i, 'queried']]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], '--', color=C_REAL, lw=1.0, label='rota real projetada')
        ax.set_title(NAMES[k], fontweight='bold'); ax.set_xlabel('coluna'); ax.set_ylim(-36, 36); ax.set_yticks([-30, -15, 0, 15, 30]); ax.grid(True, lw=0.3, alpha=0.5)
    axs[0].set_ylabel('deslocamento lateral (nm)'); axs[0].legend(loc='lower left', frameon=False, fontsize=6)
    fig.tight_layout(); fig.savefig(OUT + 'fig_offsets.png', dpi=300); plt.close(fig)


if __name__ == '__main__':
    fig_pareto(); fig_offsets(); fig_maps(); print('figuras geradas em', OUT)
