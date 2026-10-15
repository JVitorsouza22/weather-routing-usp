# -*- coding: utf-8 -*-
"""Teste de sensibilidade: devolve à malha os nós de mar aberto que tinham sido retirados na triagem inicial.

Contexto. Dos 150 nós candidatos (10 colunas x 5 deslocamentos x 3 viagens), 18 ficaram fora da malha:
13 em terra, 1 com folga da costa menor que 3 nm e 4 em mar aberto (folga de 19 a 59 nm). A triagem inicial foi
manual e o motivo de retirar esses 4 nós de mar aberto não foi registrado. Este script refaz toda a análise com
esses 4 nós de volta (mais 1, de folga 2,4 nm, que a regra de folga descarta de qualquer forma) e compara com
o caso principal do artigo (Tabela 7, linha "Com os nós de mar aberto").

Os 5 nós estão em data/nos_retirados.csv (mesmo formato dos arquivos E*_nos_consultados.csv).
O teste roda numa cópia temporária do projeto e não altera nada em out/.

Uso (depois de rodar_tudo.py):
    python src/teste_nos_retirados.py
"""
import json, pathlib, shutil, subprocess, sys, tempfile
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
NAVIOS = ('E1', 'E2', 'E3')


def rodar(tmp):
    for script in ('src/make_mesh.py', 'src/analysis.py'):
        r = subprocess.run([sys.executable, script], cwd=tmp, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-2000:]); raise SystemExit(f'falhou: {script}')


def resumo(R, k):
    f = R[k]['frontier']; s = R[k]['sens']['principal']; fs = R[k]['frontier_summary']
    real = R[k]['real']['proj']
    domina = all(x['T'] < real['time'] and x['H'] < real['rh'] for x in f)
    dts = [100 * (x['T'] / real['time'] - 1) for x in f]
    return dict(n=fs['n'], sup=fs['n_supported'], dT_minT=s['minT']['dT'], dH_minT=s['minT']['dH'],
                dT_minH=s['minH']['dT'], dH_minH=s['minH']['dH'], domina=domina, menor_reducao=max(dts),
                nos=R[k]['mesh']['nodes_q'], cortados=R[k]['mesh']['nodes_dropped'])


def main():
    base_file = ROOT / 'out' / 'results.json'
    if not base_file.exists():
        raise SystemExit('Rode antes:  python rodar_tudo.py')
    extra = pd.read_csv(ROOT / 'data' / 'nos_retirados.csv')
    tmp = pathlib.Path(tempfile.mkdtemp(prefix='wr_teste_'))
    try:
        shutil.copytree(ROOT / 'src', tmp / 'src', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy(ROOT / 'config.json', tmp / 'config.json')
        shutil.copytree(ROOT / 'data', tmp / 'data')
        (tmp / 'out').mkdir()
        grade = ROOT / 'out' / 'clearance_grid.pkl'
        if grade.exists():
            shutil.copy(grade, tmp / 'out' / 'clearance_grid.pkl')
        for k in NAVIOS:
            f = tmp / 'data' / f'{k}_nos_consultados.csv'
            d = pd.read_csv(f)
            novos = extra[extra.navio == k].drop(columns='navio')
            pd.concat([d, novos], ignore_index=True).to_csv(f, index=False)
            print(f'{k}: {len(d)} nós no caso principal + {len(novos)} devolvidos = {len(d) + len(novos)}')
        print('Refazendo a malha e a análise numa cópia temporária (alguns segundos)...')
        rodar(tmp)
        A = json.load(open(base_file, encoding='utf-8'))
        B = json.load(open(tmp / 'out' / 'results.json', encoding='utf-8'))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print()
    print('Variação em relação à rota real projetada (%), soluções extremas da fronteira')
    print(f'{"":4}{"caso":22}{"soluções (suport.)":20}{"tempo mín.: tempo / risco-h":29}{"risco-h mín.: tempo / risco-h":31}{"todas dominam":14}')
    tot = {'A': [0, 0], 'B': [0, 0]}
    for k in NAVIOS:
        for nome, R in (('principal', A), ('com nós de mar aberto', B)):
            x = resumo(R, k)
            print(f'{k:4}{nome:22}{str(x["n"]) + " (" + str(x["sup"]) + ")":20}'
                  f'{x["dT_minT"]:+6.1f} / {x["dH_minT"]:+6.1f}{"":14}{x["dT_minH"]:+6.1f} / {x["dH_minH"]:+6.1f}{"":16}{"sim" if x["domina"] else "NÃO"}')
            t = tot['A' if R is A else 'B']; t[0] += x['n']; t[1] += x['sup']
    print()
    print(f'Total de soluções eficientes (suportadas): principal {tot["A"][0]} ({tot["A"][1]}); com os nós de volta {tot["B"][0]} ({tot["B"][1]})')
    menor = min(abs(resumo(B, k)['menor_reducao']) for k in NAVIOS)
    print(f'Menor redução de tempo entre todas as soluções, com os nós de volta: {menor:.1f}%')
    saida = {'principal': {k: resumo(A, k) for k in NAVIOS}, 'com_nos_de_mar_aberto': {k: resumo(B, k) for k in NAVIOS}}
    json.dump(saida, open(ROOT / 'out' / 'teste_nos_retirados.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('Resumo gravado em out/teste_nos_retirados.json')


if __name__ == '__main__':
    main()
