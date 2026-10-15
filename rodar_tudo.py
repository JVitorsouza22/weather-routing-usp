# -*- coding: utf-8 -*-
"""Roda todas as etapas do trabalho, uma depois da outra.

Etapas, em ordem:
  1. src/make_mesh.py      monta os nós e os arcos de cada viagem (regra de folga da costa)
  2. src/analysis.py       calcula rota real, fronteira de Pareto, sensibilidades (gera out/results.json)
  3. src/export_gurobi.py  grava o modelo em arquivos, para conferir com o Gurobi
  4. src/figures.py        desenha as figuras (pode ser pulada)
  5. verificacao/verificar_gurobi.py   confere os resultados com o Gurobi (só com --gurobi)

Uso:
    python rodar_tudo.py              # etapas 1 a 4
    python rodar_tudo.py --gurobi     # etapas 1 a 5 (a licença gratuita do gurobipy basta)
    python rodar_tudo.py --sem-figuras
Em Windows, use  py rodar_tudo.py  no lugar de  python rodar_tudo.py.
Os resultados ficam em out/ (results.json, mesh/, fig/); os arquivos do Gurobi, em verificacao/.
"""
import argparse, pathlib, subprocess, sys, time
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
ROOT = pathlib.Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--gurobi', action='store_true', help='inclui a conferência com o Gurobi')
ap.add_argument('--sem-figuras', action='store_true', help='não desenha as figuras')
a = ap.parse_args()


def run(*cmd, cwd=ROOT):
    t0 = time.time(); print('>>', ' '.join(cmd), flush=True)
    subprocess.run([sys.executable, *cmd], cwd=cwd, check=True); print(f'   ok ({time.time() - t0:.0f} s)', flush=True)


run('src/make_mesh.py'); run('src/analysis.py'); run('src/export_gurobi.py')
if not a.sem_figuras: run('src/figures.py')
if a.gurobi: run('verificar_gurobi.py', cwd=ROOT / 'verificacao')
print('Concluído. Os resultados estão em out/results.json')
