# -*- coding: utf-8 -*-
"""Exemplo numérico simples do artigo (Tabela 1): cinco rotas e dois critérios.

As contas cabem à mão e usam os mesmos conceitos do modelo: risco-hora, dominância, solução eficiente,
ponto ideal e nadir, normalização, soma ponderada, solução não suportada e ε-restrição.
Os números são inventados para ilustrar e não vêm dos dados do artigo. O script só usa a biblioteca padrão.

Uso:
    python src/exemplo_didatico.py
"""
import sys
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')

# A rota real tem dois trechos: (duração em horas, risco do trecho entre 0 e 1).
REAL_TRECHOS = [(60.0, 0.30), (42.0, 0.46)]
# Em cada rota: (tempo total F_T em horas, risco-hora total F_H). O risco-hora de um trecho é risco x duração.
ROTAS = {
    'Real': (sum(t for t, _ in REAL_TRECHOS), sum(t * r for t, r in REAL_TRECHOS)),
    'A': (96.0, 34.0),
    'B': (98.4, 32.4),
    'C': (100.0, 30.0),
    'D': (99.0, 33.4),
}


def f(x, nd=1):
    """número com vírgula decimal, como nas tabelas do artigo"""
    return f'{x:.{nd}f}'.replace('.', ',')


def domina(p, q):
    """p domina q: p não é pior em nenhum critério e é melhor em pelo menos um (minimização)"""
    return p[0] <= q[0] and p[1] <= q[1] and (p[0] < q[0] or p[1] < q[1])


def analisar(rotas=ROTAS):
    r = {}
    r['dominada_por'] = {a: [b for b in rotas if b != a and domina(rotas[b], rotas[a])] for a in rotas}
    r['eficientes'] = [a for a in rotas if not r['dominada_por'][a]]
    ef = {a: rotas[a] for a in r['eficientes']}
    # tabela de payoff por ordem lexicográfica: tempo mínimo (desempate: risco-hora) e risco-hora mínimo (desempate: tempo)
    a_t = min(ef, key=lambda a: (ef[a][0], ef[a][1]))
    a_h = min(ef, key=lambda a: (ef[a][1], ef[a][0]))
    r['rota_tempo_min'], r['rota_risco_min'] = a_t, a_h
    r['ideal'] = (ef[a_t][0], ef[a_h][1])
    r['nadir'] = (ef[a_h][0], ef[a_t][1])
    (tI, hI), (tN, hN) = r['ideal'], r['nadir']
    r['faixa'] = (tN - tI, hN - hI)
    r['norm'] = {a: ((rotas[a][0] - tI) / (tN - tI), (rotas[a][1] - hI) / (hN - hI)) for a in rotas}
    # soluções suportadas: vértices da envoltória convexa inferior dos pontos eficientes (varredura de pesos)
    venc = set()
    for k in range(1, 1000):
        wT = k / 1000
        z = {a: wT * r['norm'][a][0] + (1 - wT) * r['norm'][a][1] for a in rotas}
        m = min(z.values())
        venc |= {a for a in rotas if abs(z[a] - m) < 1e-12}
    r['suportadas'] = sorted(a for a in r['eficientes'] if a in venc)
    r['nao_suportadas'] = sorted(a for a in r['eficientes'] if a not in venc)
    return r


def soma_ponderada(r, wT, wH, rotas=ROTAS):
    z = {a: wT * r['norm'][a][0] + wH * r['norm'][a][1] for a in rotas}
    return z, min(z, key=z.get)


def epsilon(eps, rotas=ROTAS):
    """menor tempo entre as rotas com risco-hora <= eps"""
    ok = {a: v for a, v in rotas.items() if v[1] <= eps}
    return min(ok, key=lambda a: ok[a][0]) if ok else None


def rho(r, wT, wH):
    """taxa de substituição implícita (Eq. 11 do artigo): risco-hora que vale uma hora de viagem"""
    return (wT / r['faixa'][0]) / (wH / r['faixa'][1])


def tabela_exemplo():
    r = analisar()
    cab = ['Rota', 'Tempo (h)', 'Risco-hora', 'R méd.', 'Tempo normalizado', 'Risco-hora normalizado', 'Situação']
    L = []
    for a, (T, H) in ROTAS.items():
        if a in r['suportadas']:
            sit = 'Eficiente, suportada'
        elif a in r['nao_suportadas']:
            sit = 'Eficiente, não suportada'
        else:
            sit = 'Dominada'
        L.append(['Real' if a == 'Real' else f'Rota {a}', f(T), f(H, 2), f(H / T, 3),
                  f(r['norm'][a][0], 2), f(r['norm'][a][1], 2), sit])
    return cab, L


def verificar(r):
    """confere os números citados no artigo (Seção 4.6); falha com AssertionError se algum mudar"""
    assert r['eficientes'] == ['A', 'B', 'C'], r['eficientes']
    assert r['suportadas'] == ['A', 'C'] and r['nao_suportadas'] == ['B']
    assert r['ideal'] == (96.0, 30.0) and r['nadir'] == (100.0, 34.0)
    assert set(r['dominada_por']['Real']) == {'A', 'B', 'C', 'D'} and r['dominada_por']['D'] == ['B']
    z1, a1 = soma_ponderada(r, 0.6, 0.4); z2, a2 = soma_ponderada(r, 0.4, 0.6)
    assert a1 == 'A' and a2 == 'C'
    assert [round(z1[a], 2) for a in 'ABC'] == [0.40, 0.60, 0.60] and [round(z2[a], 2) for a in 'ABC'] == [0.60, 0.60, 0.40]
    assert epsilon(33.0) == 'B'
    assert abs(rho(r, 0.6, 0.4) - 1.5) < 1e-12 and abs(rho(r, 0.4, 0.6) - 2 / 3) < 1e-12
    return True


def main():
    r = analisar()
    print('Rota real, por trecho: duração (h) x risco = risco-hora')
    for t, k in REAL_TRECHOS:
        print(f'   {f(t, 0)} h x {f(k, 2)} = {f(t * k, 2)}')
    T, H = ROTAS['Real']
    print(f'   total: F_T = {f(T)} h, F_H = {f(H, 2)}, risco médio = F_H / F_T = {f(H / T, 3)}\n')
    cab, L = tabela_exemplo()
    larg = [max(len(str(x[j])) for x in [cab] + L) for j in range(len(cab))]
    print('--- Tabela 1 ---')
    for x in [cab] + L:
        print('  '.join(str(c).ljust(larg[j]) for j, c in enumerate(x)).rstrip())
    print()
    print('Dominância (quem domina quem):')
    for a, d in r['dominada_por'].items():
        print(f'   rota {a}: ' + ('eficiente (nenhuma rota a domina)' if not d else 'dominada por ' + ', '.join(d)))
    print(f'\nTabela de payoff: tempo mínimo = rota {r["rota_tempo_min"]}; risco-hora mínimo = rota {r["rota_risco_min"]}')
    print(f'   ideal = ({f(r["ideal"][0])} h; {f(r["ideal"][1])}),  nadir = ({f(r["nadir"][0])} h; {f(r["nadir"][1])}),  '
          f'faixas = ({f(r["faixa"][0])} h; {f(r["faixa"][1])})')
    for wT, wH in ((0.6, 0.4), (0.5, 0.5), (0.4, 0.6)):
        z, a = soma_ponderada(r, wT, wH)
        txt = '  '.join(f'{x}={f(z[x], 2)}' for x in ROTAS)
        print(f'\nSoma ponderada, pesos ({f(wT, 1)}; {f(wH, 1)}): Z = {txt}')
        print(f'   menor Z: rota {a}' + ('  (empate entre A e C)' if abs(z["A"] - z["C"]) < 1e-12 else '') +
              f';  taxa de substituição ρ = {f(rho(r, wT, wH), 2)} risco-hora por hora')
    print(f'\nVarrendo os pesos de 0,001 a 0,999: só as rotas {", ".join(r["suportadas"])} vencem alguma vez '
          f'(suportadas); a rota {", ".join(r["nao_suportadas"])}, eficiente, nunca vence (não suportada).')
    print('   Troca de A para B: +2,4 h e -1,6 de risco-hora (0,67 por hora); de B para C: +1,6 h e -2,4 (1,5 por hora).')
    print('   A taxa de troca cresce ao longo da fronteira; por isso nenhum peso escolhe B.')
    print(f'\nε-restrição: menor tempo com risco-hora <= 33,0  ->  rota {epsilon(33.0)}')
    verificar(r)
    print('\nConferência dos números citados no artigo: OK')


if __name__ == '__main__':
    main()
