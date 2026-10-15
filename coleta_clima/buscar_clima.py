# -*- coding: utf-8 -*-
"""Busca vento, onda e corrente históricos (Open-Meteo, gratuito, sem cadastro) para uma lista de pontos e horas.

Entrada: um CSV com as colunas  Latitude, Longitude, Sample time  (hora em UTC, no formato AAAA-MM-DD HH:MM:SS+00:00).
Saída:   o mesmo CSV com 6 colunas a mais:
         vento_vel_kmh, vento_dir_graus, onda_altura_m, onda_dir_graus, corrente_vel_kmh, corrente_dir_graus

Uso (Windows):   py buscar_clima.py entrada.csv saida.csv
Antes:           py -m pip install pandas requests

Fontes (a Open-Meteo escolhe o modelo, "best match"):
  vento     https://archive-api.open-meteo.com/v1/archive   (wind_speed_10m, wind_direction_10m; km/h e graus)
  onda      https://marine-api.open-meteo.com/v1/marine     (wave_height, wave_direction; m e graus)
  corrente  https://marine-api.open-meteo.com/v1/marine     (ocean_current_velocity, ocean_current_direction; km/h e graus
            -- direção PARA onde a água flui, convenção oceanográfica)

Este arquivo reúne, sem mudar a lógica, os dois scripts que foram usados para obter o clima da rota real e o clima dos nós
da malha. Ele não foi executado de novo na revisão do artigo, porque o ambiente de revisão não alcança a API.
A Open-Meteo limita o número de chamadas por minuto e por dia; por isso há uma pausa entre as chamadas.
"""
import sys, time
# No Windows, a saída redirecionada para um arquivo usa a codificação cp1252 e quebra com símbolos como ε e →.
# Nesse caso, força UTF-8 (num terminal normal nada muda).
for _f in (sys.stdout, sys.stderr):
    if _f is not None and not _f.isatty() and hasattr(_f, 'reconfigure'):
        _f.reconfigure(encoding='utf-8')
import pandas as pd
import requests

ARCHIVE_URL = 'https://archive-api.open-meteo.com/v1/archive'
MARINE_URL = 'https://marine-api.open-meteo.com/v1/marine'
PAUSA_S = 0.3


def pedir(url, lat, lon, dia, variaveis):
    r = requests.get(url, params={'latitude': lat, 'longitude': lon, 'start_date': dia, 'end_date': dia,
                                  'hourly': variaveis, 'timezone': 'UTC'}, timeout=15)
    r.raise_for_status()
    return r.json()['hourly']


def buscar(lat, lon, instante):
    dia, hora = instante.strftime('%Y-%m-%d'), instante.strftime('%Y-%m-%dT%H:%M')
    out = dict(vento_vel_kmh=None, vento_dir_graus=None, onda_altura_m=None, onda_dir_graus=None,
               corrente_vel_kmh=None, corrente_dir_graus=None)
    try:
        h = pedir(ARCHIVE_URL, lat, lon, dia, 'wind_speed_10m,wind_direction_10m')
        if hora in h['time']:
            i = h['time'].index(hora)
            out.update(vento_vel_kmh=h['wind_speed_10m'][i], vento_dir_graus=h['wind_direction_10m'][i])
    except Exception as e:
        print(f'  [vento] falhou em {lat},{lon} {dia}: {e}')
    time.sleep(PAUSA_S)
    try:
        h = pedir(MARINE_URL, lat, lon, dia, 'wave_height,wave_direction,ocean_current_velocity,ocean_current_direction')
        if hora in h['time']:
            i = h['time'].index(hora)
            out.update(onda_altura_m=h['wave_height'][i], onda_dir_graus=h['wave_direction'][i],
                       corrente_vel_kmh=h['ocean_current_velocity'][i], corrente_dir_graus=h['ocean_current_direction'][i])
    except Exception as e:
        print(f'  [mar] falhou em {lat},{lon} {dia}: {e}')
    time.sleep(PAUSA_S)
    return out


def main(entrada, saida):
    d = pd.read_csv(entrada)
    d['Sample time'] = pd.to_datetime(d['Sample time'])
    linhas = []
    for i, r in d.iterrows():
        linhas.append({**r.to_dict(), **buscar(r['Latitude'], r['Longitude'], r['Sample time'])})
        if i % 10 == 0:
            print(f'  {i}/{len(d)} pontos processados...')
    pd.DataFrame(linhas).to_csv(saida, index=False)
    print(f'Salvo: {saida} ({len(linhas)} linhas)')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Uso: py buscar_clima.py entrada.csv saida.csv')
    main(sys.argv[1], sys.argv[2])
