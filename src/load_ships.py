import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import *
from wr_model import *
SOG = json.load(open(DATA / 'parametros.json', encoding='utf-8'))


def load(k, wp='idw', use_current=True, mesh=True):
    S = Ship(k, SOG[k], MESHDIR / f'{k}_nos.csv', MESHDIR / f'{k}_arcos.csv', wp_treatment=wp)
    S.calibrate(SOG[k]['v_obs'], use_current); S.build_arcs()
    return S
