import numpy as np
from global_land_mask import globe
from wr_model import haversine

def sample_segment(la1, lo1, la2, lo2, step_nm=1.0):
    d = float(haversine(la1, lo1, la2, lo2)); n = max(int(np.ceil(d/step_nm)), 1)
    f = np.linspace(0, 1, n+1)
    return la1 + f*(la2-la1), lo1 + f*(lo2-lo1), f*d

def path_samples(nodes, path, step_nm=1.0):
    LA, LO, DD = [], [], []; off = 0.0
    for i, j in zip(path[:-1], path[1:]):
        la, lo, dd = sample_segment(nodes.at[i,'Latitude'], nodes.at[i,'Longitude'], nodes.at[j,'Latitude'], nodes.at[j,'Longitude'], step_nm)
        LA.append(la); LO.append(lo); DD.append(dd + off); off += float(dd[-1])
    return np.concatenate(LA), np.concatenate(LO), np.concatenate(DD)

def clearance_nm(lat, lon, rmax=6.0, dr=0.25, ndir=48):
    """menor distância (nm) até um pixel de terra da máscara GLOBE (1 km), por anéis; rmax se não achar"""
    lat = np.asarray(lat); lon = np.asarray(lon)
    clr = np.full(lat.shape, rmax)
    ang = np.linspace(0, 2*np.pi, ndir, endpoint=False)
    for r in np.arange(dr, rmax + 1e-9, dr):
        dlat = (r/60.0) * np.cos(ang)[None, :]
        dlon = (r/60.0) * np.sin(ang)[None, :] / np.cos(np.radians(lat))[:, None]
        hit = globe.is_land(np.clip(lat[:, None] + dlat, -89.9, 89.9), lon[:, None] + dlon).any(axis=1)
        clr = np.where(hit & (clr == rmax), r, clr)
    return clr

def route_geo_report(nodes, path, port_margin_nm=15.0, step_nm=1.0, dmin=3.0):
    la, lo, dd = path_samples(nodes, path, step_nm)
    tot = dd[-1]
    land = globe.is_land(la, lo)
    core = (dd >= port_margin_nm) & (dd <= tot - port_margin_nm)          # exclui aproximações portuárias
    clr = clearance_nm(la[core], lo[core], rmax=6.0)
    return dict(length=float(tot), n_samples=int(len(la)), land_samples=int(land.sum()), land_samples_core=int(land[core].sum()),
                min_clearance=float(clr.min()), frac_ge_dmin=float((clr >= dmin).mean()), frac_lt_1=float((clr < 1.0).mean()))
