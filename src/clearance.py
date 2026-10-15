import numpy as np
from global_land_mask import globe
from scipy import ndimage

class ClearGrid:
    """Grade de distância (nm) até a terra, por faixas de latitude (anisotropia local), a partir da máscara GLOBE (1 km)."""
    def __init__(self, lat0=-27.0, lat1=5.0, lon0=-52.0, lon1=-33.0, res=0.01, band=4.0):
        self.res = res; self.bands = []
        lat = lat0
        while lat < lat1 - 1e-9:
            la1 = min(lat + band, lat1)
            pad = 1.5                                   # sobreposição para evitar efeito de borda
            gl0, gl1 = max(lat - pad, lat0), min(la1 + pad, lat1)
            lats = np.arange(gl0, gl1 + 1e-9, res); lons = np.arange(lon0, lon1 + 1e-9, res)
            LA, LO = np.meshgrid(lats, lons, indexing='ij')
            land = globe.is_land(LA, LO)
            mid = 0.5 * (gl0 + gl1)
            d = ndimage.distance_transform_edt(~land, sampling=(res * 60.0, res * 60.0 * np.cos(np.radians(mid))))
            self.bands.append((lat, la1, gl0, lats, lons, d.astype(np.float32), land))
            lat = la1
        self.lon0, self.lon1 = lon0, lon1

    def _band(self, lat):
        for b in self.bands:
            if b[0] <= lat <= b[1] + 1e-9:
                return b
        return self.bands[-1]

    def clearance(self, lat, lon):
        lat = np.atleast_1d(lat).astype(float); lon = np.atleast_1d(lon).astype(float)
        out = np.empty(lat.shape)
        for k in range(len(lat)):
            b = self._band(lat[k]); _, _, gl0, lats, lons, d, land = b
            i = int(round((lat[k] - lats[0]) / self.res)); j = int(round((lon[k] - lons[0]) / self.res))
            i = min(max(i, 0), d.shape[0] - 1); j = min(max(j, 0), d.shape[1] - 1)
            out[k] = d[i, j]
        return out

    def is_land_grid(self, lat, lon):
        return globe.is_land(lat, lon)
