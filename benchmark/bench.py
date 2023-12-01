from parinterp import Linear2DInterpolator
from scipy.interpolate import griddata
import numpy as np
from typing import Mapping
from deli import load_json


path = load_json('benchmark/path.json')

arrays: Mapping[str, np.ndarray]= dict()

for key, p in path.items():
    arrays[key] = np.load(p, allow_pickle=True)

for k, x in tqdm(arrays.items()):
    add_cols = x.shape[1] // 4
    distances = np.concatenate((x[...,-add_cols:], x, x[...,:add_cols] ), axis=1)
    int_points = np.transpose((np.isnan(distances)).nonzero())
    x_points = np.transpose((~np.isnan(distances)).nonzero())
    x_values = distances[~np.isnan(distances)]
    surface_map = distances.copy()

    my_interp = lambda: Linear2DInterpolator(x_points, 4)(int_points, x_values, 0.0)
    sp_interp = lambda: griddata(x_points, x_values, int_points, method='linear', fill_value=0.0)

    sizes.append(len(int_points))

    times['pydelaunay'].append(timeit(my_interp, number=10) / 10)
    times['scipy'].append(timeit(sp_interp, number=10) / 10)

    my_interp = my_interp()
    sp_interp = sp_interp()

    mape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.5))
    qape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.75))
    tape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.95))


