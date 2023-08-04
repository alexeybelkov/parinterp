import pydelaunay
import numpy as np
from timeit import timeit
from typing import Mapping
from deli import load_json
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay
from tqdm import tqdm
from interp.interpolator import Linear2DInterpolator
from scipy.interpolate import griddata
import os

path = load_json('benchmark/path.json')

arrays: Mapping[str, np.ndarray]= dict()

for key, p in path.items():
    arrays[key] = np.load(p, allow_pickle=True)

times = {'pydelaunay': [], 'scipy': []}
sizes = []

os.environ['PARLAY_NUM_THREADS'] = '32'

for k, x in tqdm(arrays.items()):
    points = np.transpose((~np.isnan(x)).nonzero())

    my_delaunay = lambda: pydelaunay.delaunay(points)
    sp_delaunay = lambda: Delaunay(points)

    sizes.append(len(points))

    times['pydelaunay'].append(timeit(my_delaunay, number=10) / 10)
    times['scipy'].append(timeit(sp_delaunay, number=10) / 10)

indexes = np.argsort(sizes)
sizes = [sizes[i] for i in indexes]
times['pydelaunay'] = [times['pydelaunay'][i] for i in indexes]
times['scipy'] = [times['scipy'][i] for i in indexes]

plt.plot(sizes, times['pydelaunay'], label='pydelaunay', marker='x')
plt.plot(sizes, times['scipy'], label='scipy.spatial.Delaunay', marker='x')
plt.legend()
plt.title('Triangulation\n' + '\n'.join([f'{k} :: <t>: {np.round(np.mean(times[k]),3)} | t_min: {np.round(np.min(times[k]), 3)} | t_max: {np.round(np.max(times[k]), 3)}' for k in times.keys()]))
plt.xlabel('array size')
plt.ylabel('seconds')
plt.savefig('benchmark/triangulation.jpg')
plt.clf()

times = {'pydelaunay': [], 'scipy': []}
sizes = []
mape = []
qape = []
tape = []

for k, x in tqdm(arrays.items()):
    add_cols = x.shape[1] // 4
    distances = np.concatenate((x[...,-add_cols:], x, x[...,:add_cols] ), axis=1)
    int_points = np.transpose((np.isnan(distances)).nonzero())
    x_points = np.transpose((~np.isnan(distances)).nonzero())
    x_values = distances[~np.isnan(distances)]
    surface_map = distances.copy()

    my_interp = lambda: Linear2DInterpolator(x_points, x_values)(int_points)
    sp_interp = lambda: griddata(x_points, x_values, int_points, method='linear', fill_value=0.0)

    sizes.append(len(int_points))

    times['pydelaunay'].append(timeit(my_interp, number=10) / 10)
    times['scipy'].append(timeit(sp_interp, number=10) / 10)

    my_interp = my_interp()
    sp_interp = sp_interp()

    mape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.5))
    qape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.75))
    tape.append(np.quantile(np.abs(my_interp - sp_interp) / (sp_interp + 1e-8), 0.95))

indexes = np.argsort(sizes)
sizes = [sizes[i] for i in indexes]
times['pydelaunay'] = [times['pydelaunay'][i] for i in indexes]
times['scipy'] = [times['scipy'][i] for i in indexes]
mape = [mape[i] for i in indexes]
qape = [qape[i] for i in indexes]
tape = [tape[i] for i in indexes]

plt.plot(sizes, times['pydelaunay'], label='pydelaunay', marker='x')
plt.plot(sizes, times['scipy'], label='scipy.interpolate.griddata', marker='x')
plt.legend()
plt.title('Interpolation\n' + '\n'.join([f'{k} :: <t>: {np.round(np.mean(times[k]),3)} | t_min: {np.round(np.min(times[k]), 3)} | t_max: {np.round(np.max(times[k]), 3)}' for k in times.keys()]))
plt.xlabel('array size')
plt.ylabel('seconds')
plt.savefig('benchmark/interpolation.jpg')
plt.clf()

plt.plot(sizes, mape, label='median ape', marker='x')
plt.plot(sizes, qape, label='q0.75 ape', marker='x')
plt.plot(sizes, tape, label='q0.95 ape', marker='x')
plt.legend()
plt.title('Interpolation APE')
plt.xlabel('array size')
plt.ylabel('absolute percentage error')
plt.savefig('benchmark/interpolation_error.jpg')
plt.clf()
