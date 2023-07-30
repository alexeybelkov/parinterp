import pydelaunay
import numpy as np
from timeit import timeit
from typing import Mapping
from deli import load_json
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay
from tqdm import tqdm

path = load_json('benchmark/path.json')

arrays: Mapping[str, np.ndarray]= dict()

for key, p in path.items():
    arrays[key] = np.load(p, allow_pickle=True)

times = {'pydelaunay': [], 'scipy': []}
sizes = []

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
plt.title('\n'.join([f'{k} :: <t>: {np.round(np.mean(times[k]),3)} | t_min: {np.round(np.min(times[k]), 3)} | t_max: {np.round(np.max(times[k]), 3)}' for k in times.keys()]))
plt.xlabel('array size')
plt.ylabel('seconds')
plt.savefig('benchmark/benchmark.jpg')