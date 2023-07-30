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

times = {'my': [], 'sp': []}
sizes = []

for k, x in tqdm(arrays.items()):
    points = np.transpose((~np.isnan(x)).nonzero())

    my_delaunay = lambda: pydelaunay.delaunay(points)
    sp_delaunay = lambda: Delaunay(points)

    sizes.append(len(points))

    times['my'].append(timeit(my_delaunay, number=10) / 10)
    times['sp'].append(timeit(sp_delaunay, number=10) / 10)

indexes = np.argsort(sizes)
sizes = [sizes[i] for i in indexes]
times['my'] = [times['my'][i] for i in indexes]
times['sp'] = [times['sp'][i] for i in indexes]

plt.plot(sizes, times['my'], label='pydelaunay', marker='x')
plt.plot(sizes, times['sp'], label='scipy.spatial.Delaunay', marker='x')
plt.legend()
plt.xlabel('array size')
plt.ylabel('seconds')
plt.savefig('benchmark/benchmark.jpg')