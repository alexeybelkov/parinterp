import pydelaunay
import numpy as np
from timeit import timeit, repeat
from typing import Mapping, List
from deli import load_json
from matplotlib import pyplot as plt
from tqdm import tqdm
from pydelaunay_interp.interpolator import Linear2DInterpolator
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--n-jobs", required=True, help="Set number of workers", type=int)
args = parser.parse_args()


MAX_WORKERS = args.n_jobs
WORKERS_LIST = [1] + list(range(2, MAX_WORKERS + 1, 2))

path = load_json('benchmark/path.json')

arrays: Mapping[str, List[np.ndarray]] = dict()

for key, p in path.items():
    x = np.load(p, allow_pickle=True)
    add_cols = x.shape[1] // 4
    distances = np.concatenate((x[...,-add_cols:], x, x[...,:add_cols] ), axis=1)
    int_points = np.transpose((np.isnan(distances)).nonzero())
    x_points = np.transpose((~np.isnan(distances)).nonzero())
    x_values = distances[~np.isnan(distances)]
    arrays[key] = [x_points, x_values, int_points]

for k, arrs in tqdm(arrays.items()):
    x_points, x_values, int_points = arrs
    arr_size = len(x_points) + len(int_points)
    times = []

    for workers in WORKERS_LIST:
        os.environ['PARLAY_NUM_THREADS'] = str(workers)
        interp = lambda: Linear2DInterpolator(x_points, x_values)(int_points, fill_value=0.0, workers=1)
        times.append(repeat(interp, number=1, repeat=16))

    plt.boxplot(times, labels=WORKERS_LIST, vert=True)
    plt.title(f'Array size: {len(x_points) + len(int_points)}')
    plt.xlabel('num workers')
    plt.ylabel('')
    plt.savefig(f'benchmark/1kd_{arr_size}_num_workers.jpg')
    plt.clf()
