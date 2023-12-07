# %%
from deli import load_numpy, save_numpy
from parinterp import Linear2DInterpolator
from scipy.interpolate import griddata
import numpy as np
import os
from tqdm.notebook import tqdm
from matplotlib import pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--n-jobs", required=True, help="Set number of workers", type=int)
args = parser.parse_args()

n_jobs = args.n_jobs

arrays = {}
for npy_name in tqdm(list(os.walk('/shared/personal/leonmeon/griddata'))[0][-1]):
    arrays[npy_name[:-4]] = load_numpy('/shared/personal/leonmeon/griddata/' + npy_name)

q99_err = []
max_err = []
parinterp_time = []
scipy_time = []
sizes = []

for k, x in tqdm(arrays.items()):
    add_cols = x.shape[1] // 4
    distances = np.concatenate((x[...,-add_cols:], x, x[...,:add_cols] ), axis=1)
    int_points = np.transpose((np.isnan(distances)).nonzero())
    x_points: np.ndarray = np.transpose((~np.isnan(distances)).nonzero())
    x_values = distances[~np.isnan(distances)]
    parinterp_map = distances.copy()
    scipy_map = distances.copy()
    parinterp_t = time.time()
    parinterp_values = Linear2DInterpolator(x_points, n_jobs)(int_points, x_values, fill_value=0.0)
    parinterp_t = time.time() - parinterp_t
    scipy_t = time.time()
    scipy_values = griddata(x_points, x_values, int_points, method='linear', fill_value=0.0)
    scipy_t = time.time() - scipy_t

    parinterp_time.append(parinterp_t)
    scipy_time.append(scipy_t)

    sizes.append(len(x_points) + len(int_points))

    for i, v in enumerate(parinterp_values):
        parinterp_map[int_points[i][0], int_points[i][1]] = v

    for i, v in enumerate(scipy_values):
        scipy_map[int_points[i][0], int_points[i][1]] = v


    delta = np.abs(parinterp_map - scipy_map)
    q99_delta = np.quantile(delta, 0.99)
    max_delta = np.max(delta)
    outliers_num = np.sum(delta > q99_delta)
    q99_err.append(q99_delta)
    max_err.append(max_delta)

indexes = np.argsort(sizes)
sizes.sort()
parinterp_time = [parinterp_time[i] for i in indexes]
scipy_time = [scipy_time[i] for i in indexes]
plt.plot(sizes, parinterp_time, marker='x', label='parinterp')
plt.plot(sizes, scipy_time, marker='x', label='scipy griddata')
plt.title('Exectution time')
plt.show()
plt.hist(q99_err, bins=32)
plt.title('q99 of delta')
plt.show()

plt.hist(max_err, bins=32)
plt.title('max of delta')
plt.show()


