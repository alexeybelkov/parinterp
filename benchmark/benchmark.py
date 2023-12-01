from parinterp import Linear2DInterpolator
from scipy.interpolate import griddata
import numpy as np
from typing import Mapping
from deli import load_json
from timeit import timeit
from tqdm import tqdm
from matplotlib import pyplot as plt


path = load_json('benchmark/path.json')

arrays: Mapping[str, np.ndarray]= dict()

for key, p in path.items():
    arrays[key] = np.load(p, allow_pickle=True)

sizes = []
times = {'parinterp': [], 'scipy': []}
max_error = []
q = 0.99
q_error = []
for k, x in tqdm(arrays.items()):
    add_cols = x.shape[1] // 4
    distances = np.concatenate((x[...,-add_cols:], x, x[...,:add_cols] ), axis=1)
    int_points = np.transpose((np.isnan(distances)).nonzero())
    x_points: np.ndarray = np.transpose((~np.isnan(distances)).nonzero())
    x_values = distances[~np.isnan(distances)]
    surface_map = distances.copy()

    my_interp = lambda: Linear2DInterpolator(x_points, -1)(int_points, x_values, 0.0)
    sp_interp = lambda: griddata(x_points, x_values, int_points, method='linear', fill_value=0.0)

    sizes.append(x_points.size + int_points.size)

    times['parinterp'].append(timeit(my_interp, number=4) / 4)
    times['scipy'].append(timeit(sp_interp, number=4) / 4)

    my_interp = my_interp()
    sp_interp = sp_interp()

    ae = np.abs(my_interp - sp_interp)
    max_error.append(ae.max())
    q_error.append(np.quantile(ae, q))

indexes = np.argsort(sizes)
sizes = [sizes[i] for i in indexes]
times['parinterp'] = [times['parinterp'][i] for i in indexes]
times['scipy'] = [times['scipy'][i] for i in indexes]
max_error = [max_error[i] for i in indexes]
q_error = [q_error[i] for i in indexes]

print('array sizes:', sizes)
print('max error:', max_error)
print('q99 error:', q_error)
plt.plot(sizes, max_error, label='max_error', marker='x')
plt.plot(sizes, q_error, label=f'q{q}_error', marker='x')

plt.legend()
plt.title('Interpolation error')
plt.xlabel('array size')
plt.ylabel('abs error')
plt.savefig('benchmark/interpolation_error.jpg')
# plt.show()
plt.clf()

plt.plot(sizes, times['parinterp'], label='parinterp', marker='x')
plt.plot(sizes, times['scipy'], label='scipy.interpolate.griddata', marker='x')
plt.legend()
plt.title('Interpolation\n' + '\n'.join([f'{k} :: <t>: {np.round(np.mean(times[k]),3)} | t_min: {np.round(np.min(times[k]), 3)} | t_max: {np.round(np.max(times[k]), 3)}' for k in times.keys()]))
plt.xlabel('array size')
plt.ylabel('seconds')
plt.savefig('benchmark/interpolation_time.jpg')
# plt.show()

plt.clf()
