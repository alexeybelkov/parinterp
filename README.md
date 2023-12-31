# parinterp

Python library for parallel linear 2D interpolation, it consists of two main parts:
- [Delaunay Triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) from [delaunator-cpp](https://github.com/abellgithub/delaunator-cpp) lib, still sequential part        
- Linear 2D interpolation via barycentric coordinates, parallel part

Being almost fully written in C++ and binded via [pybind11](https://pybind11.readthedocs.io/en/stable/index.html), library achieved very high perfomance

### Installation             
```console
git clone https://github.com/alexeybelkov/parinterp
cd parinterp
pip install -r requirements.txt
./install.sh
```
### Usage      
```python
import numpy as np
from parinterp import Linear2DInterpolator
n, m = 1024, 2
points = np.random.randint(low=0, high=1024, size=(n, m))
points = np.unique(points, axis=0)
x_points = points[: n // 2]
values = np.random.uniform(low=0.0, high=1.0, size=(len(x_points),))
interp_points = points[n // 2:]
n_jobs = -1 # will be equal to num of CPU cores
interpolator = Linear2DInterpolator(x_points, values, n_jobs)
# Also you can pass values to __call__ and rewrite the ones that were passed to __init__
interp_values = interpolator(interp_points, values + 1.0, fill_value=0.0)
```

### TODO
* Currently, there are a lot of presumably unnecessary reallocations, try to find a way to remove them        
* Using Python's scipy KDTree is more like a crutch, find fast C++ one-nearest-neighbour algorithm implementation         
* Find a way to build a fast bijection $f : S \rightarrow \{0, 1, ..., n\}$, where $S \subset \mathbf{N}$       
* Maybe consider this [Parallel Nearest Neighbors in Low Dimensions with Batch Updates](https://arxiv.org/pdf/2111.04182.pdf)      
* Compare perfomance with and without safety

### Usefull links
Point location algorithm was highly inspired by paper [Fast randomized point location without preprocessing in two- and
three-dimensional Delaunay triangulations](https://web.cs.ucdavis.edu/~amenta/w07/jump.and.walk.pdf)
<!-- ### Benchmarking
#### Triangulation
![](/benchmark/triangulation.jpg)
#### Interpolation
![](/benchmark/interpolation.jpg)
#### Interpolation error
![](/benchmark/interpolation_error.jpg) -->
