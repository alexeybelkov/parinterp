# pydelaunay
Python binding of Parallel [Delaunay Triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) from **[ParGeo](https://github.com/ParAlg/ParGeo)** lib       
ParGeo is heavily based on **[parlaylib](https://cmuparlay.github.io/parlaylib/)**, so you'll need to install it as in [this guide](https://cmuparlay.github.io/parlaylib/installation.html) and add it to CMakeLists.txt      
### Installation       
You'll have to manually put all required pathes into CMakeLists.txt       
```console
git clone https://github.com/alexeybelkov/pydelaunay.git
cd pydelaunay
python3 setup.py develop
```
For the last command you may need to manually write install directory which is visible by python to --install-dir arg, or just run installation with sudo     
### Usage      
```python
import numpy as np
from pydelaunay import delaunay
points = np.random.randint(low=0, high=1024, size=(512, 2))
points = np.unique(points, axis=0)
triangulation = delaunay(points)
```
### Benchmarking       
![](/benchmark/benchmark.jpg)
