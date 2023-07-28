# pydelaunay
Python binding of Parallel Delaunay Triangulation from **[ParGeo](https://github.com/ParAlg/ParGeo)** lib       
ParGeo is heavily based on **[parlaylib](https://cmuparlay.github.io/parlaylib/)**, so you'll need to install it as in [this guide](https://cmuparlay.github.io/parlaylib/installation.html)      
### Installation       
```console
cd pydelaunay
python3 setup.py develop
```
For the last command you may need to manually write install directory which is visible by python to --install-dir arg, or just installation with sudo     
### Usage      
```python
import numpy as np
from pydelaunay import delaunay
points = np.random.binomial(1024, 0.5, (1000, 2))
triangulation = delaunay(points)
```
