import numpy as np
from pydelaunay import BiLinearInterpolator
from scipy.spatial import KDTree
from typing import Union

class Linear2DInterpolator(BiLinearInterpolator):
    def __init__(self, points: np.ndarray[Union[np.int32, np.double]], values: np.ndarray[np.float32], **kwargs):
        print('starting triangulation ....')
        super().__init__(points.astype(np.double), values)
        print('fininshed triangulation')
        print('starting tree building')
        self.kdtree = KDTree(data=self.triangulation.vertices, **kwargs)
        print('finished tree building')

    def __call__(self, points: np.ndarray[Union[np.int32, np.double]], fill_value: float = 0.0, workers: int = -1) -> np.ndarray[np.float32]:
        print('starting query')
        _, neighbors = self.kdtree.query(points, 1, workers=workers)
        print('finished query')
        print('starting interpolation')
        int_values: np.ndarray[np.float32] = super().__call__(points.astype(np.double), neighbors, fill_value)
        print('finished interpolation')
        return int_values.astype(np.float32)