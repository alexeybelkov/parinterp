import numpy as np
from pydelaunay import BiLinearInterpolator
from scipy.spatial import KDTree


class Linear2DInterpolator(BiLinearInterpolator):
    def __init__(self, points: np.ndarray[np.int32], values: np.ndarray[np.float32]):
        super().__init__(points.astype(np.double), values)
        self.kdtree = KDTree(points)

    def __call__(self, points: np.ndarray[np.int32], fill_value: float = 0.0, workers: int = -1) -> np.ndarray[np.float32]:
        _, neighbors = self.kdtree.query(points, 1, workers=workers)
        int_values: np.ndarray[np.float32] = super().__call__(points.astype(np.double), neighbors, fill_value)
        return int_values.astype(np.float32)