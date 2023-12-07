import numpy as np
from cpp.build.parinterp import Linear2DInterpolatorCpp
from scipy.spatial import KDTree, Delaunay

class Linear2DInterpolator(Linear2DInterpolatorCpp):
    def __init__(self, points: np.ndarray, n_jobs: int = 1, **kwargs):
        simplices = Delaunay(points).simplices
        super().__init__(points, simplices, n_jobs)
        self.kdtree = KDTree(data=points, **kwargs)
        self.n_jobs = n_jobs

    def __call__(self, points: np.ndarray, values: np.ndarray, fill_value: float = 0.0) -> np.ndarray[np.float32]:
        _, neighbors = self.kdtree.query(points, 1, workers=self.n_jobs)
        int_values = super().__call__(points, values, neighbors, fill_value)
        return int_values
    
    def naive_points_location(self, points: np.ndarray):
        _, neighbors = self.kdtree.query(points, 1, workers=self.n_jobs)
        return super().naive_points_location(points, neighbors)
    