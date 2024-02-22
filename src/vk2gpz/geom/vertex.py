from abc import ABCMeta

import numpy as np
from numpy import ndarray

class IManifold:
    pass

class Vertex(metaclass=ABCMeta):

    def __init__(self, x: int, y: int):
        self.visited: bool = False
        self.x: int = x
        self.y: int = y
        self.manifold: IManifold
        self.data = None
        self.color = None
        self.id: int = -1
        self.coord: ndarray = np.array([0.0, 0.0, 0.0])

    def get_neighbours_in_distance(self, src_vertex, distance):
        if self.manifold is not None:
            return self.manifold.get_neighbours_in_distance(src_vertex, distance)
        return None

    def set_data(self, data):
        self.data = data