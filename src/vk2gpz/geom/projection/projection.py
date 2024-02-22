from abc import ABCMeta, abstractmethod

import numpy
from numpy import array, ndarray

from vk2gpz.geom import util


class IProjection:
    pass


class Projection(IProjection, metaclass=ABCMeta):
    @abstractmethod
    def _latlong_to_2d(self, latlon: array) -> array:
        raise NotImplemented()

    def xyz_to_2d(self, coord: array) -> array:
        return self._latlong_to_2d(util.xyz_to_latlong(coord))

    def build(self, dome) -> ndarray:
        triangles = dome.get_faces()
        ver_per_face = dome.get_number_of_vertices_per_face()

        for v in dome.get_all_vertices():
            v.projected_coord = self.xyz_to_2d(v.coord)

        tmp2Dtri = []
        for i in range(len(triangles)):
            if i % ver_per_face == 0:
                # check the triangle is facing you or not.
                p1 = triangles[i].projected_coord
                p2 = triangles[i + 1].projected_coord
                p3 = triangles[i + 2].projected_coord
                if util.facing(p1, p2, p3):
                    tmp2Dtri.append([triangles[i + 0].id, triangles[i + 1].id, triangles[i + 2].id])

        return numpy.array(tmp2Dtri)

