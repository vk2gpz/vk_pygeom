import numpy as np
import numpy.linalg as la
from numpy import array

_2pi = np.pi * 2
x_axis = np.array([1, 0, 0])
y_axis = np.array([0, 1, 0])
z_axis = np.array([0, 0, 1])


def angle(a: array, b: array) -> (float, float):
    """
    Computes the angle between two vectors.

    :param a: 1st vector
    :param b: 2nd vector
    :return: angle in radian and degree
    """
    inner = np.inner(a, b)
    norms = la.norm(a) * la.norm(b)

    cos = inner / norms
    rad = np.arccos((np.clip(cos, -1.0, 1.0)))
    deg = np.rad2deg(rad)

    return rad, deg


def facing(p1: array, p2: array, p3: array) -> bool: # List[GeodesicVertex], index: int) -> bool:
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = np.cross(v1, v2)
    return v3 > 0


def xyz_to_latlong(coord: array) -> array:
    lat = np.arcsin(coord[2])
    lon = np.arctan2(coord[1], coord[0])
    return np.array([lat, lon])


def spherical_to_xyz(latitude, longitude) -> array:
    y = np.cos(latitude)
    x = np.sin(latitude) * np.sin(longitude)
    z = np.sin(latitude) * np.cos(longitude)
    return np.array([x, y, z])


def xyz_to_spherical(coord: array) -> array:
    latlong: array = np.array([0, 0])
    if coord[1] > 0.9999:
        latlong = np.array([0, 0])
    elif coord[1] < -0.9999:
        latlong = np.array([np.pi, 0])
    else:
        latlong[0] = angle(coord, y_axis)[0]
        temp = np.array([coord[0], 0, coord[2]])
        latlong[1] = angle(temp, z_axis)[0]
        if coord[0] < 0:
            latlong[1] = _2pi - latlong[1]
    return latlong


def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    mat = np.eye(3, 3)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.)
    b, c, d = -axis * np.sin(theta / 2.)

    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c), 0],
                     [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b), 0],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c, 0],
                     [0, 0, 0, 1]])


def rotation_matrix4(axis: np.ndarray, theta: float) -> np.ndarray:
    matrix3 = rotation_matrix(axis, theta)
    return np.array([
        [[matrix3[0][0], matrix3[0][1], matrix3[0][2], 0],
         [matrix3[1][0], matrix3[1][1], matrix3[1][2], 0],
         [matrix3[2][0], matrix3[2][1], matrix3[2][2], 0],
         [0, 0, 0, 1.]]
    ])
