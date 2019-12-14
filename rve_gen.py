# -*- coding: utf-8 -*-
"""
Creating random micromechanical RVE of composites
author Hossein Ghayoor

"""

import random
import matplotlib.pyplot as plt
from scipy.spatial import distance
import scipy
import numpy as np
#import csv
#import copy


WIDTH = 100
HEIGHT = 100
CIRCLE_RADIUS = 2


class Point: # not sure if classes need docstring
    """
    Center points of non-colliding circles
    """
    def __init__(self):
        self.coordinates = np.random.rand(1, 2) * np.array([WIDTH, HEIGHT])
        
    def move(self):
        
        move_val = np.random.rand(1, 2)* np.array([3, 3]) # 3 is the move intensity
        self.coordinates += move_val
    
all_points = []

        
def colliding(points, new_point):
    radius = CIRCLE_RADIUS
    for point in points:
        if _distance(point, new_point) < 2*radius:
            return True

    return False

        
def _distance(point_1, point_2):
    dist = distance.euclidean(point_1.coordinates, point_2.coordinates)
    return dist


def show_points(points):
    all_xs, all_ys = [], []
    for point in points:
        all_xs.append(point.coordinates[0][0])
        all_ys.append(point.coordinates[0][1])


    plt.plot(all_xs, all_ys, "ro")
    plt.show()


def calculate_all_distances(points):
    # TODO: list comprehension to optimize the speed here
    all_distances = np.empty((0, len(points)), float)

    for point in points:
        one_point_to_rest = np.array([])
        for second_point in points:
            dist = _distance(point, second_point)
            one_point_to_rest = np.append(one_point_to_rest, dist)

        all_distances = np.vstack([all_distances, one_point_to_rest])
            
    return all_distances

def sum_n_closest_point_distance(n_points: int = 3, distance_arrays):
    """
    based on the result of this we find the secluded points to move them
    towards their neighbors
    """
    sum_n_close_points = []
    for dist_array in distance_arrays:
        close_point_ind = np.argpartition(distance_arrays, n_points)
        distance_sum = sum(distance_arrays[close_point_ind[:n_points]])
        sum_n_close_points.append(distance_sum)

    return sum_n_close_points
    
if __name__ == "__main__":
    point_collection = []
    for i in range(100):
        new_point = Point()
        if not colliding(point_collection, new_point):
            point_collection.append(new_point)

    all_dist = calculate_all_distances(point_collection)
    pass
#    for point in points:
#        point.move()

    show_points(point_collection)
    

    

