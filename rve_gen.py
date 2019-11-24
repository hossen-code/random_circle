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
        #self.x = random.random()*10 # this should be the limits
        #self.y = random.random()*10
        
    def move(self):
        
        move_val = np.random.rand(1, 2)* np.array([3, 3]) # 3 is the move intensity
#        move_x = random.random()*2 + 1
#        move_y = random.random()*3 + 1
#        self.x += move_x
#        self.y += move_y
        self.coordinates += move_val
    

    
    
    
all_points = []

        
def colliding(points, new_point):
    radius = CIRCLE_RADIUS
    for point in points:
        if _distance(point, new_point) < 2*radius:
            return True

    return False

def show_points(points):
    all_xs, all_ys = [], []
    for point in points:
        all_xs.append(point.coordinates[0][0])
        all_ys.append(point.coordinates[0][1])


    plt.plot(all_xs, all_ys, "ro")
    plt.show()

        
def _distance(point_1, point_2):

    dist = distance.euclidean(point_1.coordinates, point_2.coordinates)

    return dist
    

def calculate_all_distances(points):
    # TODO: two numpy array have distance of all of them
    # TODO: in respect to each other
    all_distances = np.array([])

    for point in points:
        one_point_to_rest = np.array([])
        for second_point in points:
            dist = _distance(point, second_point)
            all_distances
            
    return all_distances
        
    
if __name__ == "__main__":
    point_collection = []
    for i in range(1000):
        new_point = Point()
        if not colliding(point_collection, new_point):
            point_collection.append(new_point)

    calculate_all_distances(point_collection)
#    for point in points:
#        point.move()

    show_points(point_collection)
    

    

