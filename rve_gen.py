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


WIDTH = 10
HEIGHT = 10
CIRCLE_RADIUS = 2
#TODO: convert all the lists of points to NdArray for algebraic operations
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

        
def _distance(point_1, point_2):
    
    dist = distance.euclidean(point_1.coordinates, point_2.coordinates)
    
    return dist

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
    

def calculate_all_distances(points):
    all_distances = np.array([])
    for point in points:
        #distance_(point, points)
        
        
    
if __name__ == "__main__":
    point_collection = []
    for i in range(100):
        new_point = Point()
        if not colliding(point_collection, new_point):
            point_collection.append(new_point)

#    for point in points:
#        point.move()
        
    show_points(point_collection)
    

    

