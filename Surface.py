import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib
import Quadtree as QT

class Surface:
    """Class for defining a custom surface"""
    def __init__(self, xy):
        mins    = np.array([np.min(xy[:,0]),np.min(xy[:,1])]);
        maxs    = np.array([np.max(xy[:,0]),np.max(xy[:,1])]);
        self.xy = xy;
        self.qt = QT.QuadTree(xy,mins,maxs,depth=10,bucket=5);
        
