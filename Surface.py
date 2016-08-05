import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib
import Quadtree as QT

class Surface:
    """Class for defining a custom surface"""
    def __init__(self, xy, uv, R):
        self.calculateQuadtree(xy);
        self.uv     = uv;
        self.dudv   = np.array([0.0,0.0]);
        self.R      = R;
        self.mass   = np.pi*np.power(R,2);
        self.xycent = np.average(xy,0);
        
    def calculateQuadtree(self, xy):
        mins      = np.array([np.min(xy[:,0]),np.min(xy[:,1])]);
        maxs      = np.array([np.max(xy[:,0]),np.max(xy[:,1])]);
        samp      = np.shape(xy)[0];
        tx        = np.diff(xy[:,0]); tx = np.append(tx,tx[samp-2]);
        ty        = np.diff(xy[:,1]); ty = np.append(ty,ty[samp-2]);
        NORM      = np.sqrt(np.power(tx,2) + np.power(ty,2));
        tx        = tx/NORM;
        ty        = ty/NORM;
        nx        = -ty;
        ny        = tx;
        self.xy   = xy;
        self.tang = np.transpose(np.vstack([tx,ty]));
        self.norm = np.transpose(np.vstack([nx,ny]));
        values    = np.hstack([self.tang,self.norm]);
        self.qt   = QT.QuadTree(xy,mins,maxs,values,depth=10,bucket=2);
        
    def translateQuadtree(self,dx,dy):
        self.qt.translateQuadtree(dx,dy);

    def set_xy(self,x,y):
        self.xy = np.transpose(np.array([x,y]));

    def set_uv(self,u,v):
        self.uv = np.array([u,v]);

    def increment_dudv(self,du,dv):
        self.dudv += np.array([du,dv]);

    def clear_dudv(self):
        self.dudv = np.array([0.0,0.0]);

    def calculateXYcent(self):
        self.xycent = np.average(self.xy,0);
