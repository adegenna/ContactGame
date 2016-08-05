import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib

class QuadTree:
    """Simple QuadTree class"""

    # class initialization function
    def __init__(self, data, mins, maxs, values, depth, bucket):
        self.data    = np.asarray(data)
        self.num     = np.shape(data)[0];
        self.bucket  = bucket;
        self.depth   = depth;
        self.values  = values;
        self.contact = np.zeros(self.num);
        
        # data should be two-dimensional
        assert self.data.shape[1] == 2

        self.mins = np.asarray(mins)
        self.maxs = np.asarray(maxs)
        self.sizes = self.maxs - self.mins

        self.children = []

        mids = 0.5*(self.mins + self.maxs)
        xmin, ymin = self.mins
        xmax, ymax = self.maxs
        xmid, ymid = mids
        self.mids = mids;

        if depth > 0:
            # split the data into four quadrants
            ind1 = (data[:, 0] < mids[0])  & (data[:, 1] < mids[1]);
            ind2 = (data[:, 0] < mids[0])  & (data[:, 1] >= mids[1]);
            ind3 = (data[:, 0] >= mids[0]) & (data[:, 1] < mids[1]);
            ind4 = (data[:, 0] >= mids[0]) & (data[:, 1] >= mids[1]);
            data_q1 = data[ind1]; values1 = values[ind1,:];
            data_q2 = data[ind2]; values2 = values[ind2,:];
            data_q3 = data[ind3]; values3 = values[ind3,:];
            data_q4 = data[ind4]; values4 = values[ind4,:];
            num1 = np.size(ind1);
            num2 = np.size(ind2);
            num3 = np.size(ind3);
            num4 = np.size(ind4);

            # recursively build a quad tree on each quadrant which has data
            if ((data_q1.shape[0] > 0) & (num1 > self.bucket)):
                self.children.append(QuadTree(data_q1,
                                              [xmin, ymin], [xmid, ymid], values1,
                                              depth - 1,self.bucket))
            if ((data_q2.shape[0] > 0) & (num2 > self.bucket)):
                self.children.append(QuadTree(data_q2,
                                              [xmin, ymid], [xmid, ymax], values2, 
                                              depth - 1,self.bucket))
            if ((data_q3.shape[0] > 0) & (num3 > self.bucket)):
                self.children.append(QuadTree(data_q3, 
                                              [xmid, ymin], [xmax, ymid], values3,
                                              depth - 1,self.bucket))
            if ((data_q4.shape[0] > 0) & (num4 > self.bucket)):
                self.children.append(QuadTree(data_q4,
                                              [xmid, ymid], [xmax, ymax], values4,
                                              depth - 1,self.bucket))
                    
    def draw_rectangle(self, ax):
        """Recursively plot a visualization of the quad tree region"""
        if self.num <= self.bucket:
            rect = plt.Rectangle(self.mins, *self.sizes, zorder=2,
                                 ec='k', fc='red', alpha = 0.2);
            ax.add_patch(rect)
        for child in self.children:
            num = child.num;
            child.draw_rectangle(ax)

    def search_QT(self, x, y):
        """Search quadtree for a query point"""
        qt = self;
        count = 0;
        # First check that (x,y) is even in bounding QT box
        flag = ( (x >= qt.mins[0]) & (x < qt.maxs[0]) & (y >= qt.mins[1]) & (y < qt.maxs[1]) );
        while ((qt.children != []) & (qt.num != 0) & (flag == True)):
            for child in qt.children:
                flag = ( (x >= child.mins[0]) & (x < child.maxs[0]) & (y >= child.mins[1]) & (y < child.maxs[1]) );
                if (flag == True):
                    qt = child;
                    break;
        if (flag == True):
            mids = qt.mids;
            if (mids.ndim == 2):
                xmid = np.average(mids[:,0]);
                ymid = np.average(mids[:,1]);
                vals = np.average(qt.values,0);
            else:
                xmid = mids[0];
                ymid = mids[1];
                vals = np.average(qt.values,0);
            data = np.hstack([xmid,ymid,vals]);
        else:
            data = [];
        return data;

    
    def translateChildren(self,dx,dy):
        # Function to translate quadtree child nodes via recursion
        for child in self.children:
            child.data += np.array([dx,dy]);
            child.mins += np.array([dx,dy]);
            child.maxs += np.array([dx,dy]);
            child.mids += np.array([dx,dy]);
            child.translateChildren(dx,dy);
            
    def translateQuadtree(self,dx,dy):
        # Function to translate quadtree

        # Initialization: translate mother node
        qt       = self;
        qt.data += np.array([dx,dy]);
        qt.mins += np.array([dx,dy]);
        qt.maxs += np.array([dx,dy]);
        qt.mids += np.array([dx,dy]);
        # Recursive calls to translateChildren
        qt.translateChildren(dx,dy);


def draw_grid(ax, xlim, ylim, Nx, Ny, **kwargs):
    """ draw a background grid for the quad tree"""
    for x in np.linspace(xlim[0], xlim[1], Nx):
        ax.plot([x, x], ylim, **kwargs)
    for y in np.linspace(ylim[0], ylim[1], Ny):
        ax.plot(xlim, [y, y], **kwargs)

#------------------------------------------------------------
