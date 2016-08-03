import numpy as np
import Quadtree as QT
import Surface as Surf

def contact(body1, body2):
    XY = np.array([]);
    for i in range(0,body1.qt.num):
        xy1 = body2.qt.search_QT(body1.qt.data[i,0],body1.qt.data[i,1]);
        XY = np.append(XY,xy1);

    XY = np.reshape(XY,(np.size(XY)/2,2));
    return XY;
