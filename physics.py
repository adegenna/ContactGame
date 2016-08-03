import numpy as np
import Quadtree as QT
import Surface as Surf

def contact(body1, body2):
    # Function to detect whether there is contact between two bodies
    XY   = np.array([]);
    cols = np.shape(body1.qt.values)[1] + np.shape(body1.qt.data)[1];
    for i in range(0,body1.qt.num):
        xy1 = body2.qt.search_QT(body1.qt.data[i,0],body1.qt.data[i,1]);
        XY = np.append(XY,xy1);

    XY = np.reshape(XY,(np.size(XY)/cols,cols));
    return XY;

def timestep(body, wall, dt):
    # Function to update body's position based on contact
    xy = contact(body,wall);
    diff = 1.0;
    if (xy.any()):
        xmean       = np.average(xy[:,0]);
        ymean       = np.average(xy[:,1]);
        tx          = np.average(xy[:,2]);
        ty          = np.average(xy[:,3]);
        nx          = np.average(xy[:,4]);
        ny          = np.average(xy[:,5]);
        u           = body.uv[0];
        v           = body.uv[1];
        tvel        = u*tx + v*ty;
        nvel        = u*nx + v*ny;
        tvelBounce  = diff*tvel;
        nvelBounce  = -diff*nvel;
        uBounce     = tx*tvelBounce + nx*nvelBounce;
        vBounce     = ty*tvelBounce + ny*nvelBounce;
        dx          = uBounce*dt;
        dy          = vBounce*dt;
        body.set_uv(uBounce,vBounce);
        print(tvel,nvel)
    else:
        dx          = body.uv[0]*dt;
        dy          = body.uv[1]*dt;
    body.set_xy(body.xy[:,0] + dx, body.xy[:,1] + dy);
    body.set_uv(body.uv[0],        body.uv[1] -9.81*dt);
    body.calculateQuadtree(body.xy,body.uv);
        

    return;
