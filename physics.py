import numpy as np
import Quadtree as QT
import Surface as Surf

diff = 0.7;   # Global dissipation factor
GRAV = -9.81; # Gravity

def contact(body1, body2):
    # Function to detect whether there is contact between two bodies
    XY   = np.array([]);
    cols = np.shape(body1.qt.values)[1] + np.shape(body1.qt.data)[1];
    for i in range(0,body1.qt.num):
        xy1 = body2.qt.search_QT(body1.qt.data[i,0],body1.qt.data[i,1]);
        XY = np.append(XY,xy1);
    XY = np.reshape(XY,(np.size(XY)/cols,cols));
    return XY;

def wallCollision(bodies, wall, dt):
    # Function to calculate collision between bodies and the wall
    
    num = np.size(bodies);
    for i in range(0,num):    # Index through all bodies
        body = bodies[i];
        xy   = contact(body,wall);
        if (xy.any()):
            # Calculate average collision point properties
            xmean       = np.average(xy[:,0]);
            ymean       = np.average(xy[:,1]);
            tx          = np.average(xy[:,2]);
            ty          = np.average(xy[:,3]);
            nx          = np.average(xy[:,4]);
            ny          = np.average(xy[:,5]);
            # Elastic wall collision + wall does not move
            u           = body.uv[0];
            v           = body.uv[1];
            tvel        = u*tx + v*ty;
            nvel        = u*nx + v*ny;
            tvelBounce  = tvel;
            if (nvel < 0):
                nvelBounce  = -diff*nvel;
            else:
                nvelBounce  = diff*nvel;
            uBounce     = tx*tvelBounce + nx*nvelBounce;
            vBounce     = ty*tvelBounce + ny*nvelBounce;
            du          = uBounce - u;
            dv          = vBounce - v;
        else:
            # No collision
            du          = 0;
            dv          = 0;
        # Update body position/velocity
        bodies[i].increment_dudv(du,dv);

def bodyCollision(bodies, dt):
    # Function to calculate collision between bodies
    
    num = np.size(bodies);
    for i in range(0,num-1):    # Index through all bodies
        for j in range(i+1,num):
            body1 = bodies[i];
            body2 = bodies[j];
            xy    = contact(body1,body2);
            if (xy.any()):
                # Calculate average collision point properties
                xmean       = np.average(xy[:,0]);
                ymean       = np.average(xy[:,1]);
                tx          = np.average(xy[:,2]);
                ty          = np.average(xy[:,3]);
                nx          = np.average(xy[:,4]);
                ny          = np.average(xy[:,5]);
                # Elastic collision
                mass1       = body1.mass;
                mass2       = body2.mass;
                v21         = body2.uv - body1.uv;
                x21         = body2.xycent - body1.xycent;
                VEL1F       = body1.uv - 2*mass2/(mass1+mass2)*np.inner(v21,x21)/np.power(np.linalg.norm(x21),2.0)*(-x21);
                VEL2F       = body2.uv - 2*mass1/(mass1+mass2)*np.inner(v21,x21)/np.power(np.linalg.norm(x21),2.0)*(x21);
                du1         = VEL1F[0] - body1.uv[0];
                dv1         = VEL1F[1] - body1.uv[1];
                du2         = VEL2F[0] - body2.uv[0];
                dv2         = VEL2F[1] - body2.uv[1];
            else:
                # No collision
                du1         = 0;
                dv1         = 0;
                du2         = 0;
                dv2         = 0;
            # Update body position/velocity
            bodies[i].increment_dudv(du1,dv1);
            bodies[j].increment_dudv(du2,dv2);
            
def forwardEuler(bodies, wall, dt):
    # Function to update body's position based on contact
    
    # Calculate body-wall collisions
    wallCollision(bodies, wall, dt);
    # Calculate body-body collisions
    bodyCollision(bodies, dt);
    # Update body positions
    num = np.size(bodies);
    for i in range(0,num):
        u  = bodies[i].dudv[0] + bodies[i].uv[0];
        v  = bodies[i].dudv[1] + bodies[i].uv[1];
        dx = u*dt;
        dy = v*dt;
        bodies[i].set_xy(bodies[i].xy[:,0] + dx , bodies[i].xy[:,1] + dy);
        bodies[i].set_uv(u                 , v + GRAV*dt);
        bodies[i].calculateXYcent();
        bodies[i].clear_dudv();
        bodies[i].calculateQuadtree(bodies[i].xy);
