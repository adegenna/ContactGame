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
        if ( (xy.any()) ):
            mass        = body.mass;
            # Calculate average collision point properties
            xmean       = np.average(xy[:,0]);
            ymean       = np.average(xy[:,1]);
            tx          = np.average(xy[:,2]);
            ty          = np.average(xy[:,3]);
            nx          = np.average(xy[:,4]);
            ny          = np.average(xy[:,5]);
            # Wall condition: reflection of normal velocity
            u           = body.uv[0];
            v           = body.uv[1];
            tvel        = u*tx + v*ty;
            nvel        = u*nx + v*ny;
            tvelBounce  = tvel;
            if (nvel <= 0):
                nvelBounce  = -diff*nvel;
            else:
                nvelBounce  = diff*nvel;
            uBounce     = tx*tvelBounce + nx*nvelBounce;
            vBounce     = ty*tvelBounce + ny*nvelBounce;
            du          = uBounce - u;
            dv          = vBounce - v;
            # Calculate impulse from wall condition and update J_i
            J_contact   = body.J;
            J_grav      = GRAV*np.array([0,1])*dt;
            Jt_contact  = J_contact[0]*tx + J_contact[1]*ty;
            Jn_contact  = J_contact[0]*nx + J_contact[1]*ny;
            Jt_grav     = J_grav[0]*tx + J_grav[1]*ty;
            Jn_grav     = J_grav[0]*nx + J_grav[1]*ny;
            dut_reflect = 0.0;
            dun_reflect = 2.0*nvelBounce;
            Jt_wall     = dut_reflect;
            Jn_wall     = dun_reflect - (Jn_contact + Jn_grav);
            Jt          = Jt_contact + Jt_grav + Jt_wall;
            Jn          = Jn_contact + Jn_grav + Jn_wall;
            Jx          = tx*Jt + nx*Jn;
            Jy          = ty*Jt + ny*Jn;
            J           = np.array([Jx,Jy]);
            bodies[i].set_J(J);
            # Apply impulse from wall condition to bodies colliding simultaneously with J_i
            numContact = np.size(bodies[i].bodyContacts);
            for k in range(0,numContact):
                print("COND")
                ind        = int(bodies[i].bodyContacts[k]);
                J0         = bodies[ind].J;
                xij        = (bodies[i].xycent - bodies[ind].xycent);
                xij        = xij/np.linalg.norm(xij);
                Jx_wall    = tx*Jt_wall + nx*Jn_wall;
                Jy_wall    = ty*Jt_wall + ny*Jn_wall;
                J_wallProj = np.inner(np.array([Jx_wall,Jy_wall]),xij);
                J_wall     = (bodies[i].mass/bodies[ind].mass)*J_wallProj*xij;
                Jnew       = J0 + J_wall;
                bodies[ind].set_J(Jnew);
            du          = 0;
            dv          = 0;
        else: # No collision
            J_grav      = GRAV*np.array([0,1])*dt;
            bodies[i].increment_J(J_grav);
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
            # Calculate average collision point properties
            xmean       = np.average(xy[:,0]);
            ymean       = np.average(xy[:,1]);
            tx          = np.average(xy[:,2]);
            ty          = np.average(xy[:,3]);
            nx          = np.average(xy[:,4]);
            ny          = np.average(xy[:,5]);
            uv2         = body2.uv;
            tvel        = uv2[0]*tx + uv2[1]*ty;
            nvel        = uv2[0]*nx + uv2[1]*ny;
            if ( (xy.any()) & (nvel < 0.0) ) :
                # Elastic collision
                mass1       = body1.mass;
                mass2       = body2.mass;
                v21         = body2.uv - body1.uv;
                x21         = body2.xycent - body1.xycent;
                x21         = (x21/np.linalg.norm(x21))*(body1.R + body2.R); # Rescale x21 to fix glitching behavior
                J1          = diff*(-2*mass2/(mass1+mass2)*np.inner(v21,x21)/np.power(np.linalg.norm(x21),2.0)*(-x21));
                J2          = diff*(-2*mass1/(mass1+mass2)*np.inner(v21,x21)/np.power(np.linalg.norm(x21),2.0)*(x21));
                VEL1F       = body1.uv + J1;
                VEL2F       = body2.uv + J2;
                du1         = diff*(VEL1F[0] - body1.uv[0]);
                dv1         = diff*(VEL1F[1] - body1.uv[1]);
                du2         = diff*(VEL2F[0] - body2.uv[0]);
                dv2         = diff*(VEL2F[1] - body2.uv[1]);
                bodies[i].append_contacts(j);
                bodies[j].append_contacts(i);
            else:
                # No collision
                du1         = 0;
                dv1         = 0;
                du2         = 0;
                dv2         = 0;
                J1          = np.array([0,0]);
                J2          = np.array([0,0]);
            # Update body position/velocity
            bodies[i].increment_dudv(du1,dv1);
            bodies[j].increment_dudv(du2,dv2);
            bodies[i].increment_J(J1);
            bodies[j].increment_J(J2);
            
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
        bodies[i].set_uv(u                      , v + GRAV*dt);
        bodies[i].calculateXYcent();
        bodies[i].clear_dudv();
        bodies[i].translateQuadtree(dx,dy);
        #bodies[i].calculateQuadtree(bodies[i].xy);

def verletIntegration(bodies, wall, dt):
    # Verlet integration scheme

    # Calculate body-body collisions
    bodyCollision(bodies, dt);
    # Calculate body-wall collisions
    wallCollision(bodies, wall, dt);
    # Update body positions
    num    = np.size(bodies);
    sizeXY = np.shape(bodies[0].xy)[0];
    for i in range(0,num):
        J              = bodies[i].J;
        a_tTimesDT     = J;
        xy_tPlusDT     = 2*bodies[i].xy - bodies[i].xyPrev + np.tile(a_tTimesDT*dt,[sizeXY,1]);
        bodies[i].set_xyPrev(bodies[i].xy[:,0],bodies[i].xy[:,1]);
        bodies[i].set_xy(xy_tPlusDT[:,0],xy_tPlusDT[:,1]);
        xy0            = bodies[i].xycent;
        bodies[i].calculateXYcent();
        dx             = bodies[i].xycent[0] - xy0[0];
        dy             = bodies[i].xycent[1] - xy0[1];
        u              = dx/dt;
        v              = dy/dt;
        bodies[i].set_uv(u,v);
        bodies[i].clear_dudv();
        bodies[i].clear_J();
        bodies[i].translateQuadtree(dx,dy);
        bodies[i].clear_contacts();
        #bodies[i].calculateQuadtree(bodies[i].xy);
