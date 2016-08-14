import numpy as np
import Quadtree as QT
import Surface as Surf
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=132)

diff = 0.7;   # Global dissipation factor
GRAV = -9.81; # Gravity

def contact(body1, body2, CONT, I, J):
    # Function to detect whether there is contact between two bodies
    XY   = np.array([]);
    cols = np.shape(body1.qt.values)[1] + np.shape(body1.qt.data)[1];
    for i in range(0,body1.qt.num):
        xy1 = body2.qt.search_QT(body1.qt.data[i,0],body1.qt.data[i,1]);
        XY  = np.append(XY,xy1);
    # Rule out collisions with positive normal velocity
    if XY.any():
        XY = np.reshape(XY,(np.size(XY)/cols,cols));
        tx          = np.average(XY[:,2]);
        ty          = np.average(XY[:,3]);
        nx          = np.average(XY[:,4]);
        ny          = np.average(XY[:,5]);
        uv2         = body2.uv;
        tvel        = uv2[0]*tx + uv2[1]*ty;
        nvel        = uv2[0]*nx + uv2[1]*ny;
        # Remember contacting bodies
        if (nvel <= 0):
            CONT[I,J] = 1;
            CONT[J,I] = 1;
    return XY,CONT;

def contact_V2(body1,body2,CONT,I,J):
    # Function to detect whether there is contact between two bodies
    XY   = np.array([]);
    dx   = body1.xycent[0] - body2.xycent[0];
    dy   = body1.xycent[1] - body2.xycent[1];
    dist = np.sqrt(np.power(dx,2) + np.power(dy,2));
    if ( (dist < (body1.R+body2.R)) ):
        CONT[I,J] = 1;
        CONT[J,I] = 1;
    return XY,CONT;

def contactWall(body, wall):
    # Function to detect whether there is contact between two bodies
    XY   = np.array([]);
    cols = np.shape(body.qt.values)[1] + np.shape(body.qt.data)[1];
    for i in range(0,body.qt.num):
        xy1 = wall.qt.search_QT(body.qt.data[i,0],body.qt.data[i,1]);
        XY  = np.append(XY,xy1);
    # Rule out collisions with positive normal velocity
    if XY.any():
        XY = np.reshape(XY,(np.size(XY)/cols,cols));
        tx          = np.average(XY[:,2]);
        ty          = np.average(XY[:,3]);
        nx          = np.average(XY[:,4]);
        ny          = np.average(XY[:,5]);
        uv2         = body.uv;
        tvel        = uv2[0]*tx + uv2[1]*ty;
        nvel        = uv2[0]*nx + uv2[1]*ny;
        # Remember contacting bodies
        if (nvel <= 0):
            XY_return = XY;
        else:
            XY_return = np.array([]);
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
    
    num  = np.shape(bodies)[0];
    # Calculate all contacts between all bodies
    CONT = np.zeros([num,num]);
    for i in range(0,num-1):
        for j in range(i+1,num):
            body1   = bodies[i];
            body2   = bodies[j];
            XY,CONT = contact_V2(body1,body2,CONT,i,j);
    return CONT;
            
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



# Jacobian calculation
def jacobian(x,E,bodies,num):
    # Calculate initial/final energies
    E0 = 0.0; Ef = 0.0;
    M = np.zeros(num);
    X = np.zeros([2,num]);
    V = np.zeros([2,num]);
    for i in range(0,num):
        M[i]   = bodies[i].mass;
        X[:,i] = bodies[i].xycent;
        V[:,i] = bodies[i].uv;
    for i in range(0,num):
        E0 += 0.5*M[i]*(np.power(V[0,i],2) + np.power(V[1,i],2));
        Ef += 0.5*M[i]*(np.power(x[2*i],2) + np.power(x[2*i+1],2));
    # Construct mass matrix
    MASS = np.zeros([2*num,2*num]);
    for i in range(0,num):
        MASS[2*i,2*i]     = M[i];
        MASS[2*i+1,2*i+1] = M[i];
    # Construct impulse line-of-action matrix
    rows = 2*num; cols = num*(num-1)/2;
    CONT = np.zeros([rows,cols]);
    count = 0;
    for i in range(0,num):
        for j in range(i+1,num):
            K = E[i,j];
            e = np.squeeze((X[:,i]-X[:,j])/np.linalg.norm(X[:,i]-X[:,j]));
            CONT[2*i:2*i+2,count] = K*e;
            CONT[2*j:2*j+2,count] = -K*e;
            count += 1;
    CONT = CONT[:,~(CONT==0).all(0)]; # Remove zero columns from line-of-action matrix (bodies not in contact)
    # Construct matrix derivative of energy equation
    ENERGY = np.zeros(2*num);
    for i in range(0,num):
        ENERGY[2*i]   = M[i]*x[2*i];
        ENERGY[2*i+1] = M[i]*x[2*i+1];
    DF1 = np.concatenate([MASS,CONT],axis=1);
    DF2 = np.concatenate([ENERGY,np.zeros(np.shape(CONT)[1])]);
    DF  = np.vstack([DF1,DF2]);
    # Calculate function
    f = np.zeros(2*num+1);
    for i in range(0,num):
        f[2*i]   = M[i]*x[2*i]   - M[i]*V[0,i] + np.dot(CONT[2*i,:]  , x[2*num:]);
        f[2*i+1] = M[i]*x[2*i+1] - M[i]*V[1,i] + np.dot(CONT[2*i+1,:], x[2*num:]);
    f[-1] = Ef - E0;
    # Calculate DJ
    DJ   = np.dot(np.transpose(DF),f);
    fval = np.dot(np.transpose(f),f);
    #print(DJ)
    return DJ,fval;

def potentialFunction(r,R):
    # Potential function for collision force calculation

    return np.power(10.0,6.0)*np.abs(R-r);
    #return np.power(10.0,11)*np.power(R-r,1.5);

def wallBoundaryCondition(bodies,wall,dt):
    # Wall reflection boundary condition
    num = np.size(bodies);
    for i in range(0,num):    # Index through all bodies
        body = bodies[i];
        xy   = contactWall(body,wall);
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
            bodies[i].set_uv(uBounce,vBounce);
    
def potentialMethod(bodies,wall,dt):
    num      = np.size(bodies);
    CONT     = bodyCollision(bodies,dt);
    contacts = int(np.sum(np.triu(CONT,1)));
    if (contacts > 0):
        I,J      = CONT.nonzero();
        bodInd   = np.unique(np.array([I,J]));
        for i in range(0,contacts):
            xij = bodies[J[i]].xycent - bodies[I[i]].xycent;
            rij = np.linalg.norm(xij);
            eij = xij/rij;
            Rij = 1.0*(bodies[I[i]].R + bodies[J[i]].R);
            Fij = potentialFunction(rij,Rij)*eij;
            dvj = Fij*dt/bodies[J[i]].mass;
            dvi = -Fij*dt/bodies[I[i]].mass;
            uvi = bodies[I[i]].uv + dvi;
            uvj = bodies[J[i]].uv + dvj;
            dxi = uvi[0]*dt;
            dyi = uvi[1]*dt;
            dxj = uvj[0]*dt;
            dyj = uvj[1]*dt;
            bodies[I[i]].set_uv(uvi[0],uvi[1]);
            bodies[J[i]].set_uv(uvj[0],uvj[1]);
    # Wall reflection boundary condition        
    wallBoundaryCondition(bodies,wall,dt);            
    for i in range(0,num):
        dx = bodies[i].uv[0]*dt;
        dy = bodies[i].uv[1]*dt;
        bodies[i].set_xy(bodies[i].xy[:,0] + dx , bodies[i].xy[:,1] + dy);
        bodies[i].calculateXYcent();
        bodies[i].clear_dudv();
        bodies[i].translateQuadtree(dx,dy);

def forwardEuler_V2(bodies,dt):
    # Function to update body's position based on contact

    num      = np.size(bodies);
    # Calculate body collision matrix
    CONT     = bodyCollision(bodies, dt);
    # Solve simultaneous collision physics
    steps    = 100000;
    eps      = 0.1;
    contacts = int(np.sum(np.triu(CONT,1)));
    if (contacts > 0):
        # Set up collision calculation
        x        = np.zeros(2*num + contacts);
        I,J      = CONT.nonzero();
        bodInd   = np.unique(np.array([I,J]));
        # Check velocities of colliding bodies
        flag = True;
        for i in range(0,contacts):
            flag = flag & np.all(bodies[bodInd[i]].uv==0);
        if (flag == False):
            for i in range(0,num):
                x[2*i]   = -bodies[i].uv[0];
                x[2*i+1] = -bodies[i].uv[1];
            for i in range(0,int(contacts)):
                x[2*I[i]]   = bodies[J[i]].uv[0];
                x[2*I[i]+1] = bodies[J[i]].uv[1];
                x[2*J[i]]   = bodies[I[i]].uv[0];
                x[2*J[i]+1] = bodies[I[i]].uv[1];
            for i in range(0,steps):
                DF,f = jacobian(x,CONT,bodies,num);
                x   += -eps*DF;
            print(x)
        else:
            for i in range(0,num):
                x[2*i]   = bodies[i].uv[0];
                x[2*i+1] = bodies[i].uv[1];
        # Update body velocities, positions
        for i in range(0,num):
            uv = np.array([x[2*i] , x[2*i+1]]);
            dx = uv[0]*dt;
            dy = uv[1]*dt;
            bodies[i].set_xy(bodies[i].xy[:,0] + dx , bodies[i].xy[:,1] + dy);
            bodies[i].set_uv(uv[0] , uv[1]);
            bodies[i].calculateXYcent();
            bodies[i].clear_dudv();
            bodies[i].translateQuadtree(dx,dy);
            #bodies[i].calculateQuadtree(bodies[i].xy);
    else:
        for i in range(0,num):
            dx = bodies[i].uv[0]*dt;
            dy = bodies[i].uv[1]*dt;
            bodies[i].set_xy(bodies[i].xy[:,0] + dx , bodies[i].xy[:,1] + dy);
            bodies[i].calculateXYcent();
            bodies[i].clear_dudv();
            bodies[i].translateQuadtree(dx,dy);
            #bodies[i].calculateQuadtree(bodies[i].xy);

    
            



        
