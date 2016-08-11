import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact
import physics

# ****************************************************
# GEOMETRY
# ****************************************************

dt    = 0.001;
# Generate cone of close-packed, equal sized balls
# Generate first ball at bottom of cone
levels  = 2;
samp    = 1000;
R       = 0.05;
X0      = 0;
Y0      = 2*R*1.1;
x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy1     = np.transpose(np.array([x,y]));
uv1     = np.array([0.0,0.0]);
surf1   = Surf.Surface(xy1,uv1,R,dt);
bodies  = [surf1];

for i in range(0,levels-1):
    # Generate balls left and right of leftmost ball
    ind     = (i)*(i+1)/2;
    X0      = bodies[ind].xycent[0] - R*1.1;
    Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*1.1;
    x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
    y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
    xy      = np.transpose(np.array([x,y]));
    uv      = np.array([0.0,0.0]);
    surf    = Surf.Surface(xy,uv,R,dt);
    bodies.append(surf);
    X0      = bodies[ind].xycent[0] + R*1.1;
    Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*1.1;
    x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
    y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
    xy      = np.transpose(np.array([x,y]));
    uv      = np.array([0.0,0.0]);
    surf    = Surf.Surface(xy,uv,R,dt);
    bodies.append(surf);
    # Generate remainder of triangle
    for j in range(0,i):
        ind += 1;
        X0      = bodies[ind].xycent[0] + R*1.1;
        Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*1.1;
        x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
        y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
        xy      = np.transpose(np.array([x,y]));
        uv      = np.array([0.0,0.0]);
        surf    = Surf.Surface(xy,uv,R);
        bodies.append(surf);

# Generate a single ball above cone
X0      = 0.1;
Y0      = bodies[-1].xycent[1] + 4*R*1.1;
x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy      = np.transpose(np.array([x,y]));
uv      = np.array([0.0,-4.0]);
surf    = Surf.Surface(xy,uv,R,dt);
bodies.append(surf);

# Generate wall
# xw     = np.linspace(-1.2*R*(levels+1),1.2*R*(levels+1),2*samp);
# yw     = np.abs(xw)*np.sqrt(3.0);
# xyw    = np.transpose(np.array([xw,yw]));
# uvw    = np.array([0.0,0.0]);
# wall   = Surf.Surface(xyw,uvw,0);
# np.savetxt('OUT/WALL.csv',xyw,delimiter=',');

# ****************************************************
# MECHANICS
# ****************************************************

# Physics timestepping
steps = 1000;
for i in range(0,steps):
    print(str(i*dt) + '/' + str(steps*dt));
    physics.forwardEuler_V2(bodies,dt);
    for j in range(0,np.size(bodies)):
        np.savetxt('OUT/T' + str((i+1)*dt) + '_' + str(j) + '.csv',bodies[j].xy,delimiter=',');


