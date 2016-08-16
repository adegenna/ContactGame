import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact
import physics

OUTDIR = "/home/adegenna/ContactGame/OUT";

# ****************************************************
# GEOMETRY
# ****************************************************

dt    = 1.0e-4;
# Generate cone of close-packed, equal sized balls
# Generate first ball at bottom of cone
levels  = 2;
samp    = 50;
R       = 0.05;
x       = R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy1     = np.transpose(np.array([x,y]));
uv1     = np.array([0.0,0.0]);
sepR    = 1.0;
bodies  = [];
for i in range(0,levels):
    for j in range(-4,4):
        X0   = (2*j+1)*R;
        Y0   = (2*i+1)*R*1.01 + 0.02;
        xc   = X0 + x;
        yc   = Y0 + y;
        xy   = np.transpose(np.array([xc,yc]));
        uv   = np.array([0.0,0.0]);
        surf = Surf.Surface(xy,uv,R,dt);
        bodies.append(surf);
print(np.shape(bodies))

# Generate balls above floor
X0      = 2*R;
Y0      = (2*(levels+1)+1)*R*1.01;
x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy      = np.transpose(np.array([x,y]));
uv      = np.array([-1.0,-4.0]);
surf    = Surf.Surface(xy,uv,R,dt);
bodies.append(surf);

# Generate wall
samp   = 500;
x1     = 10*R*np.linspace(-1,1,samp);
y1     = 10*R*np.linspace(-0,2,samp);
xw     = np.hstack([x1,10*R*np.ones(samp),-x1,-10*R*np.ones(samp)])
yw     = np.hstack([np.zeros(samp),y1,20*R*np.ones(samp),y1[::-1]]) + 0.02;
xyw    = np.transpose(np.array([xw,yw]));
print(np.shape(xyw))
uvw    = np.array([0.0,0.0]);
wall   = Surf.Surface(xyw,uvw,0,0);
np.savetxt(OUTDIR + '/WALL.csv',xyw,delimiter=',');

# ****************************************************
# MECHANICS
# ****************************************************

# Physics timestepping
steps = 50000;
for i in range(0,steps):
    print(str(i*dt) + '/' + str(steps*dt));
    physics.potentialMethod(bodies,wall,dt);
    for j in range(0,np.size(bodies)):
        np.savetxt(OUTDIR + '/T' + str((i+1)*dt) + '_' + str(j) + '.csv',bodies[j].xycent,delimiter=',');
