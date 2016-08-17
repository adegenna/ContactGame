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
levels  = 3;
samp    = 50;
R       = 0.01;
x       = R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy1     = np.transpose(np.array([x,y]));
uv1     = np.array([0.0,0.0]);
sepR    = 1.0;
bodies  = [];
for i in range(0,levels):
    for j in range(-24,24):
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
R       = 5*R;
X0      = 3*R;
Y0      = 3*R;
x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy      = np.transpose(np.array([x,y]));
uv      = np.array([-3.0,-3.0]);
surf    = Surf.Surface(xy,uv,R,dt);
bodies.append(surf);
# X0      = 18*R;
# Y0      = 10*R;
# x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
# y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
# xy      = np.transpose(np.array([x,y]));
# uv      = np.array([-3.0,-3.0]);
# surf    = Surf.Surface(xy,uv,R,dt);
# bodies.append(surf);
# X0      = 20*R;
# Y0      = 10*R;
# x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
# y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
# xy      = np.transpose(np.array([x,y]));
# uv      = np.array([-3.0,-3.0]);
# surf    = Surf.Surface(xy,uv,R,dt);
# bodies.append(surf);

# Generate wall
samp   = 500;
x1     = 0.5*np.linspace(-1,1,samp);
y1     = 0.5*np.linspace(0,2,samp);
xw     = np.hstack([x1,0.5*np.ones(samp),-x1,-0.5*np.ones(samp)])
yw     = np.hstack([np.zeros(samp),y1,1*np.ones(samp),y1[::-1]]) + 0.02;
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
STATE = np.zeros([np.size(bodies),2]);
for i in range(0,steps):
    print(str(i*dt) + '/' + str(steps*dt));
    physics.potentialMethod(bodies,wall,dt);
    for j in range(0,np.size(bodies)):
        STATE[j,:] = bodies[j].xycent;
    np.savetxt(OUTDIR + '/T' + str((i+1)*dt) + '.csv',STATE,delimiter=',');
