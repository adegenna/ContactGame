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
# Generate balls
levels  = 3;
samp    = 50;
R       = 0.05;
x       = R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy1     = np.transpose(np.array([x,y]));
uv1     = np.array([0.0,0.0]);
sepR    = 1.0;
bodies  = [];
for i in range(0,levels):
    for j in range(-20,20):
        X0   = (2*j+1)*R*1.1;
        Y0   = (2*i+1)*R*1.1 + 1.95;
        xc   = X0 + x;
        yc   = Y0 + y;
        xy   = np.transpose(np.array([xc,yc]));
        uv   = np.array([0.0,0.0]);
        surf = Surf.Surface(xy,uv,R,dt);
        bodies.append(surf);
print(np.shape(bodies))

# Generate balls above floor
# for i in range(0,1):
#     for j in range(0,1):
#         R       = 5*R;
#         X0      = -0.4; #-0.25 + (2*j+1)*R*1.02;
#         Y0      = 0.7; #9*R  + (2*i+1)*R*1.02;
#         x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
#         y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
#         xy      = np.transpose(np.array([x,y]));
#         uv      = np.array([-2.0,-2.0]);
#         surf    = Surf.Surface(xy,uv,R,dt);
#         bodies.append(surf);

# Generate wall
wall   = [];
samp   = 500;
x1     = 0.5*np.linspace(-1,1,samp);
y1     = 0.5*np.linspace(0,2,samp);
xw     = np.hstack([x1,0.5*np.ones(samp),-x1,-0.5*np.ones(samp)])
yw     = np.hstack([np.zeros(samp),y1,1*np.ones(samp),y1[::-1]]);
xyw    = np.transpose(5*np.array([xw,yw]));
uvw    = np.array([0.0,0.0]);
wall1   = Surf.Surface(xyw,uvw,0,0);
wall.append(wall1);
np.savetxt(OUTDIR + '/WALL.csv',xyw,delimiter=',');
x1     = 0.0 + 0.4*np.cos(np.linspace(-np.pi,np.pi,samp));
y1     = 1.5 + 0.4*np.sin(np.linspace(-np.pi,np.pi,samp));
xyw    = np.transpose(np.array([x1,y1]));
uvw    = np.array([0.0,0.0]);
wall2   = Surf.Surface(xyw,uvw,0,0);
wall.append(wall2);
np.savetxt(OUTDIR + '/WALL2.csv',xyw,delimiter=',');
x2     = -1.5 + 0.4*np.cos(np.linspace(-np.pi,np.pi,samp));
y2     = 1.5 + 0.4*np.sin(np.linspace(-np.pi,np.pi,samp));
xyw    = np.transpose(np.array([x2,y2]));
uvw    = np.array([0.0,0.0]);
wall3   = Surf.Surface(xyw,uvw,0,0);
wall.append(wall3);
np.savetxt(OUTDIR + '/WALL3.csv',xyw,delimiter=',');
x3     = 1.5 + 0.4*np.cos(np.linspace(-np.pi,np.pi,samp));
y3     = 1.5 + 0.4*np.sin(np.linspace(-np.pi,np.pi,samp));
xyw    = np.transpose(np.array([x3,y3]));
uvw    = np.array([0.0,0.0]);
wall4   = Surf.Surface(xyw,uvw,0,0);
wall.append(wall4);
np.savetxt(OUTDIR + '/WALL4.csv',xyw,delimiter=',');
print(np.size(wall));


# ****************************************************
# MECHANICS
# ****************************************************

# Physics timestepping
steps  = 50000;
for i in range(0,steps):
    print(str(i*dt) + '/' + str(steps*dt));
    STATE = np.zeros([np.size(bodies),2]);
    physics.potentialMethod(bodies,wall,dt);
    for j in range(0,np.size(bodies)):
        STATE[j,:] = bodies[j].xycent;
    np.savetxt(OUTDIR + '/T' + str((i+1)*dt) + '.csv',STATE,delimiter=',');
