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
levels  = 4;
samp    = 500;
R       = 0.05;
X0      = 0;
Y0      = 2*R*1.1;
x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy1     = np.transpose(np.array([x,y]));
uv1     = np.array([0.0,-4.0]);
surf1   = Surf.Surface(xy1,uv1,R,dt);
bodies  = [surf1];

sepR    = 1.1;

for i in range(0,levels-1):
    # Generate balls left and right of leftmost ball
    ind     = (i)*(i+1)/2;
    X0      = bodies[ind].xycent[0] - R*sepR;
    Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*sepR;
    x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
    y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
    xy      = np.transpose(np.array([x,y]));
    uv      = np.array([0.0,-4.0]);
    surf    = Surf.Surface(xy,uv,R,dt);
    bodies.append(surf);
    X0      = bodies[ind].xycent[0] + R*sepR;
    Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*sepR;
    x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
    y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
    xy      = np.transpose(np.array([x,y]));
    uv      = np.array([0.0,-4.0]);
    surf    = Surf.Surface(xy,uv,R,dt);
    bodies.append(surf);
    # Generate remainder of triangle
    for j in range(0,i):
        ind += 1;
        X0      = bodies[ind].xycent[0] + R*sepR;
        Y0      = bodies[ind].xycent[1] + np.sqrt(3)*R*sepR;
        x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
        y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
        xy      = np.transpose(np.array([x,y]));
        uv      = np.array([0.0,-4.0]);
        surf    = Surf.Surface(xy,uv,R,dt);
        bodies.append(surf);

print(np.shape(bodies))
# Generate a single ball above cone
# X0      = 2*R;
# Y0      = bodies[-1].xycent[1] - 7*R*1.1;
# x       = X0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
# y       = Y0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
# xy      = np.transpose(np.array([x,y]));
# uv      = np.array([-1.0,-4.0]);
# surf    = Surf.Surface(xy,uv,R,dt);
# bodies.append(surf);

# Generate wall
xw     = np.hstack([np.linspace(-10*R,10*R,samp),10*R*np.ones(samp),np.linspace(10*R,-10*R,samp),-10*R*np.ones(samp)])
yw     = np.hstack([np.zeros(samp),np.linspace(0,20*R,samp),20*R*np.ones(samp),np.linspace(20*R,0,samp)]) + 0.02;
xyw    = np.transpose(np.array([xw,yw]));
print(np.shape(xyw))
uvw    = np.array([0.0,0.0]);
wall   = Surf.Surface(xyw,uvw,0,0);
np.savetxt(OUTDIR + '/WALL.csv',xyw,delimiter=',');

# ****************************************************
# MECHANICS
# ****************************************************

# Physics timestepping
steps = 10000;
for i in range(0,steps):
    print(str(i*dt) + '/' + str(steps*dt));
    physics.potentialMethod(bodies,wall,dt);
    for j in range(0,np.size(bodies)):
        np.savetxt(OUTDIR + '/T' + str((i+1)*dt) + '_' + str(j) + '.csv',bodies[j].xycent,delimiter=',');


