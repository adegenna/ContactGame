import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact
import physics

ax   = plt.gca();

# Initialize objects
samp   = 1000;
R      = 0.02;
x      = 0.2 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y      = 1.0 + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy     = np.array([x,y]);
xy     = np.transpose(xy);
uv     = np.array([2.0,0.0]);
surf   = Surf.Surface(xy,uv,R);
R      = 0.05;
x      = 0.0 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y      = 1.0  + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy     = np.array([x,y]);
xy     = np.transpose(xy);
uv     = np.array([-1.0,0.0]);
surf2  = Surf.Surface(xy,uv,R);
bodies = [surf,surf2];
R      = 0.1;
x      = -0.2 + R*np.cos(np.linspace(-np.pi,np.pi,samp));
y      = 1.0  + R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy     = np.array([x,y]);
xy     = np.transpose(xy);
uv     = np.array([-2.0,0.0]);
surf3  = Surf.Surface(xy,uv,R);
bodies = [surf,surf2,surf3];

# Define wall boundaries
xw    = np.hstack(np.array([-1*np.ones(samp),np.linspace(-1,1,samp),1*np.ones(samp)]));
yw    = np.hstack(np.array([np.linspace(1,0,samp),np.zeros(samp),np.linspace(0,1,samp)]));
#xw    = np.linspace(-1,1,2*samp);
#yw    = np.power(xw,2);
xyw   = np.array([xw,yw]);
xyw   = np.transpose(xyw);
uvw   = np.array([0.0,0.0]);
wall  = Surf.Surface(xyw,uvw,0);
np.savetxt('OUT/WALL.csv',xyw,delimiter=',');

# Physics timestepping
plt.plot(xw,yw,'k');
plt.plot(x,y,'k');
steps = 500;
dt    = 0.01;
for i in range(0,steps):
    print((i+1)*dt)
    physics.forwardEuler(bodies,wall,dt);
    for j in range(0,np.size(bodies)):
        plt.scatter(bodies[j].xycent[0], bodies[j].xycent[1], s=50, c='b');
        np.savetxt('OUT/T' + str((i+1)*dt) + '_' + str(j) + '.csv',bodies[j].xy,delimiter=',');

# Calculate contact
# xy1 = contact(surf,wall);

# if (xy1 != []):
#     print np.shape(xy1)
#     for i in range(0,np.shape(xy1)[0]):
#         plt.scatter(xy1[i,0],xy1[i,1],s=100,c='b');
plt.xlim([-1.2,1.2]);
plt.ylim([-0.2,1.2]);
ax.set_aspect('equal');

plt.show()
