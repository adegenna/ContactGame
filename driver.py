import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact
import physics

ax   = plt.gca();

# Initialize object
samp = 1000;
x    = 0.0 + 0.1*np.cos(np.linspace(-np.pi,np.pi,samp));
y    = 0.5 + 0.1*np.sin(np.linspace(-np.pi,np.pi,samp));
xy   = np.array([x,y]);
xy   = np.transpose(xy);
uv   = np.array([1.0,1.0]);
surf = Surf.Surface(xy,uv);

# Define wall boundaries
#xw    = np.hstack(np.array([-1*np.ones(samp),np.linspace(-1,1,samp),1*np.ones(samp)]));
#yw    = np.hstack(np.array([np.linspace(1,0,samp),np.zeros(samp),np.linspace(0,1,samp)]));
xw     = np.linspace(-1,1,2*samp);
yw     = np.abs(xw);
xyw   = np.array([xw,yw]);
xyw   = np.transpose(xyw);
uvw   = np.array([0.0,0.0]);
surfw = Surf.Surface(xyw,uvw);
np.savetxt('OUT/WALL.csv',xyw,delimiter=',');

# Physics timestepping
plt.plot(xw,yw,'k');
plt.plot(x,y,'k');
steps = 200;
dt    = 0.01;
for i in range(0,steps):
    print((i+1)*dt)
    physics.timestep(surf,surfw,dt);
    plt.scatter(np.average(surf.xy[:,0]),np.average(surf.xy[:,1]),s=50,c='r');
    np.savetxt('OUT/T' + str((i+1)*dt) + '.csv',surf.xy,delimiter=',');

# Calculate contact
# xy1 = contact(surf,surfw);

# if (xy1 != []):
#     print np.shape(xy1)
#     for i in range(0,np.shape(xy1)[0]):
#         plt.scatter(xy1[i,0],xy1[i,1],s=100,c='b');
plt.xlim([-1.2,1.2]);
plt.ylim([-0.2,1.2]);
ax.set_aspect('equal');

plt.show()
