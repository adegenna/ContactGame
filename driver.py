import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact

ax   = plt.gca();

# Define object
samp = 1000;
x    = 0.25*np.cos(np.linspace(-np.pi,np.pi,samp));
y    = 0.5 + 0.25*np.sin(np.linspace(-np.pi,np.pi,samp));
xy   = np.array([x,y]);
xy   = np.transpose(xy);
surf = Surf.Surface(xy);

# Define wall boundaries
xw    = np.hstack(np.array([-1*np.ones(samp),np.linspace(-1,1,samp),1*np.ones(samp)]));
yw    = np.hstack(np.array([np.linspace(1,0,samp),np.zeros(samp),np.linspace(0,1,samp)]));
xyw   = np.array([xw,yw]);
xyw   = np.transpose(xyw);
surfw = Surf.Surface(xyw);
surf.qt.draw_rectangle(ax);
surfw.qt.draw_rectangle(ax);
plt.scatter(x,y,edgecolor='k',facecolor='none');
plt.scatter(xw,yw,edgecolor='k',facecolor='none');

# Calculate contact
xy1 = contact(surf,surfw);

if (xy1 != []):
    print np.shape(xy1)
    for i in range(0,np.shape(xy1)[0]):
        plt.scatter(xy1[i,0],xy1[i,1],s=100,c='b');
plt.xlim([-1.2,1.2]);
plt.ylim([-0.2,1.2]);
ax.set_aspect('equal');

plt.show()
