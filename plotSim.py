from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import random
import Quadtree as QT
import Surface as Surf
from contact import contact
import physics
import time as TIME

plt.ion()

steps = 200;
dt    = 0.01;

xyw   = np.genfromtxt('OUT/WALL.csv',delimiter=',');
plt.plot(xyw[:,0],xyw[:,1],'k');

for i in range(0,steps):
    print(i+1)
    plt.gca().clear();
    time = (i+1)*dt;
    xy   = np.genfromtxt('OUT/T' + str((i+1)*dt) + '.csv',delimiter=',');
    plt.plot(xyw[:,0],xyw[:,1],'k');
    plt.plot(xy[:,0],xy[:,1],'b');
    plt.xlim([-1.2,1.2]);
    plt.ylim([-0.2,1.2]);
    plt.gca().set_aspect('equal');
    plt.draw();
    plt.show();
    TIME.sleep(0.1);
