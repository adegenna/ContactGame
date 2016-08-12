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

steps = 1000;
dt    = 0.001;
num   = 10;

#xyw   = np.genfromtxt('OUT/WALL.csv',delimiter=',');
#plt.plot(xyw[:,0],xyw[:,1],'k');

for i in range(0,steps):
    print(i+1)
    plt.gca().clear();
    time = (i+1)*dt;
    for j in range(0,num):
        try:
            xy   = np.genfromtxt('OUT/T' + str((i+1)*dt) + '_' + str(j) + '.csv',delimiter=',');
            plt.plot(xy[:,0],xy[:,1],'b');
        except:
            pass;
    #plt.plot(xyw[:,0],xyw[:,1],'k');
    plt.xlim([-0.5,0.5]);
    plt.ylim([-0.5,0.5])
    plt.gca().set_aspect('equal');
    plt.draw();
    plt.show();
    TIME.sleep(0.0001);
