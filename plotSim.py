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
dt    = 0.01;
num   = 144;
R     = 0.01;
samp  = 30;
x     = R*np.cos(np.linspace(-np.pi,np.pi,samp));
y     = R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy    = np.transpose(np.array([x,y]));

xyw   = np.genfromtxt('OUT/WALL.csv',delimiter=',');
plt.plot(xyw[:,0],xyw[:,1],'k');

for i in range(0,steps):
    print('OUT/T' + str((i+1)*dt))
    plt.gca().clear();
    time = (i+1)*dt;
    try:
        xy0   = np.genfromtxt('OUT/T' + str((i+1)*dt) + '.csv',delimiter=',');
        for j in range(0,num):
            xt    = xy0[j,0] + x;
            yt    = xy0[j,1] + y;
            plt.plot(xt,yt,'b');
        xt    = xy0[-1,0] + 5*x;
        yt    = xy0[-1,1] + 5*y;
        plt.plot(xt,yt,'b');
    except:
        pass;
    plt.plot(xyw[:,0],xyw[:,1],'k');
    plt.xlim([-1,1]);
    plt.ylim([0,1.1])
    plt.gca().set_aspect('equal');
    plt.draw();
    plt.show();
    #TIME.sleep(0.5);
