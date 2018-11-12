import numpy as np
import matplotlib.pyplot as plt
import sys

R = 1.0
testdir = sys.argv[1]
tsave   = int(sys.argv[2])
tfinal  = int(sys.argv[3])

NT    = tfinal/tsave;
theta = np.linspace(0,2*np.pi,100);
xcirc = R*np.cos(theta);
ycirc = R*np.sin(theta);

for i in range(NT):
    xy = np.genfromtxt(testdir + 'XY_' + str(tsave*(i+1)) + '.csv', delimiter=',')
    for j in range(xy.shape[0]):
        plt.plot(xcirc + xy[j,0] , ycirc + xy[j,1] , 'b');

plt.gca().set_aspect('equal');
plt.show();
