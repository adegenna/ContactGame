import numpy as np
import matplotlib.pyplot as plt
import sys
import time

R = 1.0
testdir = sys.argv[1]
filebase= sys.argv[2]
wallfile = sys.argv[3]
tsave   = int(sys.argv[4])
tfinal  = int(sys.argv[5])

NT    = tfinal/tsave;
theta = np.linspace(0,2*np.pi,100);
xcirc = R*np.cos(theta);
ycirc = R*np.sin(theta);

xyw = np.genfromtxt(testdir + wallfile, delimiter=',');
fig = plt.figure(1);
plt.ion();
plt.show()
for i in range(NT):
    xy = np.genfromtxt(testdir + filebase + str(tsave*(i+1)) + '.csv', delimiter=',')
    for j in range(xy.shape[0]):
        plt.plot(xcirc + xy[j,0] , ycirc + xy[j,1] , 'b');
    for j in range(xyw.shape[0]):
        plt.plot(xcirc + xyw[j,0] , ycirc + xyw[j,1] , 'r');
    plt.gca().set_aspect('equal');
    plt.xlim([-20,20])
    plt.ylim([-4,36])
    fig.canvas.draw()
    time.sleep(0.01)
    plt.clf()

plt.ioff(); plt.show();
fig = plt.figure(2);
for i in range(NT):
    em = np.genfromtxt(testdir + filebase + "EM_" + str(tsave*(i+1)) + '.csv', delimiter=',')
    plt.plot(i,em[0],'bo'); plt.plot(i,em[1],'ro'); plt.plot(i,em[2],'go');

plt.show();
