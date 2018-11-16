import numpy as np
import matplotlib.pyplot as plt

levels = 15
X      = np.zeros([levels*(levels+1)/2 + 1 , 5]);
R      = 1.0;
count  = 0;
X[:,-1] = R;
for i in range(1,levels+1):
    rackX = R*np.arange(-(i-1),i,2)
    for j in range(1,i+1):
        X[count,0] = rackX[j-1];
        X[count,1] = i*2*R*np.sqrt(3)/2.0
        count += 1;

X[-1,1] = -4*R;
X[-1,3] = 1.0;

np.savetxt('billiards.csv',X,delimiter=',');
