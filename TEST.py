import numpy as np
import matplotlib.pyplot as plt
import random

# Initialization
x1  = np.array([0.0,0.0]);
x2  = np.array([1.0,0.0]);
m1  = 1.0;
m2  = 1.0;
v1x = 1.0;  v1y = 1.0; v1 = np.array([v1x,v1y]);
v2x = -np.sqrt(2.0); v2y = 0.0; v2 = np.array([v2x,v2y]);

# Jacobian calculation
def jacobian(x1,x2,x,m1,m2,v01,v02):
    e  = (x2-x1)/np.linalg.norm(x2-x1);
    ex = e[0];  ey = e[1];
    s1 = x[0]; s2  = x[1]; s3 = x[2]; s4 = x[3]; s5 = x[4];
    E0 = 0.5*m1*np.power(np.linalg.norm(v01),2) + 0.5*m2*np.power(np.linalg.norm(v02),2);
    Ef = 0.5*m1*(np.power(s1,2) + np.power(s2,2)) + 0.5*m2*(np.power(s3,2) + np.power(s4,2));
    DF = np.array([[m1    ,  0    ,  0    ,  0    ,  ex],
                   [0     ,  m1   ,  0    ,  0    ,  ey],
                   [0     ,  0    ,  m2   ,  0    , -ex],
                   [0     ,  0    ,  0    ,  m2   , -ey],
                   [m1*s1 , m1*s2 , m2*s3 , m2*s4 , 0  ]]);
    f = np.array([m1*s1 - m1*v01[0] + s5*ex,
                  m1*s2 - m1*v01[1] + s5*ey,
                  m2*s3 - m2*v02[0] - s5*ex,
                  m2*s4 - m2*v02[1] - s5*ey,
                  Ef - E0]);
    DJ = np.dot(np.transpose(DF),f);
    return DJ;

def computeF(x1,x2,x,m1,m2,v01,v02):
    e  = (x2-x1)/np.linalg.norm(x2-x1);
    ex = e[0];  ey = e[1];
    s1 = x[0]; s2  = x[1]; s3 = x[2]; s4 = x[3]; s5 = x[4];
    E0 = 0.5*m1*np.power(np.linalg.norm(v01),2) + 0.5*m2*np.power(np.linalg.norm(v02),2);
    Ef = 0.5*m1*(np.power(s1,2) + np.power(s2,2)) + 0.5*m2*(np.power(s3,2) + np.power(s4,2));
    f = np.array([m1*s1 - m1*v01[0] + s5*ex,
                  m1*s2 - m1*v01[1] + s5*ey,
                  m2*s3 - m2*v02[0] - s5*ex,
                  m2*s4 - m2*v02[1] - s5*ey,
                  Ef - E0]);
    return(np.dot(np.transpose(f),f));

steps = 3000;
eps   = 0.01;
x     = np.ones(5);
for i in range(0,steps):
    print(i+1);
    DF = jacobian(x1,x2,x,m1,m2,v1,v2);
    x += -eps*DF;
    f = computeF(x1,x2,x,m1,m2,v1,v2);
    plt.scatter(i,f,c='b');

print(x);

plt.show();
