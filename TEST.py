import numpy as np
import matplotlib.pyplot as plt
import random

# Initialization
x1  = np.array([-1.0,1.0]);
x2  = np.array([1.0,1.0]);
m1  = 1.0;   m2  = 2.0;
v1x = 1.0;   v1y = 1.0; v1 = np.array([v1x,v1y]);
v2x = -1.0;  v2y = 0.0; v2 = np.array([v2x,v2y]);

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

steps = 1000;
eps   = 0.01;
#x     = np.ones(5);
x = np.array([-v1x,-v1y,-v2x,-v2y,0])
for i in range(0,steps):
    print(i+1);
    DF = jacobian(x1,x2,x,m1,m2,v1,v2);
    x += -eps*DF;
    f = computeF(x1,x2,x,m1,m2,v1,v2);
    #plt.scatter(i,f,c='b');

# Visualization
plt.ion()
dt    = 0.4;
steps = 40;
x01   = x1[0] - steps/2*v1x*dt; y01   = x1[1] - steps/2*v1y*dt;
x02   = x2[0] - steps/2*v2x*dt; y02   = x2[1] - steps/2*v2y*dt;
th    = np.linspace(-np.pi,np.pi,100);
R     = 0.5*np.linalg.norm(x1-x2);
cx    = R*np.cos(th);
cy    = R*np.sin(th);
for i in range(0,steps/2+1):
    plt.gca().clear()
    x1t = x01 + i*v1x*dt; y1t = y01 + i*v1y*dt;
    x2t = x02 + i*v2x*dt; y2t = y02 + i*v2y*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.xlim([-5*R,5*R]);
    plt.ylim([-5*R,5*R])
    plt.gca().set_aspect('equal');
    plt.draw()
    plt.show()

for i in range(0,steps/2+1):
    plt.gca().clear()
    x1t = x1[0] + i*x[0]*dt; y1t = x1[1] + i*x[1]*dt;
    x2t = x2[0] + i*x[2]*dt; y2t = x2[1] + i*x[3]*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.xlim([-5*R,5*R]);
    plt.ylim([-5*R,5*R])
    plt.gca().set_aspect('equal');
    plt.draw()
    plt.show()
    
print(x);


