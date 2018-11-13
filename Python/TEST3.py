import numpy as np
import matplotlib.pyplot as plt
import random

# Initialization
x1  = np.array([0.0,0.0]);
x2  = np.array([-1.0,1.0]);
x3  = np.array([1.0,1.0]);
m1  = 1.0;   m2  = 1.0;  m3 = 1.0;
v1x = 3.0;   v1y = 1.0;   v1 = np.array([v1x,v1y]);
v2x = 5.0;   v2y = -5.0;  v2 = np.array([v2x,v2y]);
v3x = -1.0;  v3y = -1.0;  v3 = np.array([v3x,v3y]);

# Jacobian calculation
def jacobian(x1,x2,x3,x,m1,m2,m3,v01,v02,v03):
    e12  = (x2-x1)/np.linalg.norm(x2-x1);
    e13  = (x3-x1)/np.linalg.norm(x3-x1);
    s1 = x[0]; s2  = x[1]; s3 = x[2]; s4 = x[3]; s5 = x[4]; s6 = x[5]; s7 = x[6]; s8 = x[7];
    E0 = 0.5*m1*np.power(np.linalg.norm(v01),2) + 0.5*m2*np.power(np.linalg.norm(v02),2) + 0.5*m3*np.power(np.linalg.norm(v03),2);
    Ef = 0.5*m1*(np.power(s1,2) + np.power(s2,2)) + 0.5*m2*(np.power(s3,2) + np.power(s4,2)) + 0.5*m3*(np.power(s5,2) + np.power(s6,2));
    DF = np.array([[m1    ,  0    ,  0    ,  0    ,  0   ,  0   ,  e12[0]   ,  e13[0]],
                   [0     ,  m1   ,  0    ,  0    ,  0   ,  0   ,  e12[1]   ,  e13[1]],
                   [0     ,  0    ,  m2   ,  0    ,  0   ,  0   , -e12[0]   ,  0     ],
                   [0     ,  0    ,  0    ,  m2   ,  0   ,  0   , -e12[1]   ,  0     ],
                   [0     ,  0    ,  0    ,  0    ,  m3  ,  0   ,  0        , -e13[0]],
                   [0     ,  0    ,  0    ,  0    ,  0   ,  m3  ,  0        , -e13[1]],
                   [m1*s1 , m1*s2 , m2*s3 , m2*s4 , m3*s5, m3*s6,  0        ,  0     ]]);
    f = np.array([m1*s1 - m1*v01[0] + s7*e12[0] + s8*e13[0],
                  m1*s2 - m1*v01[1] + s7*e12[1] + s8*e13[1],
                  m2*s3 - m2*v02[0] - s7*e12[0]            ,
                  m2*s4 - m2*v02[1] - s7*e12[1]            ,
                  m3*s5 - m2*v03[0] -             s8*e13[0],
                  m3*s6 - m2*v03[1] -             s8*e13[1],
                  Ef - E0]);
    DJ = np.dot(np.transpose(DF),f);
    return DJ;

def computeF(x1,x2,x3,x,m1,m2,m3,v01,v02,v03):
    e12  = (x2-x1)/np.linalg.norm(x2-x1);
    e13  = (x3-x1)/np.linalg.norm(x3-x1);
    s1 = x[0]; s2  = x[1]; s3 = x[2]; s4 = x[3]; s5 = x[4]; s6 = x[5]; s7 = x[6]; s8 = x[7];
    E0 = 0.5*m1*np.power(np.linalg.norm(v01),2) + 0.5*m2*np.power(np.linalg.norm(v02),2) + 0.5*m3*np.power(np.linalg.norm(v03),2);
    Ef = 0.5*m1*(np.power(s1,2) + np.power(s2,2)) + 0.5*m2*(np.power(s3,2) + np.power(s4,2)) + 0.5*m3*(np.power(s5,2) + np.power(s6,2));
    f = np.array([m1*s1 - m1*v01[0] + s7*e12[0] + s8*e13[0],
                  m1*s2 - m1*v01[1] + s7*e12[1] + s8*e13[1],
                  m2*s3 - m2*v02[0] - s7*e12[0]            ,
                  m2*s4 - m2*v02[1] - s7*e12[1]            ,
                  m3*s5 - m2*v03[0] -             s8*e13[0],
                  m3*s6 - m2*v03[1] -             s8*e13[1],
                  Ef - E0]);
    return(np.dot(np.transpose(f),f));

steps = 2000;
eps   = 0.01;
#x     = np.ones(5);
x = np.array([-v1x,-v1y,-v2x,-v2y,-v3x,-v3y,0,0])
for i in range(0,steps):
    print(i+1);
    DF = jacobian(x1,x2,x3,x,m1,m2,m3,v1,v2,v3);
    x += -eps*DF;
    f = computeF(x1,x2,x3,x,m1,m2,m3,v1,v2,v3);
    #plt.scatter(i,f,c='b');
#plt.show()

# Visualization
plt.ion()
dt    = 0.01;
steps = 100;
x01   = x1[0] - steps/2*v1x*dt; y01   = x1[1] - steps/2*v1y*dt;
x02   = x2[0] - steps/2*v2x*dt; y02   = x2[1] - steps/2*v2y*dt;
x03   = x3[0] - steps/2*v3x*dt; y03   = x3[1] - steps/2*v3y*dt;
th    = np.linspace(-np.pi,np.pi,100);
R     = 0.5*np.linalg.norm(x1-x3);
cx    = R*np.cos(th);
cy    = R*np.sin(th);
for i in range(0,steps/2+1):
    plt.gca().clear()
    x1t = x01 + i*v1x*dt; y1t = y01 + i*v1y*dt;
    x2t = x02 + i*v2x*dt; y2t = y02 + i*v2y*dt;
    x3t = x03 + i*v3x*dt; y3t = y03 + i*v3y*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.plot(x3t+cx,y3t+cy,'g');
    plt.xlim([-5*R,5*R]);
    plt.ylim([-5*R,5*R])
    plt.gca().set_aspect('equal');
    plt.draw()
    plt.show()

for i in range(0,steps/2+1):
    plt.gca().clear()
    x1t = x1[0] + i*x[0]*dt; y1t = x1[1] + i*x[1]*dt;
    x2t = x2[0] + i*x[2]*dt; y2t = x2[1] + i*x[3]*dt;
    x3t = x3[0] + i*x[4]*dt; y3t = x3[1] + i*x[5]*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.plot(x3t+cx,y3t+cy,'g');
    plt.xlim([-5*R,5*R]);
    plt.ylim([-5*R,5*R])
    plt.gca().set_aspect('equal');
    plt.draw()
    plt.show()
    
print(x);
print(f);

