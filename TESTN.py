import numpy as np
import matplotlib.pyplot as plt
import random

np.set_printoptions(linewidth=132)

# Initialization
x1  = np.array([0.0,0.0]);
x2  = np.array([0.0,1.0]);
x3  = np.array([0.5 , -np.sqrt(3.0)*0.5]);
x4  = np.array([-0.5, -np.sqrt(3.0)*0.5]);
m1  = 1.0;   m2  = 1.0;    m3 = 1.0;  m4 = 1.0;
v1x = 0.0;   v1y = 0.0;    v1 = np.array([v1x,v1y]);
v2x = 0.0;   v2y = -2.0;   v2 = np.array([v2x,v2y]);
v3x = 0.0;   v3y = 0.0;   v3 = np.array([v3x,v3y]);
v4x = 0.0;   v4y = 0.0;   v4 = np.array([v4x,v4y]);

# Jacobian calculation
def jacobian(x,X,V,E,M,num):
    # Calculate initial/final energies
    E0 = 0.0; Ef = 0.0;
    for i in range(0,num):
        E0 += 0.5*M[i]*(np.power(V[0,i],2) + np.power(V[1,i],2));
        Ef += 0.5*M[i]*(np.power(x[2*i],2) + np.power(x[2*i+1],2));
    # Construct mass matrix
    MASS = np.zeros([2*num,2*num]);
    for i in range(0,num):
        MASS[2*i,2*i]     = M[i];
        MASS[2*i+1,2*i+1] = M[i];
    # Construct impulse line-of-action matrix
    rows = 2*num; cols = num*(num-1)/2;
    CONT = np.zeros([rows,cols]);
    count = 0;
    for i in range(0,num):
        for j in range(i+1,num):
            K = E[i,j];
            e = np.squeeze((X[:,i]-X[:,j])/np.linalg.norm(X[:,i]-X[:,j]));
            CONT[2*i:2*i+2,count] = K*e;
            CONT[2*j:2*j+2,count] = -K*e;
            count += 1;
    CONT = CONT[:,~(CONT==0).all(0)]; # Remove zero columns from line-of-action matrix  
    # Construct matrix derivative of energy equation
    ENERGY = np.zeros(2*num);
    for i in range(0,num):
        ENERGY[2*i]   = M[i]*x[2*i];
        ENERGY[2*i+1] = M[i]*x[2*i+1];
    DF1 = np.concatenate([MASS,CONT],axis=1);
    DF2 = np.concatenate([ENERGY,np.zeros(np.shape(CONT)[1])]);
    DF  = np.vstack([DF1,DF2]);
    # Calculate function
    f = np.zeros(2*num+1);
    for i in range(0,num):
        f[2*i]   = M[i]*x[2*i]   - M[i]*V[0,i] + np.dot(CONT[2*i,:]  , x[2*num:]);
        f[2*i+1] = M[i]*x[2*i+1] - M[i]*V[1,i] + np.dot(CONT[2*i+1,:], x[2*num:]);
    f[-1] = Ef - E0;
    # Calculate DJ
    DJ   = np.dot(np.transpose(DF),f);
    fval = np.dot(np.transpose(f),f);
    return DJ,fval;


# Simulation
num   = 4;
steps = 2000;
eps   = 0.01;
x = np.zeros(11);
#x = np.array([-v1x,-v1y,-v2x,-v2y,-v3x,-v3y,-v4x,-v4y,0,0,0]);
X = np.transpose(np.concatenate([[x1],[x2],[x3],[x4]],axis=0));
V = np.transpose(np.concatenate([[v1],[v2],[v3],[v4]],axis=0));
M = np.array([m1,m2,m3,m4]);
E = np.zeros([4,4]);
E[0,1] = 1; E[0,2] = 1; E[0,3] = 1;
VELTOT = 0;
for j in range(0,num):
    VELTOT += np.linalg.norm(V[:,j]);
for i in range(0,num):
    x[2*i]   = -V[0,i]*4*np.linalg.norm(V[:,i])/VELTOT;
    x[2*i+1] = -V[1,i]*4*np.linalg.norm(V[:,i])/VELTOT;
x[2*num:] = 0;
for i in range(0,steps):
    print(i+1);
    DF,f = jacobian(x,X,V,E,M,num);
    x += -eps*DF;
    #plt.scatter(i,f,c='b');
plt.show()

# Visualization
plt.ion()
dt    = 0.05;
steps = 40;
x01   = x1[0] - steps/2*v1x*dt; y01   = x1[1] - steps/2*v1y*dt;
x02   = x2[0] - steps/2*v2x*dt; y02   = x2[1] - steps/2*v2y*dt;
x03   = x3[0] - steps/2*v3x*dt; y03   = x3[1] - steps/2*v3y*dt;
x04   = x4[0] - steps/2*v4x*dt; y04   = x4[1] - steps/2*v4y*dt;
th    = np.linspace(-np.pi,np.pi,100);
R     = 0.5*np.linalg.norm(x1-x3);
cx    = R*np.cos(th);
cy    = R*np.sin(th);
for i in range(0,steps/2+1):
    plt.gca().clear()
    x1t = x01 + i*v1x*dt; y1t = y01 + i*v1y*dt;
    x2t = x02 + i*v2x*dt; y2t = y02 + i*v2y*dt;
    x3t = x03 + i*v3x*dt; y3t = y03 + i*v3y*dt;
    x4t = x04 + i*v4x*dt; y4t = y04 + i*v4y*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.plot(x3t+cx,y3t+cy,'g');
    plt.plot(x4t+cx,y4t+cy,'m');
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
    x4t = x4[0] + i*x[6]*dt; y4t = x4[1] + i*x[7]*dt;
    plt.plot(x1t+cx,y1t+cy,'b');
    plt.plot(x2t+cx,y2t+cy,'r');
    plt.plot(x3t+cx,y3t+cy,'g');
    plt.plot(x4t+cx,y4t+cy,'m');
    plt.xlim([-5*R,5*R]);
    plt.ylim([-5*R,5*R])
    plt.gca().set_aspect('equal');
    plt.draw()
    plt.show()
    
print(x);
print(f);

