# This example uses a MovieWriter directly to grab individual frames and
# write them to a file. This avoids any event loop integration, but has
# the advantage of working with even the Agg backend. This is not recommended
# for use in an interactive setting.
# -*- noplot -*-

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)
fig = plt.figure()
l, = plt.plot([], [], 'k-o')
plt.xlim([-1,1]);
plt.ylim([0,1.1])
x0, y0 = 0, 0

# Set up simulation parameters
steps = 1000;
dt    = 0.02;
num   = 27;
R     = 0.05;
samp  = 30;
x     = R*np.cos(np.linspace(-np.pi,np.pi,samp));
y     = R*np.sin(np.linspace(-np.pi,np.pi,samp));
xy    = np.transpose(np.array([x,y]));
xyw   = np.genfromtxt('OUT/WALL.csv',delimiter=',');


with writer.saving(fig, "Particle_2D.mp4",100):
    # First time frame
    xy0   = np.genfromtxt('OUT/T0.0001.csv',delimiter=',');
    for j in range(0,num):
        xt    = xy0[j,0] + x;
        yt    = xy0[j,1] + y;
        plt.plot(xt,yt,'b');
    plt.plot(xyw[:,0],xyw[:,1],'k');
    plt.xlim([-1,1]);
    plt.ylim([0,1.1])
    plt.gca().set_aspect('equal');
    writer.grab_frame();
    plt.gca().clear();
    # Plot different times
    for i in range(0,steps):
        time = (i+1)*dt;
        print(time)
        try:
            xy0   = np.genfromtxt('OUT/T' + str((i+1)*dt) + '.csv',delimiter=',');
            for j in range(0,num):
                xt    = xy0[j,0] + x;
                yt    = xy0[j,1] + y;
                plt.plot(xt,yt,'b');
        except:
            print('OUT/T' + str((i+1)*dt) + '.csv')
            pass;
        plt.plot(xyw[:,0],xyw[:,1],'k');
        plt.xlim([-1,1]);
        plt.ylim([0,1.1])
        plt.gca().set_aspect('equal');
        writer.grab_frame();
        plt.gca().clear();




# for i in range(100):
#         x0 += 0.1 * np.random.randn()
#         y0 += 0.1 * np.random.randn()
#         l.set_data(x0, y0)
#         writer.grab_frame()
