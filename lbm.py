import os
import numpy as np
import pandas as pd
import lbm_step
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


# loading simulation object
path = os.getcwd()
sim_object = pd.read_csv(path+'/geometries/sphere.csv').values.astype(bool)


# initalizing parameters
height = sim_object.shape[0]  # lattice y length 
width = sim_object.shape[1]   # lattice x length
u0 = 0.20                     # flow velocity
interval = 10                 # timesteps between updating animation
omega = 1.5                   # relaxation constant    


# intializing number density
w = [1/36., 1/9., 1/36., 1/9., 4/9., 1/9., 1/36., 1/9., 1/36.] # probability weights

f = np.array(w * width * height).reshape(height, width, 9) # initial number densities

f[:,:,0] *= (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)  # updating for flow speed 
f[:,:,1] *= (1 - 1.5*u0**2)
f[:,:,2] *= (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
f[:,:,3] *= (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
f[:,:,4] *= (1 - 1.5*u0**2)
f[:,:,5] *= (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
f[:,:,6] *= (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
f[:,:,7] *= (1 - 1.5*u0**2)
f[:,:,8] *= (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)


# initalizing macroscopic quantities
rho = np.ones((height, width))
ux = np.array((f[:,:,2] + f[:,:,5] + f[:,:,8] - (f[:,:,0] + f[:,:,3] + f[:,:,6])) / rho)
uy = np.array((f[:,:,0] + f[:,:,1] + f[:,:,2] - (f[:,:,6] + f[:,:,7] + f[:,:,8])) / rho)
u = np.sqrt(ux**2 + uy**2)


# animation
def animation() :
    fig = plt.figure(figsize=(10,5))
    ax = fig.gca()
    surf = ax.imshow(u, cmap = 'jet', interpolation='none')
    plt.xticks([])
    plt.yticks([])
    
    def update_data(i) :
        global f 
        ax.clear()
        plt.xticks([])
        plt.yticks([])
        for j in range(interval) : 
            f, rho, ux, uy, u = lbm_step.step(f, height, width, omega, u0, w, sim_object)
        u[sim_object] = u.mean()   # setting object to visible color
        surf = ax.imshow(u, cmap = 'jet', interpolation='none')
        
        if i%10000 == 0 :
            # need to save a figure for the animation to work for some 
            plt.savefig('placeholder.png')   
        return surf,
    
    anim = manimation.FuncAnimation(fig, update_data, interval = 1, blit = True, repeat = False)
    plt.show()
       
animation()
