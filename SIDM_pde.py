import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--output',type = str, default = None,
    help = 'names (prefix) of output files to write')
parser.add_argument('--r_len', type = int, default=50,
    help = 'size of array along r direction')
parser.add_argument('--z_len', type = int, default=50,
    help = 'size of array along r direction')
parser.add_argument('--scale', type = float, default=-4,
    help = 'exponential term to describe innermost point r_inner = r1 * 10^scale')
parser.add_argument('--r1', type = float, default=5.,
    help = 'the boundary of the array in kpc')
parser.add_argument('--rho0', type = float, default=5.e6,
    help = 'scale density of the solution in M_sun/kpc')
parser.add_argument('--num_iter', type = int, default=1000000,
    help = 'number of iterations in the relaxation')
parser.add_argument('--no_baryons', action='store_true',
    help = 'type --no_baryons if you want an isothermal core')
args = parser.parse_args()


def interp(grid):
    grid_shape = grid.shape
    new_grid_shape = grid_shape*2
    return 5
    
def restr(grid):
    return 5

def baryon_profile(r, z):
    # r is a 2-d array of the r-values used in the grid of points
    # z is a 2-d array of the z-values used in the grid of points     
    zthin = 0.3 #in kpc
    zthick = 0.9 #in kpc
    rthin = 2.9 #in kpc
    rthick = 3.3 #in kpc
    rcut = 2.1 #in kpc
    r0 = 0.075 #in kpc
    rhobulge = 9.9e10 #in M_sun / kpc^3
    sthin = 8.17e8 #in M_sun / kpc^2
    sthick = 2.09e9 #in M_sun / kpc^2    
    q = 0.5     
    a = 1.8
    rho_baryon = rhobulge*np.exp(-(r**2 + z**2/q**2)/rcut**2) / (1 + np.sqrt(r**2 + z**2/q**2)/r0)**a + 0.5*sthin*np.exp(-z/zthin-r/rthin)/zthin + 0.5*sthick*np.exp(-z/zthick-r/rthick)/zthick
    return rho_baryon                

def NFW_template(r,z):
    rhos = 4.589e6 #in Msun/kpc^3
    rho0 = 5.e6 #in Msun/kpc^3
    Rs = 23. #in kpc
    R = np.sqrt(r**2+z**2)
    return np.log(rhos/rho0 * Rs/R * 1/(1+R/Rs)**2)

       
def relax(u, r, z, r1, rho0, iterations):
    #u is the solution to the diffeq, given in terms of log(rho/rho0)
    g= 2.71*(rho0/5.e6)*(r1/10.)**2 # 4pi * GN * rho0 *r1**2/sigma**2 for sigma 100 km/s, it is dimensionless  
    delta_r = np.gradient(r/r1)[0]  # for a 2-d array, gradient returns gradient as a list 'vector' defined at each point, 0 axis for r and 1 axis for z   
    delta_z = np.gradient(z/r1)[1]
    #time_step = 1/10000.
    for i in range(iterations):
        u[0,:] = u[1,:] 
        u[:,0] = u[:,1] 
        #enforces the derivative goes to 0 at r,z = 0, which is required by symmetry / smoothness
        residual = r1/r*np.gradient(u)[0]/delta_r + np.gradient(u,2)[0]/delta_r**2 + np.gradient(u,2)[1]/delta_z**2 + g*(baryon_profile(r,z)/rho0+np.exp(u))
        time_step = 1e-3/residual.max()
        u += time_step*residual
        u[-1,:] = NFW_template(r[-1,:],z[-1,:])
        u[:,-1] = NFW_template(r[:,-1],z[:,-1]) #used to enforce 'NFW for large R' boundary condition 
    return u, residual

def relax_no_baryons(u, r, z, r1, rho0, iterations):
    #u is the solution to the diffeq, given in terms of log(rho/rho0)
    g= 2.71*(rho0/5.e6)*(r1/10.)**2 # 4pi * GN * rho0 *r1**2/sigma**2 for sigma 100 km/s, it is dimensionless  
    delta_r = np.gradient(r/r1)[0]  # for a 2-d array, gradient returns gradient as a list 'vector' defined at each point, 0 axis for r and 1 axis for z   
    delta_z = np.gradient(z/r1)[1]
    #time_step = 1/10000.
    for i in range(iterations):
        u[0,:] = u[1,:] 
        u[:,0] = u[:,1] 
        #enforces the derivative goes to 0 at r,z = 0, which is required by symmetry / smoothness
        residual = r1/r*np.gradient(u)[0]/delta_r + np.gradient(u,2)[0]/delta_r**2 + np.gradient(u,2)[1]/delta_z**2 + g*(np.exp(u))
        time_step = 1e-3/residual.max()
        u += time_step*residual
        u[-1,:] = NFW_template(r[-1,:],z[-1,:])
        u[:,-1] = NFW_template(r[:,-1],z[:,-1]) #used to enforce 'NFW for large R' boundary condition 
    return u, residual



#r_len = 51
#z_len = 50
#scale = -4. #exponential term to describe innermost point r_inner = r1 * 10^scale
#r1 = 5. # in kpc
#rho0 = 5.e6 # scale density of the solution in M_sun/kpc^3
#num_iter = 4000000
r_len = args.r_len
z_len = args.z_len
scale = args.scale
r1 = args.r1
rho0 = args.rho0
num_iter = args.num_iter
output = args.output

f1=open(output+'_params.txt', 'w')
f1.write('len in r direction is '+str(r_len))
f1.write('len in z direction is '+str(z_len))
f1.write('magnitude of innermost point is '+str(scale))
f1.write('outer boundary is '+str(r1)+' kpc')
f1.write('number of iterations is '+str(num_iter))
f1.write('the no baryon case is '+str(args.no_baryons))
f1.close()


#arrays (r,z,u) have the shape (r_len,z_len)  
  
r = r1*np.logspace(scale,0.,r_len,endpoint=True) #the polar radius; varies from r1*( e^-scale to 1), the [::-1] is there to reverse the array to make the index start at the inner part of the halo
r = np.tile(r[:,np.newaxis],(1,z_len))

z = r1*np.logspace(scale,0.,z_len,endpoint=True) #the height; varies from r1*( e^-scale to 1), the [::-1] is there to reverse the array to make the index start at the inner part of the halo
z = np.tile(z[np.newaxis,:],(r_len,1))


u = NFW_template(r, z) #u is the solution to the diffeq, given in terms of log(rho/rho0); this line sets up the initial condition for the pde


#### PLOTS BEFORE RELAXATION ####
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(np.log(r), np.log(z), u, rstride=2, cstride=2)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'wireplot_before.png')
plt.clf()

CS = plt.contour(np.log(r),np.log(z),u)
plt.clabel(CS, inline=1,fontsize=10)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'log_density_before.png')
plt.clf()

CS = plt.contour(r,z,u)
plt.clabel(CS, inline=1,fontsize=10)
plt.xlabel('Radius (kpc)')
plt.ylabel('Height (z) (kpc)')
plt.savefig(output + 'density_before.png')
plt.clf()

plt.contourf(r,z,u)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'log_density_filled_before.png')
plt.clf()




#### RELAXATION ####
if args.no_baryons:
    u, residual = relax_no_baryons(u, r, z, r1, rho0, num_iter)
else:
    u, residual = relax(u, r, z, r1, rho0, num_iter)


#### PLOTS AFTER RELAXATION ####
CS = plt.contour(np.log(r),np.log(z),u)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'log_density_after.png')
plt.clf()

CS = plt.contour(r,z,u)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel('Radius (kpc)')
plt.ylabel('Height (z) (kpc)')
plt.savefig(output + 'density_after.png')
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(np.log(r), np.log(z), u, rstride=2, cstride=2)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'output + wireplot_after.png')
plt.clf()

plt.contourf(r,z,u)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'log_density_filled_after.png')
plt.clf()



CS = plt.contour(np.log(r),np.log(z),residual)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'final_residual.png')
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(np.log(r), np.log(z), residual, rstride=2, cstride=2)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'wireplot_residual.png')
plt.clf()

CS = plt.contour(np.log(r),np.log(z),residual/u)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'final_residual_contrast.png')
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(np.log(r), np.log(z), residual/u, rstride=2, cstride=2)
plt.xlabel(r'$\log(r)$')
plt.ylabel(r'$\log(z)$')
plt.savefig(output + 'wireplot_residual_contrast.png')
plt.clf()
                
