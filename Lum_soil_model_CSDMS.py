# -*- coding: utf-8 -*-
"""

This code generates the luminescence versus soil depth using an equation derived from the Fokker-Plank Equation. 
Copyright (C) 2020 Harrison Gray
Developer can be contacted by hgray@usgs.gov and at 1 Denver Federal Center, MS 963, BLDG 95, Lakewood, CO, 80225
Please let me know if you have any problems!
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


"""

import numpy as np
import matplotlib.pyplot as plt

# Set up the arrays for the for loops.

L_line_styles = ['b:','b--','b-.','b-.','b-']
De_line_styles = ['r:','r--','r-.','r-.','r-']
Pe_inv_array =  [100.0, 10.0, 1., 0.1, 0.01] # note that we are using the inverse of Pe in the paper for convenience

# First, for my own piece of mind, we estimate a ballpark Dn value by 
# using order-of-magnitude values for the various parameters used in
# its calculation. 

DR = 3.0 # Dose Rate generic value [Gy/kyr]
D0 = 50. # Characteristic Dose for quartz OSL [Gy]
gamma_p = 2.4 # dimensionless scaling constant between arithmetic and geometric means
h = 1. # Soil thickness [m]
w_0 = (114.0/1000)*np.exp(-h/0.488) # globally averaged maximum soil production [m/kyr] rate from Stockmann et al 2014; note actual values do not matter for Peclet Number approach. This is used to roughly estimate the Dn value.

Dn = gamma_p*np.round((DR/D0)*(h/w_0),decimals=1) # Non-dimensional luminescence growth parameter
print("Dn = " + str(Dn))
        
# So an average Dn lies around 10 for a thin soil. 
# I use 10 for the analysis in the paper. Using higher or lower values changes 
# the scaling, but not the overall form of the results.

Dn = 10.0

##Set parameters for non dimensional numbers

for profile_type in [0, 1, 2]: # 0 is uniform, 1 is linear, 2 is exponential
    
    for i_Pe_inv in [0, 1, 2, 3, 4]:
        
        Pe_inv = Pe_inv_array[i_Pe_inv]
        print("Pe_inv = " + str(Pe_inv))

        ## Set up array variables
        
        dz = 0.01
        z = np.arange(0.,1.,dz)
        
        L = np.linspace(1.,0.,len(z)) # Start with a linear distribution for stability
        L_old = 0.
                
        dt = 0.5*(dz**2.)/(2.*Pe_inv) # rough estimate of a stable time step
        t = 0.
        
        E = 10.0**-8
        plot_counter = 0
        
        if profile_type == 0: # uniform mixing, uniform velocity "surface erosion"
            
            print("Profile Type = "+str(profile_type))
            
            while np.abs(np.sum(L - L_old)) > E:
              
              L_old = 1.0*L  
              t += dt
              plot_counter += 1
              
              dLdz = np.diff(L)/dz
              d2Ldz2 = np.diff(dLdz)/dz 
            
              # advection  
              L[1:] += -dLdz*dt
              
              # diffusion
              L[1:-1] += Pe_inv*d2Ldz2*dt 
              
              # regen
              L[1:] += Dn*(1 - L[1:])*dt
              
              ## soil surface boundary condition
              
              L[-1] = 0
              
              ## soil / saprolite boundary condition
              
              boundary_advection = -(L[0] - 1)/dz
              boundary_diffusion = Pe_inv*((L[1] - 1*L[0] + 0)/(dz**2))
              boundary_regen = Dn*(1 - L[0])
              
              dLdt_boundary = boundary_advection + boundary_diffusion + boundary_regen;
              
              L[0] += dLdt_boundary*dt
            
            De = -np.log(1 - L)
            
            De_plot = De[np.where(De<=3.)]
            z_plot = z[np.where(De<=3.)]
            
            De_scaled = (De_plot - np.min(De_plot))/(np.max(De_plot) - np.min(De_plot))
            z_scaled = (z_plot - np.min(z_plot))/(np.max(z_plot) - np.min(z_plot))
            
        elif profile_type == 1: # uniform velocity, linear mixing, "surface erosion"
            
            print("Profile Type = "+str(profile_type))
            
            while np.abs(np.sum(L - L_old)) > E:

                L_old = 1.0*L  
                t += dt
                plot_counter += 1
                
                dLdz = np.diff(L)/dz
                d2Ldz2 = np.diff(dLdz)/dz 
                
                # advection
                L[1:] += -dLdz*dt
                
                # diffusion
                L[1:-1] += Pe_inv*(d2Ldz2 + dLdz[1:])*dt 
                
                # regen
                L[1:] += Dn*(1 - L[1:])*dt
                
                ## soil surface boundary condition
                L[-1] = 0
                
                ## soil / saprolite boundary condition
                
                boundary_advection = -(L[0] - 1)/dz
                boundary_diffusion = Pe_inv*((L[1] - 1*L[0] + 0)/(dz**2.) - (L[0] - 1)/dz)*z[0]*dt 
                boundary_regen = Dn*(1. - L[0])
                
                dLdt_boundary = boundary_advection + boundary_diffusion + boundary_regen
                
                L[0] += dLdt_boundary*dt
  
            De = -np.log(1. - L)
            
            De_plot = De[np.where(De<=3.)]
            z_plot = z[np.where(De<=3.)]
            
            De_scaled = (De_plot - np.min(De_plot))/(np.max(De_plot) - np.min(De_plot))
            z_scaled = (z_plot - np.min(z_plot))/(np.max(z_plot) - np.min(z_plot))

        elif profile_type == 2: # linear velocity / linear mixing, no surface erosion
            
            print("Profile Type = "+str(profile_type))
            
            while np.abs(np.sum(L - L_old)) > E:
                                
                L_old = 1.0*L
                t += dt
                plot_counter += 1
                
                dLdz = np.diff(L)/dz
                d2Ldz2 = np.diff(dLdz)/dz 
                
                # advection
                L[1:] += -(1 - z[1:]**2.0)*dLdz*dt
                
                # diffusion
                L[1:-1] += Pe_inv*(d2Ldz2 + dLdz[1:])*dt 
                
                # regen
                L[1:] += Dn*(1. - L[1:])*dt
                
                ## soil surface boundary condition
                L[-1] = 0
                
                ## soil / saprolite boundary condition
                boundary_advection = -(L[0] - 1.)/dz
                boundary_diffusion = Pe_inv*((L[1] - 1*L[0] + 0)/(dz**2.) - (L[0] - 1)/dz)*z[0]*dt
                boundary_regen = Dn*(1. - L[0])
                
                dLdt_boundary = boundary_regen + boundary_advection + boundary_diffusion
                
                L[0] += dLdt_boundary*dt
                
                if np.any(L<0):
                    print("negative error")
                    break
            
            De = -np.log(1 - L)
            
            De_plot = De[np.where(De<=3.)]
            z_plot = z[np.where(De<=3.)]
            
            De_scaled = (De_plot - np.min(De_plot))/(np.max(De_plot) - np.min(De_plot))
            z_scaled = (z_plot - np.min(z_plot))/(np.max(z_plot) - np.min(z_plot))
            
            
        plt.subplot(1,2,1)
        plt.plot(L, z, L_line_styles[i_Pe_inv])
        plt.plot([0,0],[0,1],'k:')
        plt.xlim(-0.1,1.1)
        plt.ylim(0,1)
        plt.title("Luminescence vs. height")
        plt.ylabel("ND Soil height")
        plt.xlabel("ND Luminescence")
        
        plt.subplot(1,2,2)
        plt.plot(De_scaled, z_scaled,De_line_styles[i_Pe_inv])
        plt.plot([0,0],[0,1],'k:')
        plt.xlim(-0.1,1.1)
        plt.ylim(0,1)
        plt.title("Normalized ND De vs. height")
        plt.xlabel("ND Equivalent Dose")
        
#        np.savetxt(("L_" + str(profile_type) + "_" + str(i_Pe_inv) + ".txt"), L, fmt="%.5f")
#        np.savetxt(("De_noscale_" + str(profile_type) + "_" + str(i_Pe_inv) + ".txt"), De, fmt="%.5f")
#        np.savetxt(("De_" + str(profile_type) + "_" + str(i_Pe_inv) + ".txt"), De_scaled, fmt="%.5f")
#        np.savetxt(("z_" + str(profile_type) + "_" + str(i_Pe_inv) + ".txt"), z_scaled, fmt="%.5f")
        
    plt.show()
        



