import random as rand
import numpy as np
import sys


rand.seed()	              # starts the Mersenne Twister rand engine

try:
    num_sources = int(sys.argv[1])
except:
    num_sources = int(2e8)           # particles on non-capped sphere

##########################
# Change parameters here #
##########################

num_BS_sources = 0            # backscattered particle number
sphereR = 30                  # cm
exclusion_angle_deg = 64      # degrees

print("Generating autorun file with %i particles..." % num_sources)

exclusion_angle = exclusion_angle_deg * np.pi / 180.

# Generate all random numbers at once
randoms = np.random.rand(num_sources, 5)
energy = 100	# keV


with open('rand_cap_sphere.mac', 'w') as f:
    # init lines to set up simulation
    f.write('/run/initialize \n')
    f.write('/control/verbose 0 \n')
    f.write('/run/verbose 0 \n')
    f.write('/event/verbose 0 \n')
    f.write('/tracking/verbose 0 \n')

    for i in range(0, num_sources):

        # Rejection sampling on TRAPPED PARTICLES
        y_pos = -sphereR 
        while y_pos < np.cos(exclusion_angle)*sphereR:  # spherical cap height
            # Uniform sphere
            theta_rand = randoms[i,0]*2*np.pi # [0, 2pi)
            u_rand     = randoms[i,1]*2-1     # [-1, 1)

            # Positions on sphere of radius R
            x_pos = sphereR * np.sqrt(1 - u_rand**2) * np.cos(theta_rand)
            y_pos = sphereR * np.sqrt(1 - u_rand**2) * np.sin(theta_rand)
            z_pos = sphereR * u_rand


        position_string = str(x_pos) + ' ' + str(y_pos) + ' ' + str(z_pos) + ' cm'

        # Draw random particle directions
        x_dir = rand.uniform(0, 1)
        y_dir = rand.uniform(0, 1)
        z_dir = rand.uniform(0, 1)

        # Enforce inward particle direction
        if x_pos > 0:
            x_dir = -x_dir
        if y_pos > 0:
            y_dir = -y_dir
        if z_pos > 0:
            z_dir = -z_dir



        f.write('/gps/source/add ' + str(i+1) + '\n')
        f.write('/gps/particle e-\n')
        f.write('/gps/energy ' + str(energy) + ' keV\n')
        f.write('/gps/position ' + str(x_pos) + ' ' + str(y_pos) + ' ' + str(z_pos) + ' cm\n')
        f.write('/gps/direction ' + str(x_dir) + ' ' + str(y_dir) + ' ' + str(z_dir) + '\n')
        f.write('/gps/pos/type Point \n')


    for i in range(0, num_BS_sources):	

        # Rejection sampling on BACKSCATTERED PARTICLES
        y_pos = sphereR
        while y_pos > np.cos(exclusion_angle)*sphereR:  # spherical cap height
            # Uniform sphere
            theta_rand = rand.uniform(0, 2*np.pi)
            u_rand   = rand.uniform(-1, 1)

            # Positions on sphere of radius R
            x_pos = sphereR * np.sqrt(1 - u_rand**2) * np.cos(theta_rand)
            y_pos = sphereR * np.sqrt(1 - u_rand**2) * np.sin(theta_rand)
            z_pos = sphereR * u_rand

	
        position_string = str(x_pos) + ' ' + str(y_pos) + ' ' + str(z_pos) + ' cm'

        # Draw random particle directions
        x_dir = rand.uniform(0, 1)
        y_dir = rand.uniform(0, 1)
        z_dir = rand.uniform(0, 1)

        # Enforce inward particle direction
        if x_pos > 0:
            x_dir = -x_dir
        if y_pos > 0:
            y_dir = -y_dir
        if z_pos > 0:
            z_dir = -z_dir


        direction_string = str(x_dir) + ' ' + str(y_dir) + ' ' + str(z_dir)

        f.write('/gps/source/add ' + str(i+1) + '\n')
        f.write('/gps/particle e-\n')
        f.write('/gps/ene/min 100. keV \n')
        f.write('/gps/ene/max 10. MeV \n')
        f.write('/gps/position ' + position_string + '\n')
        f.write('/gps/direction ' + direction_string + '\n')
        f.write('/gps/pos/type Point \n')



    f.write('\n/run/beamOn ' + str(num_sources+num_BS_sources) + '\n')

print("Autorun file generated!")
