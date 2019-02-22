import random as rand
import numpy as np

rand.seed()	# starts the Mersenne Twister rand engine

num_sources = 1000             # particles on non-capped sphere
ratio = num_sources+1                     # will sample every 5th point in cap
sphereR = 30                  # cm
exclusion_angle_deg = 64      # degrees
exclusion_direction = [0,0,1] # unit direction

counter = 1

exclusion_angle = exclusion_angle_deg * np.pi / 180.


with open('rand_cap_sphere.mac', 'w') as f:
    # init lines to set up simulation
    f.write('/run/initialize \n')
    f.write('/control/verbose 0 \n')
    f.write('/run/verbose 0 \n')
    f.write('/event/verbose 0 \n')
    f.write('/tracking/verbose 0 \n')

    for i in range(0, num_sources):

        '''
        # Sampling inside a spherical cap
        z_sc   = rand.uniform(0, 1) * (1 - np.cos(exclusion_angle)) + np.cos(exclusion_angle)
        phi_sc = rand.uniform(0, 1) * 2 * np.pi
        x_sc   = np.sqrt(1 - z_sc**2)*np.cos(phi_sc)
        y_sc   = np.sqrt(1 - z_sc**2)*np.sin(phi_sc)

        position_string = str(x_sc*sphereR) + ' ' + str(y_sc*sphereR) + ' ' + str(z_sc*sphereR) + ' cm'

        # Draw random particle directions
        x_dir = rand.uniform(0, 1)
        y_dir = rand.uniform(0, 1)
        z_dir = rand.uniform(0, 1)

        # Enforce inward particle direction
        if x_sc > 0:
            x_dir = -x_dir
        if y_sc > 0:
            y_dir = -y_dir
        if z_sc > 0:
            z_dir = -z_dir

        direction_string = str(x_dir) + ' ' + str(y_dir) + ' ' + str(z_dir)

        f.write('\n# Particle source (cap) ' + str(counter) + '\n')
        f.write('/gps/source/add ' + str(i+1) + '\n')
        f.write('/gps/particle e-\n')
        f.write('/gps/ene/min 100. keV \n')
        f.write('/gps/ene/max 10. MeV \n')
        f.write('/gps/position ' + position_string + '\n')
        f.write('/gps/direction ' + direction_string + '\n')
        f.write('/gps/pos/type Point \n')

        counter += 1
        '''

        # Rejection sampling on height coordinate
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

        f.write('\n# Particle source ' + str(counter) + '\n')
        f.write('/gps/source/add ' + str(i+1) + '\n')
        f.write('/gps/particle e-\n')
        f.write('/gps/ene/min 100. keV \n')
        f.write('/gps/ene/max 10. MeV \n')
        f.write('/gps/position ' + position_string + '\n')
        f.write('/gps/direction ' + direction_string + '\n')
        f.write('/gps/pos/type Point \n')

        counter += 1

    f.write('\n/run/beamOn ' + str(num_sources) + '\n')

    f.close()
