import numpy as np

outfile_name = "particle_file"

# number of particles
n_particles = 100

y = np.linspace(0.0,1.0 - 1/2000.0,n_particles)
x = 0.1*np.cos(2*np.pi*y)

with open(outfile_name, 'w') as outfile:
    outfile.write("{}\n".format(n_particles))

    for i in range(n_particles):
        outfile.write("{} {}\n".format(x[i], y[i]))
