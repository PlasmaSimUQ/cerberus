import numpy as np

outfile_name = "particle_file"

# number of particles
n_particles = 1000

x = np.random.random(n_particles)
y = np.random.random(n_particles)

with open(outfile_name, 'w') as outfile:
    outfile.write("{}\n".format(n_particles))

    for i in range(n_particles):
        outfile.write("{} {}\n".format(x[i], y[i]))
