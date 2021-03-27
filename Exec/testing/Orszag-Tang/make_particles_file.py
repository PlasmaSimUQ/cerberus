import numpy as np

outfile_name = "particle_file"

# number of particles
n_particles = 500

xy = np.random.random((n_particles,2))

with open(outfile_name, 'w') as outfile:
    outfile.write("{}\n".format(n_particles))

    for i in range(n_particles):
        outfile.write("{} {}\n".format(xy[i,0], xy[i,1]))
