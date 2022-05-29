import h5py
import numpy
import matplotlib.pyplot as plt

h5in = h5py.File("../bin/shock_test.h5", 'r')

out_step = 2000

#open the main group
h5group = h5in['Data']

# load x and z arrays
dataset = h5group['x_pos']
x_pos = dataset[:]
dataset = h5group['z_pos']
z_pos = dataset[:]

# load datasets
# open group
groupname = "step" + str(out_step)
print(groupname)

h5group_step = h5group[groupname]
dataset = h5group_step['waves']
waves = dataset[:,:]
dataset = h5group_step['particles']
particles = dataset[:,:]
print(particles[:,2])

# Generate a related figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pos, particles[:,2], linestyle='-',color='k')
plt.ylim([0., 1.e-2])
fig.savefig('shock_test.pdf')
