import h5py
import numpy
import matplotlib.pyplot as plt

def load_step(h5group, i_step) :

    # load datasets
    # open group
    groupname = "step" + str(out_step)
    print(groupname)

    h5group_step = h5group[groupname]
    dataset = h5group_step['waves']
    waves = dataset[:,:]
    dataset = h5group_step['particles']
    particles = dataset[:,:]

    return waves, particles


h5in = h5py.File("../bin/shock_test.h5", 'r')

out_step = 950000

#open the main group
h5group = h5in['Data']

# load x and z arrays
dataset = h5group['x_pos']
x_pos = dataset[:]
dataset = h5group['z_pos']
z_pos = dataset[:]

# Load data
waves, particles = load_step(h5group, out_step)

print(particles[:,2])

# Generate a related figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pos, particles[:,2], linestyle='-',color='k')
plt.ylim([0., 5.e-3])
fig.savefig('shock_test.pdf')

# Generate a related figure
powerlaw = 8.e-3*z_pos**-2.
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(z_pos, particles[50,:], linestyle='-',color='k')
ax.plot(z_pos, powerlaw, linestyle=':',color='r')
plt.xlim([1., 1.e2])
plt.ylim([1.e-7, 1.e-2])
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig('spec_test.pdf')
