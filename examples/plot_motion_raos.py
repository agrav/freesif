# -*- coding: utf-8 -*-
"""Plot motion RAOs from G1.SIF file
"""

import freesif as fs
import matplotlib.pyplot as plt

# convert data into HDF5 format
#fs.sif2hdf5('../test_files/slowdrift_G1.SIF')

# open the hdf5 file (returns a File object)
#f = fs.open_hdf5('G1.h5')

# access the HydroData object via the dict interface
#d = f['G1']

# alternative short hand method to do above 3 steps in once:
d = fs.open_sif('../tests/files/hydro/slowdrift_G1.SIF')

# get motion data
dirs = d.get_directions('degrees')
periods = d.get_periods()
motions = d.get_motion_raos()

# plot data
resp_names = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
markers = ['.-','o-','v-','^-','<-','>-','1-','2-','3-','4-','8-','s-',
           'p-','*-','h-','H-','+-','x-','D-','d-','|-','_-']

fig, axes = plt.subplots(3,2, figsize=(10,12), sharex=False)

for i, ax in enumerate(axes.flat):
    for j in range(len(dirs)):
        ax.plot(periods, abs(motions[i,j,:]), markers[j])
    ax.set_title(resp_names[i])
    ax.grid()
    ax.set_xlabel('[s]')
    ax.set_ylabel('[m/m]' if i<3 else 'rad/m')

leg = fig.legend(ax.lines, dirs, loc = 'lower center', ncol=7)
leg.set_title('Wave directions [deg]')
fig.tight_layout(rect=(0.0,0.08,1,.95))
fig.suptitle('Motion RAOs', size=20)
#fig.savefig('motions.png', dpi=200)

plt.show()

# close data
d.close()
