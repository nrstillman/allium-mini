import pickle
import allium
import matplotlib.pyplot as plt

file = '../allium_data/sbio/4pphases_factive_1.7e-01_pairatt_7.0e-02_tau_2.3e+01_phi_9.8e-01.p'

with open(file,'rb') as f:
	obs = pickle.load(f)

print(f'pairatt = {obs.param.pairatt[0][0]}')
print(f'factive = {obs.param.factive[0]}')
print(f'tau = {obs.param.tau[0]}')
print(f'phi = {obs.param.phi}')

print(f'\nTotal number of snapshots = {obs.Nsnap}')
print(f'Total number of particles = {obs.N}')

obs.param.framerate=1

print('\nMaking animation of first 50 frames')
allium.utils.plot_sim(obs,'test_video_from_initial',time=50)

print('\nCalculating summary features')
ssvect, ssdata = allium.summstats.calculate_summary_statistics(
     obs,
     opts=['A', 'B','C','D', 'E','G'],
     log=True
     , plot=True)

print('Calculating g(r) — takes a long time')
rdist, gr = allium.summstats.calcgr(obs, resolution=0.05, verbose = True)

print('Making animation of 50 frames post-truncation')

allium.utils.plot_sim(obs,'test_video_post_trunc',time=50)
 
#plotting summary snapshots

import matplotlib.gridspec as gridspec

# Create 2x2 sub plots
gs = gridspec.GridSpec(1,3)

t = 50
fig = plt.figure(dpi=200)
ax1 = fig.add_subplot(gs[0])
ax1.scatter(obs.rval[t][:,0],obs.rval[t][:,1], s=1)
ax1.set_aspect('equal', adjustable='box')
ax1.set_axis_off()
ax1.set_title(f't={t}')

ax2 = fig.add_subplot(gs[1])
ax2.scatter(obs.rval[t+200][:,0],obs.rval[t+200][:,1], s=1)
ax2.set_aspect('equal', adjustable='box')
ax2.set_axis_off()
ax2.set_title(f't={t+200}')

ax3 = fig.add_subplot(gs[2])
ax3.scatter(obs.rval[-1][:,0],obs.rval[-1][:,1], s=1)
ax3.set_aspect('equal', adjustable='box')
ax3.set_axis_off()
ax3.set_title(f't={len(obs.rval)}')
plt.show()