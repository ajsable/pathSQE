import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

plt.rcParams.update({'font.size': 12.5})

fig, ax = plt.subplots(figsize=(6, 4))


data1 = np.load('L2Gamma_SED.npy')
data2_temp = np.load('Gamma2X_SED.npy')
data2 = data2_temp[::-1]
data = np.vstack((data1, data2))
print(data1.shape, data2.shape, data.shape)


# Plotting
norm  = SymLogNorm(linthresh=3e-3,vmin=3e-3,vmax=4e-2)
im = plt.imshow(data.T, norm=norm, aspect='auto', cmap='viridis', origin='lower')

#im = ax.imshow(np.log10(np.flipud(new_data_reshaped.T)), vmin=vmin, vmax=vmax, aspect='auto')
cbar = plt.colorbar(im)


# Customize ticks and labels
xpos = [-0.5, 20.5, 59.5]
xticks = ['L', '$\Gamma$', 'X']
plt.xticks(xpos, xticks)

# Set labels and title
plt.ylabel('Energy / meV', fontsize='12.5')
ypos = np.arange(0, 50, 14.7)
yticks = np.arange(0, 16, 5)
fmt = lambda x: "{:.0f}".format(x)
plt.yticks(ypos, [fmt(i) for i in yticks], fontsize='12.5')
plt.ylim([0, 50])

plt.savefig('300.png')
