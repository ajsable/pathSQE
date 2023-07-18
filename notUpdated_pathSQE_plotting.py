# https://docs.mantidproject.org/nightly/plotting/index.html#plotting 
# Manually plot the path for the FeSi data


plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams.update({'font.size': 15})

fig=plt.figure(figsize=(12,5))

colormesh_pars={}
colormesh_pars['norm']=LogNorm(vmin=0.00005,vmax=0.005)
colormesh_pars['cmap']='viridis'


a=4.489
m_gamma = np.sqrt(2)
gamma_r = np.sqrt(3)
r_x = np.sqrt(2)
gs = GridSpec(1, 5,width_ratios=[1,1,m_gamma,gamma_r,r_x],wspace=0)


#fig,ax=plt.subplots(subplot_kw={'projection':'mantid'})
#slice_plot=ax.pcolormesh(mtd['slice_GAMMA_to_X'])

ax1 = plt.subplot(gs[0],projection='mantid')
ax2 = plt.subplot(gs[1],sharey=ax1,projection='mantid')
ax3 = plt.subplot(gs[2],sharey=ax1,projection='mantid')
ax4 = plt.subplot(gs[3],sharey=ax1,projection='mantid')
ax5 = plt.subplot(gs[4],sharey=ax1,projection='mantid')


pcm=ax1.pcolormesh(mtd['slice_GAMMA_to_X'], **colormesh_pars)
ax2.pcolormesh(mtd['slice_X_to_M'], **colormesh_pars)
ax3.pcolormesh(mtd['slice_M_to_GAMMA'],  **colormesh_pars)#cmap='viridis', norm=LogNorm(vmin=0.00005,vmax=0.005))
ax4.pcolormesh(mtd['slice_GAMMA_to_R'], **colormesh_pars) #cmap='viridis', norm=LogNorm(vmin=0.00005,vmax=0.005))
ax5.pcolormesh(mtd['slice_R_to_X'], **colormesh_pars)


#Adjust plotting parameters
ax1.set_ylabel('E (meV)')
ax1.set_xlabel('')
ax1.set_ylim(-10.,30.)
ax1.set_xlim(0,0.5)
ax1.set_xticks([0,0.5])
ax1.set_xticklabels(['$\Gamma$','$X$'])
ax1.tick_params(direction='in')

ax2.get_yaxis().set_visible(False)
ax2.set_xlabel('')
ax2.set_xlim(0,0.5)
ax2.set_xticks([0.5])
ax2.set_xticklabels(['$M$'])
ax2.tick_params(direction='in')

ax3.get_yaxis().set_visible(False)
ax3.set_xlabel('')
ax3.set_xlim(-0.5,0)
ax3.set_xticks([0])
ax3.set_xticklabels(['$\Gamma$'])
ax3.tick_params(direction='in')

ax4.get_yaxis().set_visible(False)
ax4.set_xlabel('')
ax4.set_xlim(0,0.5)
ax4.set_xticks([0.5])
ax4.set_xticklabels(['$R$'])
ax4.tick_params(direction='in')

ax5.get_yaxis().set_visible(False)
ax5.set_xlabel('')
ax5.set_xlim(-0.5,0)
ax5.set_xticks([0])
ax5.set_xticklabels(['$X$'])
ax5.tick_params(direction='in')

fig.colorbar(pcm, ax=[ax1,ax2,ax3,ax4,ax5])

#fig.show()
#plt.close('all')