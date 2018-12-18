
import numpy as np
import matplotlib.pyplot as plt
import cgs_const as cgs
from scipy import interpolate
import InterpDiffusion
import InterpThermohaline

plt.rcParams.update({'font.serif':'Times New Roman',
                     'lines.markersize':4,
                     'lines.linewidth':1.5,
                     'text.usetex':True,
                     'font.size':10,
                     'font.family':'serif',
                     'axes.titlesize':'medium',
                     'axes.labelsize':'medium',
                     'legend.fontsize':8,
                     'legend.frameon':False,
                     'figure.dpi':300,
                     'xtick.minor.visible':True,
                     'ytick.minor.visible':True,
                     'savefig.bbox':'tight',
                     'savefig.pad_inches':0.1,
                     'savefig.dpi':300,
                     'savefig.format':'pdf',
                     'xtick.direction':'in',
                     'xtick.top':True,
                     'ytick.direction':'in',
                     'ytick.right':True,
                     'axes.formatter.use_mathtext':True,
                     'figure.autolayout':True})


# Get the data from the table of Koester & Wilken (2006)

d = np.genfromtxt('DA_data/KoesterWilken06.data')

Teff = d[:,1]
logg = d[:,2]
logCa = d[:,3]

XCa = 40.078*(10**logCa)

# Use the diffusion interpolator to infer accretion rates
mdot = 10**InterpDiffusion.log_mdots(XCa,Teff,logg)
# Use the thermohaline interpolator to infer accretion rates
mdot_th = 10**InterpThermohaline.log_mdots(XCa,Teff,logg)


# Load DB Accretion rates from Farihi et al. (2012)
jay = np.genfromtxt('DB_mdots/he.dat')
jay_teff = jay[:,0]
jay_mdot = 10**jay[:,1]


plt.figure(figsize=(3.38,3.38*3.5/4.))
ax = plt.gca()

ax.scatter(jay_teff,jay_mdot,marker='s',s=3**2,alpha=0.5,c='tab:orange',label="DBZ Farihi et al.~(2012)")
ax.scatter(Teff,mdot,label="Diffusion Only",marker='o',edgecolor='tab:gray',c='w')
ax.scatter(Teff,mdot_th,marker='+',s=6**2,label="Thermohaline")
ax.set_yscale('log')

handles,labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles[:1],labels[:1],loc=4,frameon=True,handletextpad=0.1)
ax.add_artist(first_legend)
ax.legend(handles[1:],labels[1:],loc=2,frameon=True,handletextpad=0.1)

ax.set_ylabel(r'$\dot M_{\rm acc}$ [$\rm g \, s^{-1}$]')

plt.xlabel(r'$T_{\rm eff}$ [K]')

plt.xlim(3000,21500)
plt.ylim(2e4,9e12)

# Build interpolator to map Teff to age
d = np.genfromtxt('WD_cool/cooling_0.6.data')
Teff = d[:,0]
age = d[:,1]
interp_age = interpolate.interp1d(Teff,age,fill_value='extrapolate')
interp_Teff = interpolate.interp1d(age,Teff,fill_value='extrapolate')

tickTeffs = [interp_Teff(3e9),interp_Teff(1e9),interp_Teff(3e8),interp_Teff(1e8)]
tickNames = [r'$3 \times 10^9$',r'$10^9$',r'$3 \times 10^8$',r'$10^8$']

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(tickTeffs)
ax2.set_xticklabels(tickNames)
ax2.set_xlabel(r'Estimated Cooling Time [yr] (for $0.6 \, M_\odot$ WD)',fontsize='small')
ax2.minorticks_off()

plt.savefig('Mdot.pdf')
