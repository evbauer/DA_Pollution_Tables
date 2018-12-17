
import numpy as np
import cgs_const as cgs
from scipy import interpolate

# Build interpolators for cvz mass and diffusion timescales
# for Ca based on the master tables.


# Note: column 0 for Teff in these files just represents the label of the run,
# column 1 is the actual model Teff 
# at the end of the run that started at the temperature of column 0.

master = np.genfromtxt('diffusion_timescales/WD_M_0.9.data')
fm85 = interpolate.interp1d(master[:,1],master[:,3],fill_value='extrapolate') # log(M
def Mcvz85(Teff):
    return (10**fm85(Teff))*0.9*cgs.msun # in grams
ft85 = interpolate.interp1d(master[:,1],master[:,10],fill_value='extrapolate') # log(taudiff/yr) as a function of Teff
def tau85(Teff):
    return (10**ft85(Teff))*cgs.year # in years
logg85 = interpolate.interp1d(master[:,1],master[:,2],fill_value='extrapolate')
# needed because loggs vary slightly from 7.5,8.0,8.5,
# but this can be easily corrected by just interpolating in the actual logg of the model at the Teff being used.

master = np.genfromtxt('diffusion_timescales/WD_M_0.6.data')
fm80 = interpolate.interp1d(master[:,1],master[:,3],fill_value='extrapolate') # log(M
def Mcvz80(Teff):
    return (10**fm80(Teff))*0.6*cgs.msun # in grams
ft80 = interpolate.interp1d(master[:,1],master[:,10],fill_value='extrapolate') # log(taudiff/yr) as a function of Teff
def tau80(Teff):
    return (10**ft80(Teff))*cgs.year # in years
logg80 = interpolate.interp1d(master[:,1],master[:,2],fill_value='extrapolate')
master = np.genfromtxt('diffusion_timescales/WD_M_0.38.data')

fm75 = interpolate.interp1d(master[:,1],master[:,3],fill_value='extrapolate') # log(M
def Mcvz75(Teff):
    return (10**fm75(Teff))*0.38*cgs.msun # in grams
ft75 = interpolate.interp1d(master[:,1],master[:,10],fill_value='extrapolate') # log(taudiff/yr) as a function of Teff
def tau75(Teff):
    return (10**ft75(Teff))*cgs.year # in years
logg75 = interpolate.interp1d(master[:,1],master[:,2],fill_value='extrapolate')

# Rescale Ca mdot to a total based on Bulk Earth
from EarthComposition import *

def mdot038(XCa,Teff):
    return (XCa*Mcvz75(Teff)/tau75(Teff))/Earth_MassFrac_McD['Ca']
def mdot06(XCa,Teff):
    return (XCa*Mcvz80(Teff)/tau80(Teff))/Earth_MassFrac_McD['Ca']
def mdot09(XCa,Teff):
    return (XCa*Mcvz85(Teff)/tau85(Teff))/Earth_MassFrac_McD['Ca']

def accretion_rate(XCa,Teff,logg):
    mdot = interpolate.interp1d([logg75(Teff), logg80(Teff), logg85(Teff)],
                                [mdot038(XCa,Teff), mdot06(XCa,Teff), mdot09(XCa,Teff)],
                                fill_value='extrapolate')
    return mdot(logg)

