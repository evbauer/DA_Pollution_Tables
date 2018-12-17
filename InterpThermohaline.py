
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

Teffs = [6000,
         7000,
         8000,
         9000,
         10000,
         11000,
         12000,
         13000,
         14000,
         15000,
         18000,
         20500]

# Interpolating functions at constant Teff
fs038 = []
fs06 = []
fs09 = []
# need to read the slightly different value from the MESA model that has cooled a bit over 100tau_diff
actual_Teffs038 = [] 
actual_Teffs06 = [] 
actual_Teffs09 = []
loggs038 = []
loggs06 = []
loggs09 = []

for entry in Teffs:
    label = str(entry) + 'K'

    filename = 'thermohaline_surfaceX/WD_M_0.38/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa038 = d[:,1]
    teff038 = d[0,3]
    logg038 = d[0,4]
    # Interpolating log XCa as a function of log Mdot
    fs038.append(interpolate.interp1d(np.log10(XCa038),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs038.append(teff038)
    loggs038.append(logg038)

    filename = 'thermohaline_surfaceX/WD_M_0.6/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa06 = d[:,1]
    diff = d[:,2] # unused, prediction of diffusion
    teff06 = d[0,3]
    logg06 = d[0,4]
    # Interpolating log XCa as a function of log Mdot
    fs06.append(interpolate.interp1d(np.log10(XCa06),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs06.append(teff06)
    loggs06.append(logg06)

    filename = 'thermohaline_surfaceX/WD_M_0.9/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa09 = d[:,1]
    teff09 = d[0,3]
    logg09 = d[0,4]
    # Interpolating log XCa as a function of log Mdot
    fs09.append(interpolate.interp1d(np.log10(XCa09),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs09.append(teff09)
    loggs09.append(logg09)


# Make function for logg as a function of Teff for each mass
g038 = interpolate.interp1d(actual_Teffs038,loggs038,fill_value='extrapolate')
g06 = interpolate.interp1d(actual_Teffs06,loggs06,fill_value='extrapolate')
g09 = interpolate.interp1d(actual_Teffs09,loggs09,fill_value='extrapolate')


# Interpolate to find mdot assuming fixed mass for each mass.
# Builds a temporary interpolator at fixed log_XCa each time.
# Inefficient, but seems robust, and I don't need it to be fast.

def log_mdot038(log_XCa,Teff):
    log_mdots_for_this_XCa = []
    for logmdot in fs038:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs038,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)

def log_mdot06(log_XCa,Teff):
    log_mdots_for_this_XCa = []
    for logmdot in fs06:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs06,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)

def log_mdot09(log_XCa,Teff):
    log_mdots_for_this_XCa = []
    for logmdot in fs09:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs09,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)


def log_mdot(XCa,Teff,logg):
    log_XCa = np.log10(XCa)
    # Interpolate the other three functions in logg at this Teff
    f = interpolate.interp1d([g038(Teff), g06(Teff), g09(Teff)],
                             [log_mdot038(log_XCa,Teff), log_mdot06(log_XCa,Teff), log_mdot09(log_XCa,Teff)],
                             fill_value='extrapolate')
    return f(logg)

