
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

for entry in Teffs:
    label = str(entry) + 'K'

    filename = 'thermohaline_surfaceX/WD_M_0.38/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa038 = d[:,1]
    diff = d[:,2] # unused, prediction of diffusion
    teff038 = d[0,3]
    # Interpolating log XCa as a function of log Mdot
    fs038.append(interpolate.interp1d(np.log10(XCa038),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs038.append(teff038)

    filename = 'thermohaline_surfaceX/WD_M_0.6/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa06 = d[:,1]
    diff = d[:,2] # unused, prediction of diffusion
    teff06 = d[0,3]
    # Interpolating log XCa as a function of log Mdot
    fs06.append(interpolate.interp1d(np.log10(XCa06),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs06.append(teff06)

    filename = 'thermohaline_surfaceX/WD_M_0.9/' + label + '.data'
    d = np.genfromtxt(filename)
    mdot = d[:,0]
    XCa09 = d[:,1]
    diff = d[:,2] # unused, prediction of diffusion
    teff09 = d[0,3]
    # Interpolating log XCa as a function of log Mdot
    fs09.append(interpolate.interp1d(np.log10(XCa09),np.log10(mdot),fill_value='extrapolate'))
    actual_Teffs09.append(teff09)


def log_mdot038(log_XCa,Teff):
    # Build a temporary interpolator at fixed log_XCa each time.
    # Inefficient, but seems robust, and I don't need it to be fast.
    log_mdots_for_this_XCa = []
    for logmdot in fs038:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs038,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)

def log_mdot06(log_XCa,Teff):
    # Build a temporary interpolator at fixed log_XCa each time.
    # Inefficient, but seems robust, and I don't need it to be fast.
    log_mdots_for_this_XCa = []
    for logmdot in fs06:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs06,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)

def log_mdot09(log_XCa,Teff):
    # Build a temporary interpolator at fixed log_XCa each time.
    # Inefficient, but seems robust, and I don't need it to be fast.
    log_mdots_for_this_XCa = []
    for logmdot in fs09:
        log_mdots_for_this_XCa.append(logmdot(log_XCa))
    f = interpolate.interp1d(actual_Teffs09,log_mdots_for_this_XCa,fill_value='extrapolate')
    return f(Teff)


def log_mdot(log_XCa,Teff,logg):
    # Interpolate the other three functions in logg at this Teff
    f = interpolate.interp1d([7.5,8.0,8.5],
                             [log_mdot038(log_XCa,Teff),log_mdot06(log_XCa,Teff),log_mdot09(log_XCa,Teff)],
                             fill_value='extrapolate')
    return f(logg)

