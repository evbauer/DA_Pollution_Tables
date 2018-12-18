
[![DOI](https://zenodo.org/badge/162072477.svg)](https://zenodo.org/badge/latestdoi/162072477)

This repository provides tables of diffusion timescales for polluted white dwarfs.
It also provides tables of surface Ca mass fractions as a function of accretion rate for MESA models accreting bulk earth composition including thermohaline mixing.

The files `InterpDiffusion.py` and `InterpThermohaline.py` give examples of interpolating on these tables to map observable parameters to accretion rates. They each provide a `log_mdot` function for this purpose. These functions take three arguments: surface mass fraction of Ca, Teff, and logg. Here's an example of how to use them in a simple python script:

	import InterpThermohaline
	import InterpDiffusion
	
	XCa = 1e-6 # surface mass fraction of Calcium
	Teff = 11000
	logg = 8.1
	
	InterpDiffusion.log_mdot(XCa,Teff,logg)
	InterpThermohaline.log_mdot(XCa,Teff,logg)
	
This will print log10(accretion rate) in units of g/s. The first will assume no thermohaline mixing occurs. The second accounts for thermohaline mixing. In both cases, bulk earth composition is assumed.

You can also use the `log_mdots` functions for lists or arrays to return a numpy array of accretion rates.

	import InterpThermohaline
	import InterpDiffusion
	
	XCa = [1e-7, 1e-6, 3e-6]
	Teff = [11000, 8000, 14000]
	logg = [8.1, 7.8, 8.3]
	
	InterpDiffusion.log_mdots(XCa,Teff,logg)
	InterpThermohaline.log_mdots(XCa,Teff,logg)


The file `plot.py` gives an example of using these interpolation routines to plot accretion rates based on the data from Koester & Wilken (2006, https://www.aanda.org/articles/aa/abs/2006/27/aa4843-06/aa4843-06.html) and Farihi et al. (2012, https://academic.oup.com/mnras/article/424/1/464/1009557).

