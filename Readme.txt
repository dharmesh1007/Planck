PLANCK PACKAGE

When performing analysis of astrophysical data, one may be required to evaluate the power or the photon
emission rate from an astrophysical source (a star for example). As is usually the case, the source of 
emission can be approximatedto a blackbody source, in which case integration of the Planck function will 
be required.

You may only wish to know the power or the photon emission rate within a specific energy, wavelength or
frequency range. This requires integration between an upper and lower limit. This is not a very easy thing
to do analytically, so numerical techniques are usually involved.

The functions contained within the BB_Planck module use the method of Monte Carlo integration to evaluate 
these integrations and their associated errors. These functions are:
Radiance (Power/(area * wavelength or frequency * solid angle)),
Radiant emittance (Power/(area * wavelength or frequency)),
Photon radiance (Photons/(time * area * wavelength or frequency * solid angle)), and 
Photon radiant emittance (Photons/(time * area * wavelength or frequency))

In addition to these functions is the BB_plot function which plots the radiance as a function of wavelength
or frequency.

All these function require an upper and lower limit of integration in the case of the first four functions
mentioned, and upper and lower limits of the x axis in the case of the BB_plot. These limits must be in astropy
units of wavelength (m, nm or Angstroms), frequency (Hz) or energy (J, eV, or ergs). Therefore BEFORE 
USING THE FUNCTIONS, ONE MUST DOWNLOAD THE ASTROPY PACKAGE AND FROM ASTROPY IMPORT UNITS.

One further required parameter is the temperature of the blackbody source. This should be left unitless as
units of Kelvin are assumed.

A final optional parameter, used in the first four functions, allows one to choose the sampling for the
Monte Carlo integration, which by default is set to 10,000. Lowering this value increase the speed of
calculation but also increases the associated error. Increasing the value increases the time taken for
calculation but decreases the associated error.

The additional modules To_Hz and Monte_Carlo_integration contain functions used in the above calculations.






