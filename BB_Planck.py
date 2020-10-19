import numpy as np
import scipy.constants as sc
from astropy import units as u
from matplotlib import pyplot as plt
from To_Hz import conversion_to_Hz
from Monte_Carlo_Integration import MC_int

"""INTEGRATING THE PLANCK FUNCTION---------------------------------------------"""


# Radiance of a blackbody source - Power per area per solid angle (W/m^2/sr)
def Radiance(Temp,lwr,upp,N=100000):
    
    # Run Radiance calculation if all units are correct, else return an error
    # message.
    if (lwr.unit==upp.unit) and (isinstance(Temp, float) or isinstance(Temp, int)):
        # Convert integration limits into frequency(v) units(Hz) before
        # calculating x = hv/kT.
        if (lwr.unit==u.Hz) or (lwr.unit==u.J) or (lwr.unit==u.erg) or (lwr.unit==u.eV):
            upp_a = (sc.h*conversion_to_Hz(upp))/(sc.k*Temp)
            lwr_a = (sc.h*conversion_to_Hz(lwr))/(sc.k*Temp)
        elif (lwr.unit==u.m) or (lwr.unit==u.nm) or (lwr.unit==u.Angstrom):
            upp_a = (sc.h*conversion_to_Hz(lwr))/(sc.k*Temp)
            lwr_a = (sc.h*conversion_to_Hz(upp))/(sc.k*Temp)
        print(lwr_a,upp_a)
        # Function to be integrated for Radiance calculation.
        def fn(x):
            return x**3/(np.e**x-1)
      
        # Calculate radiance in different ways depending on integration limits.
        if (lwr_a==0) and (upp_a==np.inf):
            error = 0
            return (Temp**4*(2*np.pi**4*sc.k**4)/(15*sc.h**3*sc.c**2), error)*\
                    u.W/u.m**2/u.sr
        elif (lwr_a>0) and (upp_a==np.inf):
            totrad = Temp**4*(2*np.pi**4*sc.k**4)/(15*sc.h**3*sc.c**2)
            notinc = ((2*sc.k**4*Temp**4)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,0,lwr_a,N)[0]
            error = ((2*sc.k**4*Temp**4)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,0,lwr_a,N)[1]
            return (totrad-notinc, error)*u.W/u.m**2/u.sr
        else:
            rad = ((2*sc.k**4*Temp**4)/(sc.h**3*sc.c**2))*\
                   MC_int(fn,lwr_a,upp_a,N)[0]
            error = ((2*sc.k**4*Temp**4)/(sc.h**3*sc.c**2))*\
                   MC_int(fn,lwr_a,upp_a,N)[1]
            return (rad, error)*u.W/u.m**2/u.sr
    else:
        print("ERROR!!!: Temp should be unitless (assumed to be in Kelvins)")
        print("Integration limits should be in units of: Frequency (Hz),")
        print("wavelength (m, nm, or Angstrom), or energy (J, eV, or erg)")


# Radiant emittance of a blackbody source - Power per unit area (W/m^2)
def Radiant_emit(Temp,lwr,upp,N=100000):
    if (lwr.unit==upp.unit) and (isinstance(Temp, float) or isinstance(Temp, int)):
        rademit = Radiance(Temp,lwr,upp,N)[0].value*sc.pi
        rademiterr = Radiance(Temp,lwr,upp,N)[1].value*sc.pi
        return (rademit, rademiterr)*u.W/u.m**2
    else:
        print("ERROR!!!: Temp should be unitless (assumed to be in Kelvins)")
        print("Integration limits should be in units of: Frequency (Hz),")
        print("wavelength (m, nm, or Angstrom), or energy (J, eV, or erg)")


# Photon radiance of a blackbody source (photons/s/m^2/sr)
def Photon_radiance(Temp,lwr,upp,N=100000):
    
    # Run Radiant emittance calculation if all units are correct, else return
    # an error message.
    if (lwr.unit==upp.unit) and (isinstance(Temp, float) or isinstance(Temp, int)):
        # Convert integration limits into frequency(v) units(Hz) before
        # calculating x = hv/kT.
        if (lwr.unit==u.Hz) or (lwr.unit==u.J) or (lwr.unit==u.erg) or (lwr.unit==u.eV):
            upp_a = (sc.h*conversion_to_Hz(upp))/(sc.k*Temp)
            lwr_a = (sc.h*conversion_to_Hz(lwr))/(sc.k*Temp)
        elif lwr.unit==u.m or lwr.unit==u.nm or lwr.unit==u.Angstrom:
            upp_a = (sc.h*conversion_to_Hz(lwr))/(sc.k*Temp)
            lwr_a = (sc.h*conversion_to_Hz(upp))/(sc.k*Temp)
        
        print(lwr_a,upp_a)
        # Function to be integrated for Radiant emittance calculation.
        def fn(x):
            return x**2/(np.e**x-1)
        
        # Riemann zeta function
        rz_3 = 1.202056903159594
        
        # Calculate radiant emittance in different ways depending on 
        # integration limits.
        if (lwr_a==0) and (upp_a==np.inf):
            error = 0
            return ((4*rz_3*sc.k**3)/(sc.h**3*sc.c**2)*Temp**3, error)*\
                    u.photon/u.m**2/u.sr
        elif (lwr_a>0) and (upp_a==np.inf):
            totrad = ((4*rz_3*sc.k**3)/(sc.h**3*sc.c**2))*Temp**3
            notinc = ((2*sc.k**3*Temp**3)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,0,lwr_a,N)[0]
            error = ((2*sc.k**3*Temp**3)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,0,lwr_a,N)[1]
            return (totrad-notinc, error)*u.photon/u.m**2/u.sr
        else:
            rad = ((2*sc.k**3*Temp**3)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,lwr_a,upp_a,N)[0]
            error = ((2*sc.k**3*Temp**3)/(sc.h**3*sc.c**2))*\
                      MC_int(fn,lwr_a,upp_a,N)[1]
            return (rad, error)*u.photon/u.m**2/u.sr
    else:
        print("ERROR!!!: Temp should be unitless (assumed to be in Kelvins)")
        print("Integration limits should be in units of: Frequency (Hz),")
        print("wavelength (m, nm, or Angstrom), or energy (J, eV, or erg)")
        

# Photon radiant emittance of a blackbody source - (Photons/s/m^2)
def Photon_radiant_emit(Temp,lwr,upp,N=100000):
    if (lwr.unit==upp.unit) and (isinstance(Temp, float) or isinstance(Temp, int)):
        photrad = Photon_radiance(Temp,lwr,upp,N)[0].value*sc.pi
        photraderr = Photon_radiance(Temp,lwr,upp,N)[1].value*sc.pi
        return (photrad, photraderr)*u.ph/u.m**2
    else:
        print("ERROR!!!: Temp should be unitless (assumed to be in Kelvins)")
        print("Integration limits should be in units of: Frequency (Hz),")
        print("wavelength (m, nm, or Angstrom), or energy (J, eV, or erg)")

"""
print(Photon_radiance(5778.55, 3.7474e14*u.Hz, 4.2827e14*u.Hz))
print(Photon_radiance(5778.55, 700*u.nm, 800*u.nm))
print(Photon_radiance(5778.55, 7000e-10*u.m, 8000e-10*u.m))
print(Photon_radiance(5778.55, 1.5498 *u.eV, 1.7712*u.eV))
print(Photon_radiance(5778.55, 2.4830532e-19 *u.J, 2.8377751e-19*u.J))
print(Photon_radiance(5778.55, 2.4830532e-12 *u.erg, 2.8377751e-12*u.erg))
print("\n")
print(Photon_radiant_emit(5778.55, 3.7474e14*u.Hz, 4.2827e14*u.Hz))
print(Photon_radiant_emit(5778.55, 7000*u.Angstrom, 8000*u.Angstrom))
print(Photon_radiant_emit(5778.55, 700*u.nm, 800*u.nm))
print(Photon_radiant_emit(5778.55, 7000e-10*u.m, 8000e-10*u.m))
print(Photon_radiant_emit(5778.55, 1.5498 *u.eV, 1.7712*u.eV))
print(Photon_radiant_emit(5778.55, 2.4830532e-19 *u.J, 2.8377751e-19*u.J))
print(Photon_radiant_emit(5778.55, 2.4830532e-12 *u.erg, 2.8377751e-12*u.erg))
"""


"""PLANCK FUNCTION----------------------------------------------------------"""

def BB_plot(lwr, upp, T):
    if (lwr.unit==u.Hz or lwr.unit==u.eV or lwr.unit==u.J or lwr.unit==u.erg or\
        lwr.unit==u.nm or lwr.unit==u.m or lwr.unit==u.Angstrom) and\
        (isinstance(T, float) or isinstance(T, int) and lwr.unit==upp.unit):
        spectral = np.linspace(lwr, upp, 1000)
        spectral_1=conversion_to_Hz(spectral)
        num = (2*sc.h*spectral_1**3)/sc.c**2
        den = (np.e**((sc.h*spectral_1)/(sc.k*T)))-1
        if lwr.unit==u.nm:
            pl_fnc = (num/den)*((sc.c*1e9)/spectral.value**2)*u.W/u.m**2/u.nm/u.sr
        elif lwr.unit==u.m:
            pl_fnc = (num/den)*(sc.c/spectral.value**2)*u.W/u.m**2/u.m/u.sr
        elif lwr.unit==u.Angstrom:
            pl_fnc = (num/den)*((sc.c*1e10)/spectral.value**2)*u.W/u.m**2/u.Angstrom/u.sr
        else:
            pl_fnc = (num/den)*u.W/u.m**2/u.Hz/u.sr
        
        plt.figure(figsize=(12,8)) 
        plt.grid(b=True, which='major', axis='both')
        plt.plot(spectral, pl_fnc)
        plt.xlabel(spectral.unit, fontsize=16, labelpad=15)
        plt.ylabel("B(v or lambda,T)    (%s)" %pl_fnc.unit, fontsize=16, labelpad=15)
        
    else:
        print("ERROR!!!: Temp should be unitless (assumed to be in Kelvins)")
        print("Integration limits should be in units of: Frequency (Hz),")
        print("wavelength (m, nm, or Angstrom), or energy (J, eV, or erg)")

#BB_plot(0*u.eV, 10*u.eV, 100000)

