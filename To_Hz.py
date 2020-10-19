import scipy.constants as sc
from astropy import units as u

def conversion_to_Hz(x):
    if x.unit==u.Hz:
        return x.value
    elif x.unit==u.m:
        return sc.c/x.value
    elif x.unit==u.nm:
        return sc.c/(x.value*1e-9)
    elif x.unit==u.Angstrom:
        return sc.c/(x.value*1e-10)
    elif x.unit==u.J:
        return x.value/sc.h
    elif x.unit==u.erg:
        return (x.value*1e-7)/sc.h
    elif x.unit==u.eV:
        return (x.value*sc.eV/sc.h)