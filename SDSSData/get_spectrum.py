from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse

def get_spectrum(datafile):
    '''
    Get the wavelength and flux data of an SDSS spectrum

    datafile = full path to the FITS file you want to get the spectrum for

    Output:

    The loglam and flux columns of the fits file.

    Use:

    loglam, flux = get_spectrum(datafile)
    '''
    hdu_list = fits.open(datafile)
    data = hdu_list[1].data
    hd = hdu_list[1].header
    loglam = data['loglam']
    flux = data['flux']
    return(loglam,flux)
